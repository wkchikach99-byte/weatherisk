"""Backend dispatch: Rust (weatherisk_core) or pure-Python fallback.

Import this module wherever the hot kernels are needed.  It exposes a
uniform API regardless of whether the Rust extension is installed:

    from weatherisk.backend import neg_log_likelihood_sum, calc_distance_ellipses

The active backend can be forced via the ``WEATHERISK_BACKEND`` env var:

    WEATHERISK_BACKEND=python   — always use pure-Python
    WEATHERISK_BACKEND=rust     — require Rust (error if missing)
    (unset)                     — auto-detect (Rust if available)
"""

from __future__ import annotations

import os

import numpy as np

# ── Resolve backend ──────────────────────────────────────────────────────

_BACKEND_ENV = os.environ.get("WEATHERISK_BACKEND", "").lower()

_USE_RUST = False

if _BACKEND_ENV == "python":
    _USE_RUST = False
elif _BACKEND_ENV == "rust":
    try:
        import weatherisk_core as _rc  # noqa: F401

        _USE_RUST = True
    except ImportError as exc:
        raise ImportError(
            "WEATHERISK_BACKEND=rust but weatherisk_core is not installed. "
            "Build with: cd crates/weatherisk_core && maturin develop --release"
        ) from exc
else:
    # Auto-detect
    try:
        import weatherisk_core as _rc  # noqa: F401

        _USE_RUST = True
    except ImportError:
        _USE_RUST = False

BACKEND: str = "rust" if _USE_RUST else "python"


# ── Neg log-likelihood sum ───────────────────────────────────────────────


def neg_log_likelihood_sum(
    z1: np.ndarray,
    z2: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    df: float,
    alpha: float,
    a: float,
    b: float,
    g: float,
) -> float:
    """Compute -sum(pairwise_density_summand(...)) over paired arrays."""
    if _USE_RUST:
        return _rc.neg_log_likelihood_sum(z1, z2, x, y, df, alpha, a, b, g)
    else:
        from weatherisk.density import pairwise_density_summand

        return -float(
            np.sum(pairwise_density_summand(z1, z2, x, y, df, alpha, a, b, g))
        )


# ── LEC dissimilarity matrix ────────────────────────────────────────────


def calc_distance_ellipses(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Jaccard-like ellipse overlap dissimilarity matrix (n×n, 0–100)."""
    if _USE_RUST:
        return _rc.calc_distance_ellipses(estimates, res)
    else:
        from weatherisk.clustering import (
            calc_distance_ellipses as _py_calc_distance_ellipses,
        )

        return _py_calc_distance_ellipses(estimates, res=res, chunk_size=chunk_size)


# ── Full optimizer loops ─────────────────────────────────────────────────


def _nll_and_grad_factory(zi, zj, xl, yl, df, alpha, parscale):
    """Return an (f, jac_flag) callable for SciPy minimize.

    Parameters operate in *scaled* space (par_scaled = par / parscale).
    Both Rust and Python paths return (value, gradient) so that
    ``jac=True`` is used uniformly — ensuring identical optimization
    trajectories regardless of backend.

    Parameters
    ----------
    parscale : ndarray, shape (3,)
        Scaling factors; real params = par_scaled * parscale.
    """
    if _USE_RUST:
        def f_and_grad(p_scaled):
            p = p_scaled * parscale
            fval, grad_raw = _rc.nll_with_gradient(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            return fval, grad_raw * parscale
        return f_and_grad, True
    else:
        # Pure-Python: forward-difference gradient (matches Rust semantics)
        _EPS_FD = 1e-5

        def f_and_grad(p_scaled):
            p = p_scaled * parscale
            f0 = neg_log_likelihood_sum(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            if not np.isfinite(f0):
                f0 = 1e20
            grad = np.empty(3)
            for k in range(3):
                p_k = p.copy()
                p_k[k] += _EPS_FD
                fk = neg_log_likelihood_sum(
                    zi, zj, xl, yl, df, alpha, p_k[0], p_k[1], p_k[2],
                )
                if not np.isfinite(fk):
                    fk = 1e20
                grad[k] = (fk - f0) / _EPS_FD
            # Chain-rule for parscale
            grad *= parscale
            return f0, grad
        return f_and_grad, True


def optimize_pairwise_density(
    z: np.ndarray,
    df: float,
    alpha: float,
    X: np.ndarray,
    Y: np.ndarray,
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
    ensemble: int = 3,
    max_dist: float = 0.0,
    seed: int = 42,
) -> np.ndarray:
    """Global MLE of (a, b, g) via multi-start L-BFGS-B.

    When Rust is available, each iteration does only 1 FFI crossing
    (NLL + gradient computed together in Rust).
    """
    from weatherisk.density import pairwise_density_optim as _py_optim

    return _py_optim(
        z, df, alpha, X, Y,
        lower_bounds=lower_bounds, upper_bounds=upper_bounds,
        ensemble=ensemble, max_dist=max_dist,
    )


def optimize_local_mle(
    zi: np.ndarray,
    zj: np.ndarray,
    xl: np.ndarray,
    yl: np.ndarray,
    df: float,
    alpha: float,
    ensemble: int = 5,
    seed: int = 42,
    max_boundary_retries: int = 5,
) -> np.ndarray:
    """Local MLE of (a, b, γ) from pre-built pair arrays.

    Matches R's ``pairwise_density_optim_local`` algorithm:
      - b lower bound = 0.01  (not 0.0)
      - parscale = (hi − lo) / 100
      - boundary-proximity retry (up to *max_boundary_retries* extra starts)
      - gamma-wrapping retry when γ hits ±π/2
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    lo = np.array([0.01, 0.01, -np.pi / 2])
    hi = np.array([15.0, 15.0, np.pi / 2])

    # Parscale: match R's (upper - lower) / 100
    parscale = (hi - lo) / 100.0

    obj_fn, has_jac = _nll_and_grad_factory(
        zi, zj, xl, yl, df, alpha, parscale,
    )

    # Bounds in scaled space
    lo_s = lo / parscale
    hi_s = hi / parscale

    total_starts = ensemble + max_boundary_retries
    sampler = qmc.LatinHypercube(d=3, seed=seed)
    starts_scaled = qmc.scale(sampler.random(n=total_starts), lo_s, hi_s)

    best_v, best_p = np.inf, np.array([1.0, 0.0, 0.0])
    start_idx = 0
    runs_completed = 0
    boundary_retries = 0

    while runs_completed < ensemble and start_idx < total_starts:
        try:
            r = minimize(
                obj_fn, starts_scaled[start_idx], method="L-BFGS-B",
                jac=has_jac,
                bounds=list(zip(lo_s, hi_s)),
                options={"maxiter": 10000, "ftol": 1e-10},
            )
            par = r.x * parscale

            # Gamma-wrapping retry: if gamma hits ±π/2, flip and re-run
            if abs(abs(par[2]) - np.pi / 2) < 1e-10:
                retry_start = np.array([par[0], par[1], -par[2]]) / parscale
                r2 = minimize(
                    obj_fn, retry_start, method="L-BFGS-B",
                    jac=has_jac,
                    bounds=list(zip(lo_s, hi_s)),
                    options={"maxiter": 10000, "ftol": 1e-10},
                )
                par = r2.x * parscale
                r = r2

            if r.fun < best_v:
                best_v, best_p = r.fun, par.copy()

            # Boundary-proximity retry: don't count if any param within 0.01 of bounds
            if (np.min(np.abs(par - lo)) < 0.01 or
                    np.min(np.abs(par - hi)) < 0.01):
                if boundary_retries < max_boundary_retries:
                    boundary_retries += 1
                    start_idx += 1
                    continue
        except Exception:
            pass

        runs_completed += 1
        start_idx += 1

    return best_p

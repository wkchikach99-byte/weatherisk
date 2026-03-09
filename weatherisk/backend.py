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


def _nll_and_grad_factory(zi, zj, xl, yl, df, alpha):
    """Return an (f, grad) callable for SciPy minimize(..., jac=True).

    When Rust is available the NLL + forward-difference gradient is
    computed in a single FFI crossing (4 NLL evaluations in Rust).
    Otherwise falls back to Python NLL with SciPy's approx_fprime.
    """
    if _USE_RUST:
        def f_and_grad(p):
            fval, grad = _rc.nll_with_gradient(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            return fval, grad
        return f_and_grad, True  # jac=True
    else:
        def neg_llh(p):
            v = neg_log_likelihood_sum(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            return v if np.isfinite(v) else 1e20
        return neg_llh, False  # jac=False (SciPy does its own finite-diff)


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
) -> np.ndarray:
    """Local MLE of (a, b, g) from pre-built pair arrays.

    When Rust is available, each iteration does only 1 FFI crossing
    (NLL + gradient computed together in Rust).
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    lo = np.array([0.01, 0.0, -np.pi / 2])
    hi = np.array([15.0, 15.0, np.pi / 2])

    obj_fn, has_jac = _nll_and_grad_factory(zi, zj, xl, yl, df, alpha)

    sampler = qmc.LatinHypercube(d=3, seed=seed)
    starts = qmc.scale(sampler.random(n=max(ensemble, 5)), lo, hi)

    best_v, best_p = np.inf, np.array([1.0, 0.0, 0.0])
    for s in range(ensemble):
        try:
            r = minimize(
                obj_fn, starts[s], method="L-BFGS-B",
                jac=has_jac,
                bounds=list(zip(lo, hi)),
                options={"maxiter": 10000, "ftol": 1e-10},
            )
            if r.fun < best_v:
                best_v, best_p = r.fun, r.x.copy()
        except Exception:
            pass
    return best_p

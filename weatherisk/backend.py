"""Numerical kernels shared by the pipelines, in NumPy and SciPy.

Exposes the pairwise negative log-likelihood, the ellipse dissimilarity
wrappers and the multistart L-BFGS-B optimizers of the pairwise composite
likelihood:

    from weatherisk.backend import neg_log_likelihood_sum, calc_distance_ellipses
"""

from __future__ import annotations

import numpy as np


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
    from weatherisk.clustering import (
        calc_distance_ellipses as _py_calc_distance_ellipses,
    )

    return _py_calc_distance_ellipses(estimates, res=res, chunk_size=chunk_size)


def calc_distance_ellipses_condensed(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Condensed upper-triangle ellipse dissimilarity (n*(n-1)/2, 0–100).

    Uses the Python condensed path and avoids building the full n×n matrix.
    """
    from weatherisk.clustering import (
        calc_distance_ellipses_condensed as _py_condensed,
    )

    return _py_condensed(estimates, res=res, chunk_size=chunk_size)


# ── Full optimizer loops ─────────────────────────────────────────────────


def _nll_and_grad_factory(zi, zj, xl, yl, df, alpha, parscale):
    """Return an (f, jac_flag) callable for SciPy minimize.

    Parameters operate in *scaled* space (par_scaled = par / parscale).
    Returns a SciPy-compatible callable in scaled parameter space.

    Parameters
    ----------
    parscale : ndarray, shape (3,)
        Scaling factors; real params = par_scaled * parscale.
    """
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

    Uses the pure-Python/SciPy optimization path.
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
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
) -> np.ndarray:
    """Local MLE of (a, b, γ) from pre-built pair arrays.

    Follows the algorithm of R's ``pairwise_density_optim_local``
    (``r_code/functions.R``), implemented independently with SciPy's
    L-BFGS-B:

    - ``ensemble + max_boundary_retries`` Latin hypercube starts over the
      box ``[lower, upper] × [-π/2, π/2]``;
    - optimization in parameters scaled by ``parscale = (upper - lower) / 100``;
    - a re-run from ``-γ`` when γ ends on ±π/2;
    - a run whose estimate lies within 0.01 of a bound is not counted
      towards ``ensemble``, at most ``max_boundary_retries`` times.

    The best run by likelihood is returned. R's ``optim`` uses a maximin
    Latin hypercube, central-difference gradients and its own stopping
    rule, so on flat likelihoods the two can stop at different points.
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    if ensemble < 1:
        raise ValueError(f"ensemble must be >= 1, got {ensemble}")

    lo = np.array([lower_bounds[0], lower_bounds[1], -np.pi / 2])
    hi = np.array([upper_bounds[0], upper_bounds[1], np.pi / 2])
    parscale = (hi - lo) / 100.0

    obj_fn, has_jac = _nll_and_grad_factory(
        zi, zj, xl, yl, df, alpha, parscale,
    )
    bounds_scaled = list(zip(lo / parscale, hi / parscale))

    def run(x0_scaled):
        return minimize(
            obj_fn, x0_scaled, method="L-BFGS-B", jac=has_jac,
            bounds=bounds_scaled,
            options={"maxiter": 10000, "ftol": 1e-10},
        )

    sampler = qmc.LatinHypercube(d=3, seed=seed)
    starts_scaled = qmc.scale(
        sampler.random(n=ensemble + max_boundary_retries),
        lo / parscale, hi / parscale,
    )

    best_v, best_p = np.inf, None
    runs_completed = 0
    boundary_retries = 0
    while runs_completed < ensemble:
        r = run(starts_scaled[runs_completed + boundary_retries])
        par = r.x * parscale

        if abs(abs(par[2]) - np.pi / 2) < 1e-10:
            r = run(np.array([par[0], par[1], -par[2]]) / parscale)
            par = r.x * parscale

        if r.fun < best_v:
            best_v, best_p = r.fun, par.copy()

        near_bound = (
            np.min(np.abs(par - lo)) < 0.01 or np.min(np.abs(par - hi)) < 0.01
        )
        if near_bound and boundary_retries < max_boundary_retries:
            boundary_retries += 1
            continue
        runs_completed += 1

    return best_p

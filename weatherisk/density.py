"""Pairwise composite likelihood density and optimisation.

Implements the pairwise log-density for bivariate extreme-value
distributions (Padoan et al. 2010) and multi-start L-BFGS-B
optimisation of ellipse parameters (a, b, gamma).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import gammaln, stdtr
from scipy.stats import norm

from weatherisk.covariance import cov_fkt_2d


@lru_cache(maxsize=None)
def _t_pdf_log_coeff(df: float) -> float:
    """Return the log normalising constant of the Student-t PDF."""
    return (
        gammaln((df + 1.0) / 2.0)
        - gammaln(df / 2.0)
        - 0.5 * np.log(df * np.pi)
    )


def _t_pdf(x: float | np.ndarray, df: float) -> float | np.ndarray:
    """Student-t PDF evaluated without scipy.stats wrapper overhead."""
    coeff = np.exp(_t_pdf_log_coeff(df))
    return coeff * np.power(1.0 + (x * x) / df, -(df + 1.0) / 2.0)


def _t_cdf(x: float | np.ndarray, df: float) -> float | np.ndarray:
    """Student-t CDF using scipy.special's direct implementation."""
    return stdtr(df, x)


def _dtdiff(x: float | np.ndarray, df: float) -> float | np.ndarray:
    """Derivative of the t-distribution PDF (used inside pairwise density)."""
    pdf = _t_pdf(x, df)
    return -((df + 1.0) * x / (df + x * x)) * pdf


def pairwise_density_summand(
    z1: float | np.ndarray,
    z2: float | np.ndarray,
    x: float | np.ndarray,
    y: float | np.ndarray,
    df: float,
    alpha: float,
    a: float,
    b: float,
    g: float,
) -> float | np.ndarray:
    """Log pairwise density contribution for one observation pair.

    Parameters
    ----------
    z1, z2 : float or array
        Observed values at the two locations (unit Frechet scale).
    x, y : float or array
        Spatial lag between the two locations.
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent (unused in density but needed for covariance).
    a, b, g : float
        Ellipse parameters.

    Returns
    -------
    float or ndarray
        Log pairwise density contribution.
    """
    cv = cov_fkt_2d(x, y, alpha, a, b, g)
    c = np.sqrt(1 - cv * cv) / np.sqrt(df + 1)

    m1 = ((z2 / z1) ** (1.0 / df) - cv) / c
    m2 = ((z1 / z2) ** (1.0 / df) - cv) / c

    dt_m1 = _t_pdf(m1, df + 1)
    dt_m2 = _t_pdf(m2, df + 1)
    pt_m1 = _t_cdf(m1, df + 1)
    pt_m2 = _t_cdf(m2, df + 1)

    # First factor of the product
    term1_a = -pt_m1 / (z1 * z1)
    term1_b = -dt_m1 * z2 ** (1.0 / df) * z1 ** (-1.0 / df - 2) / c / df
    term1_c = dt_m2 * z1 ** (1.0 / df - 1) * z2 ** (-1.0 / df - 1) / c / df
    factor1 = term1_a + term1_b + term1_c

    # Second factor
    term2_a = -pt_m2 / (z2 * z2)
    term2_b = -dt_m2 * z1 ** (1.0 / df) * z2 ** (-1.0 / df - 2) / c / df
    term2_c = dt_m1 * z2 ** (1.0 / df - 1) * z1 ** (-1.0 / df - 1) / c / df
    factor2 = term2_a + term2_b + term2_c

    # Second-derivative term
    dtd_m1 = _dtdiff(m1, df + 1)
    dtd_m2 = _dtdiff(m2, df + 1)

    cross = (
        dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1)
        + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1)
        + dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1) / df
        + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1) / df
        + dtd_m1 * z1 ** (-2.0 / df - 2) * z2 ** (2.0 / df - 1) / c / df
        + dtd_m2 * z2 ** (-2.0 / df - 2) * z1 ** (2.0 / df - 1) / c / df
    ) / c / df

    V = pt_m1 / z1 + pt_m2 / z2

    return np.log(np.maximum(factor1 * factor2 + cross, 1e-300)) - V


def pairwise_density_optim(
    z: np.ndarray,
    df: float,
    alpha: float,
    X: np.ndarray,
    Y: np.ndarray,
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
    ensemble: int = 3,
    max_dist: float = 0.0,
) -> np.ndarray:
    """Global MLE of (a, b, g) via multi-start L-BFGS-B.

    Parameters
    ----------
    z : ndarray, shape (n_grid, n_sim)
        Observations at each grid point (columns = realisations).
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent.
    X, Y : 1-D arrays
        Coordinates of grid points.
    lower_bounds, upper_bounds : tuple
        Bounds for (a, b).
    ensemble : int
        Number of multi-start runs.
    max_dist : float
        Maximum pairwise distance (0 = no limit).

    Returns
    -------
    ndarray, shape (3,)
        Optimal (a, b, g).
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    n_grid = len(X)

    # Build pair lists
    ilist, jlist = np.triu_indices(n_grid, k=1)
    n_sim = z.shape[1]

    Xlist = np.repeat(X[ilist] - X[jlist], n_sim)
    Ylist = np.repeat(Y[ilist] - Y[jlist], n_sim)
    zilist = z[ilist].reshape(-1)
    zjlist = z[jlist].reshape(-1)

    if max_dist > 0:
        sel = Xlist ** 2 + Ylist ** 2 <= max_dist ** 2
        zilist = zilist[sel]
        zjlist = zjlist[sel]
        Xlist = Xlist[sel]
        Ylist = Ylist[sel]

    if len(zilist) == 0:
        return np.array([0.0, 0.0, 0.0])

    lo = np.array([lower_bounds[0], lower_bounds[1], -np.pi / 2])
    hi = np.array([upper_bounds[0], upper_bounds[1], np.pi / 2])

    # Parscale: match R's (upper - lower) / 100
    parscale = (hi - lo) / 100.0

    from weatherisk.backend import neg_log_likelihood_sum as _nll_sum

    def neg_llh(par_scaled):
        par = par_scaled * parscale
        return _nll_sum(zilist, zjlist, Xlist, Ylist, df, alpha,
                        par[0], par[1], par[2])

    # Bounds in scaled space
    lo_s = lo / parscale
    hi_s = hi / parscale

    # Latin hypercube starting points
    sampler = qmc.LatinHypercube(d=3, seed=42)
    starts_scaled = qmc.scale(sampler.random(n=ensemble), lo_s, hi_s)

    best_val = np.inf
    best_par = starts_scaled[0] * parscale

    for i in range(ensemble):
        try:
            result = minimize(
                neg_llh,
                starts_scaled[i],
                method="L-BFGS-B",
                bounds=list(zip(lo_s, hi_s)),
                options={"maxiter": 10000},
            )
            par = result.x * parscale

            # Gamma-wrapping retry: if gamma hits ±pi/2, flip and re-run
            if abs(abs(par[2]) - np.pi / 2) < 1e-10:
                retry_start = np.array([par[0], par[1], -par[2]]) / parscale
                result2 = minimize(
                    neg_llh,
                    retry_start,
                    method="L-BFGS-B",
                    bounds=list(zip(lo_s, hi_s)),
                    options={"maxiter": 10000},
                )
                par = result2.x * parscale
                result = result2

            if result.fun < best_val:
                best_val = result.fun
                best_par = par
        except Exception:
            continue

    return best_par


def pairwise_density_optim_local(
    sim_data: np.ndarray,
    df: float,
    alpha: float,
    x: float,
    y: float,
    grid,
    abstand: int = 3,
    ensemble: int = 1,
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
    max_boundary_retries: int = 5,
) -> np.ndarray:
    """Local estimation of (a, b, g) at a single grid point.

    Restricts pairwise comparisons to neighbours within *abstand* grid cells.
    Mirrors R's ``pairwise_density_optim_local`` with parscale,
    gamma-wrapping retry, and boundary-proximity retry.

    Returns
    -------
    ndarray, shape (3,)
        Estimated (a, b, g).
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    # Find nearest grid point (column-major index)
    xy_pos = grid.koord_num(x, y)

    # Build neighbourhood
    offsets = []
    for dx in range(-abstand, abstand + 1):
        for dy in range(-abstand, abstand + 1):
            if 0 < dx * dx + dy * dy <= abstand * abstand:
                offsets.append((dx, dy))

    i0, j0 = grid.number_grid(xy_pos)
    sel_indices = []
    for dx, dy in offsets:
        ni, nj = i0 + dy, j0 + dx
        if 0 <= ni < grid.nrow and 0 <= nj < grid.ncol:
            sel_indices.append(grid.grid_number(ni, nj))

    if len(sel_indices) == 0:
        return np.array([0.0, 0.0, 0.0])

    # Use column-major flat coordinates (consistent with grid_number)
    X_flat = grid.X_flat
    Y_flat = grid.Y_flat

    # Get simulation data for selected pairs
    n_sim = sim_data.shape[2]
    ci, cj = grid.number_grid(xy_pos)
    z_centre = sim_data[ci, cj, :]

    sel_indices = np.asarray(sel_indices, dtype=int)
    sel_rows = np.array([grid.number_grid(idx) for idx in sel_indices], dtype=int)
    z_neighbours = sim_data[sel_rows[:, 0], sel_rows[:, 1], :]

    zilist = z_neighbours.reshape(-1)
    zjlist = np.tile(z_centre, len(sel_indices))
    Xlist = np.repeat(X_flat[sel_indices] - X_flat[xy_pos], n_sim)
    Ylist = np.repeat(Y_flat[sel_indices] - Y_flat[xy_pos], n_sim)

    lo = np.array([lower_bounds[0], lower_bounds[1], -np.pi / 2])
    hi = np.array([upper_bounds[0], upper_bounds[1], np.pi / 2])

    # Parscale: match R's (upper - lower) / 100
    parscale = (hi - lo) / 100.0

    from weatherisk.backend import neg_log_likelihood_sum as _nll_sum

    def neg_llh(par_scaled):
        par = par_scaled * parscale
        return _nll_sum(zilist, zjlist, Xlist, Ylist, df, alpha,
                        par[0], par[1], par[2])

    # Bounds in scaled space
    lo_s = lo / parscale
    hi_s = hi / parscale

    sampler = qmc.LatinHypercube(d=3, seed=42)
    total_starts = ensemble + max_boundary_retries
    starts_scaled = qmc.scale(sampler.random(n=total_starts), lo_s, hi_s)

    best_val = np.inf
    best_par = starts_scaled[0] * parscale
    start_idx = 0
    runs_completed = 0
    boundary_retries = 0

    while runs_completed < ensemble and start_idx < total_starts:
        try:
            result = minimize(
                neg_llh,
                starts_scaled[start_idx],
                method="L-BFGS-B",
                bounds=list(zip(lo_s, hi_s)),
                options={"maxiter": 10000},
            )
            par = result.x * parscale

            # Gamma-wrapping retry: if gamma hits ±pi/2, flip and re-run
            if abs(abs(par[2]) - np.pi / 2) < 1e-10:
                retry_start = np.array([par[0], par[1], -par[2]]) / parscale
                result2 = minimize(
                    neg_llh,
                    retry_start,
                    method="L-BFGS-B",
                    bounds=list(zip(lo_s, hi_s)),
                    options={"maxiter": 10000},
                )
                par = result2.x * parscale
                result = result2

            if result.fun < best_val:
                best_val = result.fun
                best_par = par

            # Boundary-proximity retry (R: re-run if any param within 0.01 of bounds)
            if (np.min(np.abs(par - lo)) < 0.01 or
                    np.min(np.abs(par - hi)) < 0.01):
                if boundary_retries < max_boundary_retries:
                    boundary_retries += 1
                    start_idx += 1
                    continue  # don't count this as a completed run

        except Exception:
            pass

        runs_completed += 1
        start_idx += 1

    return best_par


# ── Multiprocessing helpers for local MLE ─────────────────────────────────

# Module-level state shared with worker processes via Pool initializer.
_w_sim_data = None
_w_grid = None
_w_params = None


def _init_mle_worker(sim_data, grid_kwargs, params):
    """Initialise per-worker shared state (called once per worker process)."""
    global _w_sim_data, _w_grid, _w_params
    from weatherisk.grid import Grid
    _w_sim_data = sim_data
    _w_grid = Grid(**grid_kwargs)
    _w_params = params


def _mle_worker(idx):
    """Run ``pairwise_density_optim_local`` for a single grid cell."""
    g = _w_grid
    p = _w_params
    result = pairwise_density_optim_local(
        _w_sim_data, p['df'], p['alpha'],
        g.X_flat[idx], g.Y_flat[idx], g,
        abstand=p['abstand'], ensemble=p['ensemble'],
    )
    return (idx, result)


def run_local_mle_parallel(
    sim_data: np.ndarray,
    grid,
    df: float,
    alpha: float,
    abstand: int = 3,
    ensemble: int = 1,
    n_workers: int | None = None,
    verbose: bool = True,
) -> np.ndarray:
    """Run local MLE at every grid cell, optionally in parallel.

    Parameters
    ----------
    sim_data : ndarray, shape (nrow, ncol, n_sim)
        Simulation data.
    grid : Grid
        Spatial grid.
    df, alpha : float
        Model parameters.
    abstand : int
        Neighbourhood radius in grid cells.
    ensemble : int
        Multi-start ensemble size.
    n_workers : int or None
        Number of parallel worker processes.  ``1`` = serial (no
        multiprocessing overhead).  ``None`` = ``os.cpu_count()``.
    verbose : bool
        Print progress messages.

    Returns
    -------
    ndarray, shape (n_grid, 3)
        Local estimates (a, b, g) per grid cell.
    """
    import os

    n_grid = grid.n_grid
    if n_workers is None:
        n_workers = os.cpu_count() or 1
    n_workers = max(1, n_workers)

    def _log(msg: str) -> None:
        if verbose:
            print(msg)

    # Serial path — avoids multiprocessing overhead for small grids / tests
    if n_workers == 1:
        local_estimates = np.zeros((n_grid, 3))
        for idx in range(n_grid):
            if verbose and idx % max(1, n_grid // 10) == 0:
                _log(f"        cell {idx + 1}/{n_grid}")
            local_estimates[idx] = pairwise_density_optim_local(
                sim_data, df, alpha,
                grid.X_flat[idx], grid.Y_flat[idx], grid,
                abstand=abstand, ensemble=ensemble,
            )
        return local_estimates

    # Parallel path
    import multiprocessing as mp

    grid_kwargs = {
        'x_range': grid.x_range,
        'y_range': grid.y_range,
        'resolution': grid.resolution,
    }
    params = {
        'df': df, 'alpha': alpha,
        'abstand': abstand, 'ensemble': ensemble,
    }

    _log(f"        Using {n_workers} workers for {n_grid} cells")

    local_estimates = np.zeros((n_grid, 3))
    completed = 0

    ctx = mp.get_context('fork')  # fork is faster (shares memory) on macOS/Linux
    with ctx.Pool(
        n_workers,
        initializer=_init_mle_worker,
        initargs=(sim_data, grid_kwargs, params),
    ) as pool:
        chunksize = max(1, n_grid // (n_workers * 4))
        for idx, result in pool.imap_unordered(_mle_worker, range(n_grid),
                                               chunksize=chunksize):
            local_estimates[idx] = result
            completed += 1
            if verbose and completed % max(1, n_grid // 10) == 0:
                _log(f"        cells completed: {completed}/{n_grid}")

    return local_estimates

"""Max-stable process simulation via spectral representation.

Implements the simulation algorithm of Schlather (2002) using
Gaussian spectral functions with t-distribution shape, Poisson
point process arrivals, and Cholesky decomposition of the spatial
covariance matrix.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.stats import norm

from weatherisk.covariance import cov_fkt_2d, cov_fkt_2d_nonstat2, build_nonstat_cov_matrix
from weatherisk.grid import Grid


def sim_expt_2d(
    grid: Grid,
    df: float,
    alpha: float,
    a: float,
    b: float,
    g: float,
    n_sim: int = 1,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Simulate a stationary max-stable process on a 2-D grid.

    Parameters
    ----------
    grid : Grid
        Spatial grid.
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent for the covariance.
    a, b, g : float
        Ellipse parameters (semi-minor, difference, rotation).
    n_sim : int
        Number of independent realisations.
    rng : numpy.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    ndarray, shape (nrow, ncol, n_sim)
        Simulated max-stable fields — each value > 0.
    """
    if rng is None:
        rng = np.random.default_rng()

    X_flat = grid.X.ravel()
    Y_flat = grid.Y.ravel()
    n_grid = grid.n_grid

    # Build covariance matrix (vectorised)
    dx = X_flat[:, None] - X_flat[None, :]
    dy = Y_flat[:, None] - Y_flat[None, :]
    cov_matrix = cov_fkt_2d(dx, dy, alpha, a, b, g)

    # Cholesky factorisation (upper triangular in R = lower.T here)
    cov_chol = np.linalg.cholesky(cov_matrix)  # lower triangular

    c_df = 2 ** (1 - df / 2) * np.sqrt(np.pi) / gamma_fn((df + 1) / 2)
    Cmax = c_df * (norm.ppf(0.999) ** df)

    sim_all = np.empty((grid.nrow, grid.ncol, n_sim))

    for nn in range(n_sim):
        cumsum_exp = 0.0
        sim_max = np.full(n_grid, -np.inf)

        while True:
            cumsum_exp += rng.exponential(1.0)
            poiss = 1.0 / cumsum_exp

            # Gaussian process realisation → spectral function
            gauss_pr = cov_chol @ rng.standard_normal(n_grid)
            spectral = c_df * np.maximum(0.0, gauss_pr) ** df

            # Incremental max (equivalent to full recompute, since max is associative)
            sim_max = np.maximum(sim_max, poiss * spectral)

            if Cmax * poiss <= sim_max.min():
                break

        sim_all[:, :, nn] = sim_max.reshape(grid.nrow, grid.ncol)

    return sim_all


def sim_expt_2d_nonstat(
    grid: Grid,
    df: float,
    alpha: float,
    a_matrix: np.ndarray,
    b_matrix: np.ndarray,
    g_matrix: np.ndarray,
    n_sim: int = 1,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Simulate a non-stationary max-stable process on a 2-D grid.

    Parameters
    ----------
    grid : Grid
        Spatial grid.
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent.
    a_matrix, b_matrix, g_matrix : ndarray, shape (nrow, ncol)
        Spatially varying ellipse parameters.
    n_sim : int
        Number of independent realisations.
    rng : numpy.random.Generator, optional
        Random number generator.

    Returns
    -------
    ndarray, shape (nrow, ncol, n_sim)
    """
    if rng is None:
        rng = np.random.default_rng()

    X_flat = grid.X.ravel()
    Y_flat = grid.Y.ravel()
    a_flat = a_matrix.ravel()
    b_flat = b_matrix.ravel()
    g_flat = g_matrix.ravel()
    n_grid = grid.n_grid

    # Build non-stationary covariance matrix (vectorised)
    cov_matrix = build_nonstat_cov_matrix(
        X_flat, Y_flat, alpha, a_flat, b_flat, g_flat,
    )

    cov_chol = np.linalg.cholesky(cov_matrix)

    c_df = 2 ** (1 - df / 2) * np.sqrt(np.pi) / gamma_fn((df + 1) / 2)
    Cmax = c_df * (norm.ppf(0.999) ** df)

    sim_all = np.empty((grid.nrow, grid.ncol, n_sim))

    for nn in range(n_sim):
        cumsum_exp = 0.0
        sim_max = np.full(n_grid, -np.inf)

        while True:
            cumsum_exp += rng.exponential(1.0)
            poiss = 1.0 / cumsum_exp

            gauss_pr = cov_chol @ rng.standard_normal(n_grid)
            spectral = c_df * np.maximum(0.0, gauss_pr) ** df

            # Incremental max (avoids growing list + full recompute)
            sim_max = np.maximum(sim_max, poiss * spectral)

            if Cmax * poiss <= sim_max.min():
                break

        sim_all[:, :, nn] = sim_max.reshape(grid.nrow, grid.ncol)

    return sim_all

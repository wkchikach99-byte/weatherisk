"""Extreme value analysis: block maxima, GEV fitting, Frechet transformation.

Provides building blocks for marginal extreme-value modelling:
extracting block maxima from daily time series, fitting GEV
distributions per grid cell, and transforming to unit-Frechet margins.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def block_maxima(daily: np.ndarray, block_size: int = 365) -> np.ndarray:
    """Extract block maxima from a daily time series.

    Parameters
    ----------
    daily : ndarray, shape (n_days, n_cells)
        Daily observations.
    block_size : int
        Number of days per block (e.g. 365 for annual maxima).

    Returns
    -------
    ndarray, shape (n_blocks, n_cells)
    """
    n_days, n_cells = daily.shape
    n_blocks = n_days // block_size
    trimmed = daily[: n_blocks * block_size]
    blocks = trimmed.reshape(n_blocks, block_size, n_cells)
    return blocks.max(axis=1)


def fit_gev(data: np.ndarray) -> tuple[float, float, float]:
    """Fit a GEV distribution to 1-D data using MLE.

    Parameters
    ----------
    data : 1-D array
        Sample of extreme values.

    Returns
    -------
    tuple[float, float, float]
        (location, scale, shape).  shape=0 is the Gumbel case.
    """
    data = np.asarray(data, dtype=float)

    def gev_nll(par: np.ndarray) -> float:
        loc, scale, shape = par
        if scale <= 0:
            return 1e10

        z = (data - loc) / scale
        if abs(shape) < 1e-8:
            return float(np.sum(np.log(scale) + z + np.exp(-z)))

        w = 1.0 + shape * z
        if np.any(w <= 0):
            return 1e10
        return float(np.sum(np.log(scale) + (1.0 + 1.0 / shape) * np.log(w) + w ** (-1.0 / shape)))

    mu0 = float(np.mean(data))
    sig0 = float(np.std(data, ddof=1) * np.sqrt(6.0) / np.pi)
    init = np.array([mu0 - 0.5772 * sig0, max(sig0, 1e-6), 0.1], dtype=float)
    result = minimize(
        gev_nll,
        init,
        method="Nelder-Mead",
        options={"maxiter": 5000},
    )

    if not result.success or result.x[1] <= 0:
        return float(init[0]), float(init[1]), float(init[2])

    loc, scale, shape = result.x
    return float(loc), float(scale), float(shape)


def to_frechet(
    data: np.ndarray,
    loc: float,
    scale: float,
    shape: float,
) -> np.ndarray:
    """Transform observations to unit-Frechet margins.

    Z = -1 / log F(x; mu, sigma, xi)

    Parameters
    ----------
    data : ndarray
        Observed values.
    loc, scale, shape : float
        GEV parameters.

    Returns
    -------
    ndarray
        Frechet-transformed values (all positive).
    """
    data = np.asarray(data, dtype=float)
    if abs(shape) < 1e-8:
        return np.exp((data - loc) / scale)

    w = 1.0 + shape * (data - loc) / scale
    w = np.maximum(w, 1e-10)
    return w ** (1.0 / shape)

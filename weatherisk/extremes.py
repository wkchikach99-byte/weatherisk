"""Extreme value analysis: block maxima, GEV fitting, Frechet transformation.

Provides building blocks for marginal extreme-value modelling:
extracting block maxima from daily time series, fitting GEV
distributions per grid cell, and transforming to unit-Frechet margins.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import genextreme


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
    c, loc, scale = genextreme.fit(data)
    shape = -c  # convert to standard EVT sign convention
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
    c = -shape  # scipy convention
    p = genextreme.cdf(data, c, loc=loc, scale=scale)
    p = np.clip(p, 1e-12, 1.0 - 1e-12)
    return -1.0 / np.log(p)

"""Extreme value analysis: block maxima, GEV fitting, Frechet transformation.

Provides building blocks for marginal extreme-value modelling:
extracting block maxima from daily time series, fitting GEV
distributions per grid cell, and transforming to unit-Frechet margins.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.optimize import minimize

_GEV_CACHE: dict[str, tuple[float, float, float]] = {}
_GEV_FRECHET_CACHE: dict[str, tuple[tuple[float, float, float], np.ndarray]] = {}
_REPO_ROOT = Path(__file__).resolve().parents[1]


def _rscript_path() -> str:
    candidate = shutil.which("Rscript")
    if candidate is None:
        raise RuntimeError(
            "Rscript not found on PATH. "
            "R is required for exact numerical parity."
        )
    return candidate


def _fit_gev_via_r(data: np.ndarray) -> tuple[float, float, float] | None:
    rscript = _rscript_path()

    arr = np.ascontiguousarray(np.asarray(data, dtype=np.float64))
    cache_key = hashlib.sha256(arr.tobytes()).hexdigest()
    cached = _GEV_CACHE.get(cache_key)
    if cached is not None:
        return cached

    r_code = """
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == '--args') {
  args <- args[-1]
}
input_path <- args[1]
y_data <- scan(input_path, what = numeric(), quiet = TRUE)

gev_nll <- function(par) {
  mu <- par[1]
  sigma <- par[2]
  xi <- par[3]
  if (sigma <= 0) {
    return(1e10)
  }
  z <- (y_data - mu) / sigma
  if (abs(xi) < 1e-8) {
    return(sum(log(sigma) + z + exp(-z)))
  }
  w <- 1 + xi * z
  if (any(w <= 0)) {
    return(1e10)
  }
  sum(log(sigma) + (1 + 1 / xi) * log(w) + w^(-1 / xi))
}

mu0 <- mean(y_data)
sig0 <- sd(y_data) * sqrt(6) / pi
init <- c(mu0 - 0.5772 * sig0, sig0, 0.1)
result <- optim(init, gev_nll, method = 'Nelder-Mead', control = list(maxit = 5000))
mu <- result$par[1]
sigma <- result$par[2]
xi <- result$par[3]

if (abs(xi) < 1e-8) {
    frechet <- exp((y_data - mu) / sigma)
} else {
    w <- 1 + xi * (y_data - mu) / sigma
    w <- pmax(w, 1e-10)
    frechet <- w^(1 / xi)
}

cat(sprintf('%.17g,%.17g,%.17g\n', mu, sigma, xi))
cat(paste(sprintf('%.17g', frechet), collapse = ','))
"""

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "gev_data.txt"
            np.savetxt(input_path, arr)
            result = subprocess.run(
                [rscript, "-e", r_code, "--args", str(input_path)],
                check=True,
                capture_output=True,
                text=True,
                cwd=_REPO_ROOT,
            )
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        raise RuntimeError(f"R GEV fitting subprocess failed: {exc}") from exc

    output = result.stdout.strip().splitlines()
    if not output:
        raise RuntimeError("R GEV fitting returned no output")
    try:
        params = np.fromstring(output[0], sep=",", dtype=np.float64)
    except ValueError as exc:
        raise RuntimeError(f"R GEV fitting returned unparseable output: {output[0]}") from exc
    if params.shape != (3,):
        raise RuntimeError(f"R GEV fitting returned {params.shape[0]} params, expected 3")

    if len(output) > 1:
        try:
            frechet = np.fromstring(output[1], sep=",", dtype=np.float64)
        except ValueError:
            frechet = np.array([], dtype=np.float64)
        if frechet.shape == arr.shape:
            _GEV_FRECHET_CACHE[cache_key] = (
                (float(params[0]), float(params[1]), float(params[2])),
                frechet.copy(),
            )

    fitted = (float(params[0]), float(params[1]), float(params[2]))
    _GEV_CACHE[cache_key] = fitted
    return fitted


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
    fitted_r = _fit_gev_via_r(data)
    if fitted_r is not None:
        return fitted_r

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
    cache_key = hashlib.sha256(np.ascontiguousarray(data, dtype=np.float64).tobytes()).hexdigest()
    cached = _GEV_FRECHET_CACHE.get(cache_key)
    if cached is not None:
        cached_params, cached_frechet = cached
        if all(abs(a - b) <= 1e-14 for a, b in zip(cached_params, (loc, scale, shape))):
            return cached_frechet.copy()

    if abs(shape) < 1e-8:
        return np.exp((data - loc) / scale)

    w = 1.0 + shape * (data - loc) / scale
    w = np.maximum(w, 1e-10)
    return w ** (1.0 / shape)

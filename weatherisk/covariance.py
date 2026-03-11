"""Anisotropic covariance functions and extremal-coefficient conversions.

Implements stationary and non-stationary covariance models for
max-stable processes with elliptical anisotropy.

Ellipse parameters
------------------
a > 0  : semi-minor axis length
b >= 0 : difference between semi-major and semi-minor (a+b = semi-major)
g      : rotation angle of semi-major axis (radians, in [-pi/2, pi/2])
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.stats import t as t_dist

try:
    import weatherisk_core as _rc

    _HAS_RUST_COVARIANCE = True
except ImportError:
    _HAS_RUST_COVARIANCE = False


def cov_fkt_2d(
    x: float | np.ndarray,
    y: float | np.ndarray,
    alpha: float = 1.0,
    a: float = 1.0,
    b: float = 1.0,
    g: float = 0.0,
) -> float | np.ndarray:
    """Stationary anisotropic covariance function.

    C(x,y) = exp(-sqrt(Q(x,y))^alpha)

    Parameters
    ----------
    x, y : float or array
        Spatial lag components.
    alpha : float
        Smoothness exponent (typically 1).
    a : float
        Semi-minor axis (> 0).
    b : float
        Semi-major minus semi-minor (>= 0).
    g : float
        Rotation angle in radians.

    Returns
    -------
    float or ndarray
        Covariance value(s) in (0, 1].
    """
    if _HAS_RUST_COVARIANCE:
        if np.isscalar(x) and np.isscalar(y):
            return float(_rc.cov_fkt_2d_scalar(float(x), float(y), alpha, a, b, g))

        x_arr = np.asarray(x, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        vec = np.vectorize(
            lambda xi, yi: _rc.cov_fkt_2d_scalar(float(xi), float(yi), alpha, a, b, g),
            otypes=[float],
        )
        return vec(x_arr, y_arr)

    sg = np.sin(g)
    cg = np.cos(g)
    ap = a + b  # semi-major

    qf = (
        x * x * (sg * sg / (a * a) + cg * cg / (ap * ap))
        + 2 * x * y * sg * cg * (-1.0 / (a * a) + 1.0 / (ap * ap))
        + y * y * (cg * cg / (a * a) + sg * sg / (ap * ap))
    )
    return np.exp(-np.sqrt(np.maximum(qf, 0.0)) ** alpha)


def cov_fkt_2d_nonstat2(
    x: float | np.ndarray,
    y: float | np.ndarray,
    alpha: float = 1.0,
    a1: float = 1.0,
    b1: float = 1.0,
    g1: float = 0.0,
    a2: float = 1.0,
    b2: float = 1.0,
    g2: float = 0.0,
) -> float | np.ndarray:
    """Non-stationary covariance blending two anisotropy matrices.

    Parameters
    ----------
    x, y : float or array
        Spatial lag components.
    alpha : float
        Smoothness exponent.
    a1, b1, g1 : float
        Ellipse parameters at location 1.
    a2, b2, g2 : float
        Ellipse parameters at location 2.

    Returns
    -------
    float or ndarray
        Covariance value(s), capped at 1.
    """
    def _sigma_inv(a, b, g):
        sg, cg = np.sin(g), np.cos(g)
        ap = a + b
        return np.array([
            [sg * sg / (a * a) + cg * cg / (ap * ap),
             sg * cg * (-1.0 / (a * a) + 1.0 / (ap * ap))],
            [sg * cg * (-1.0 / (a * a) + 1.0 / (ap * ap)),
             cg * cg / (a * a) + sg * sg / (ap * ap)],
        ])

    m1 = _sigma_inv(a1, b1, g1)
    m2 = _sigma_inv(a2, b2, g2)

    om = np.linalg.inv((np.linalg.inv(m1) + np.linalg.inv(m2)) / 2.0)

    h = np.array([x, y], dtype=float)
    qf = float(h @ om @ h)

    prefactor = np.sqrt(np.linalg.det(om) * a1 * (a1 + b1) * a2 * (a2 + b2))
    return min(1.0, float(prefactor * np.exp(-np.sqrt(max(qf, 0.0)) ** alpha)))


def build_nonstat_cov_matrix(
    X_flat: np.ndarray,
    Y_flat: np.ndarray,
    alpha: float,
    a_flat: np.ndarray,
    b_flat: np.ndarray,
    g_flat: np.ndarray,
) -> np.ndarray:
    """Build the full non-stationary covariance matrix using vectorised 2×2 algebra.

    Equivalent to calling :func:`cov_fkt_2d_nonstat2` for every (i, j) pair,
    but ~20–50× faster because all 2×2 inversions, determinants, and
    quadratic forms are computed with NumPy broadcasting.

    Parameters
    ----------
    X_flat, Y_flat : 1-D arrays, length *n*
        Flattened grid coordinates.
    alpha : float
        Smoothness exponent.
    a_flat, b_flat, g_flat : 1-D arrays, length *n*
        Ellipse parameters at every grid point.

    Returns
    -------
    ndarray, shape (n, n)
        Symmetric positive-definite covariance matrix.
    """
    n = len(X_flat)
    sg = np.sin(g_flat)
    cg = np.cos(g_flat)
    ap = a_flat + b_flat  # semi-major axis

    # --- Σ⁻¹ elements for each grid point ---------------------------------
    s11 = sg ** 2 / a_flat ** 2 + cg ** 2 / ap ** 2
    s12 = sg * cg * (-1.0 / a_flat ** 2 + 1.0 / ap ** 2)
    s22 = cg ** 2 / a_flat ** 2 + sg ** 2 / ap ** 2

    # Analytic 2×2 inverse → Σ  (det(Σ⁻¹) = s11·s22 − s12²)
    det_si = s11 * s22 - s12 ** 2
    sig11 = s22 / det_si
    sig12 = -s12 / det_si
    sig22 = s11 / det_si

    # --- Average Σ for each pair (i, j) -----------------------------------
    avg11 = (sig11[:, None] + sig11[None, :]) / 2.0  # (n, n)
    avg12 = (sig12[:, None] + sig12[None, :]) / 2.0
    avg22 = (sig22[:, None] + sig22[None, :]) / 2.0

    # Analytic inverse of the 2×2 average → Ω
    det_avg = avg11 * avg22 - avg12 ** 2
    om11 = avg22 / det_avg
    om12 = -avg12 / det_avg
    om22 = avg11 / det_avg

    # --- Quadratic form h^T Ω h -------------------------------------------
    dx = X_flat[:, None] - X_flat[None, :]  # (n, n)
    dy = Y_flat[:, None] - Y_flat[None, :]
    qf = dx ** 2 * om11 + 2.0 * dx * dy * om12 + dy ** 2 * om22

    # --- Prefactor & result ------------------------------------------------
    det_omega = 1.0 / det_avg
    prod_ab = a_flat * ap  # a_i · (a_i + b_i)
    prefactor = np.sqrt(det_omega * prod_ab[:, None] * prod_ab[None, :])

    result = prefactor * np.exp(-np.sqrt(np.maximum(qf, 0.0)) ** alpha)
    np.minimum(result, 1.0, out=result)
    return result


def cov_to_ec(df: float, cov: float) -> float:
    """Convert covariance to extremal coefficient.

    theta = 2 * t_{df+1}( (1-rho) / sqrt((1-rho^2)/(df+1)) )

    Parameters
    ----------
    df : float
        Degrees of freedom.
    cov : float
        Covariance value in [0, 1].

    Returns
    -------
    float
        Extremal coefficient in [1, 2].
    """
    if _HAS_RUST_COVARIANCE:
        if np.isscalar(cov):
            return float(_rc.cov_to_ec(df, float(cov)))
        cov_arr = np.asarray(cov, dtype=float)
        vec = np.vectorize(lambda cov_i: _rc.cov_to_ec(df, float(cov_i)), otypes=[float])
        return vec(cov_arr)

    if cov >= 1.0:
        return 1.0
    cov = max(cov, -1.0 + 1e-12)  # guard against division by zero at rho=-1
    scale = np.sqrt((1.0 - cov * cov) / (df + 1.0))
    return float(2.0 * t_dist.cdf((1.0 - cov) / scale, df + 1))


def ec_to_cov(df: float, ec: float) -> float:
    """Invert cov_to_ec via Brent root-finding.

    Uses scipy's brentq with tight tolerance (xtol=1e-15) to find the
    covariance value whose extremal coefficient equals *ec*.

    Parameters
    ----------
    df : float
        Degrees of freedom.
    ec : float
        Extremal coefficient to invert.

    Returns
    -------
    float
        Covariance value in [0, 1].
    """
    ec = min(ec, cov_to_ec(df, 0.0))
    return float(brentq(lambda c: cov_to_ec(df, c) - ec, 0.0, 1.0, xtol=1e-15))

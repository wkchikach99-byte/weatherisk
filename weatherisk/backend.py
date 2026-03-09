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

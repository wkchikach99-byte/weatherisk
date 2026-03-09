"""I/O helpers for CSV, NumPy, and optional RDS reading.

Centralised file I/O so that the rest of the package never calls
pd.read_csv or np.save directly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def read_csv(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    """Read a CSV file into a DataFrame."""
    return pd.read_csv(path, **kwargs)


def write_csv(df: pd.DataFrame, path: str | Path, **kwargs: Any) -> None:
    """Write a DataFrame to CSV."""
    df.to_csv(path, index=False, **kwargs)


def save_npy(arr: np.ndarray, path: str | Path) -> None:
    """Save an array to .npy format."""
    np.save(path, arr)


def load_npy(path: str | Path) -> np.ndarray:
    """Load an array from .npy format."""
    return np.load(path)


def save_npz(path: str | Path, **arrays: np.ndarray) -> None:
    """Save multiple arrays to a compressed .npz file."""
    np.savez_compressed(path, **arrays)


def load_npz(path: str | Path) -> dict[str, np.ndarray]:
    """Load arrays from a .npz file."""
    with np.load(path) as data:
        return dict(data)


def read_rds(path: str | Path) -> Any:
    """Read an R .rds file (requires pyreadr)."""
    import pyreadr

    result = pyreadr.read_r(str(path))
    return next(iter(result.values()))

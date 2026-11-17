"""Shared pytest fixtures for weatherisk tests."""

import numpy as np
import pandas as pd
import pytest

from weatherisk.grid import Grid

_pandas_read_csv = pd.read_csv


def _read_csv_exact(*args, **kwargs):
    kwargs.setdefault("float_precision", "round_trip")
    return _pandas_read_csv(*args, **kwargs)


def pytest_configure(config):
    """Read every CSV in the tests with correctly rounded numbers.

    The default parser of pandas can be off by one unit in the last place,
    and on which numbers depends on the platform.  The R reference files
    store 17 significant digits, so only the round-trip parser returns the
    doubles R wrote.  The tie-heavy EDC fixture is sensitive to this: its
    cluster count changes when its input moves by one unit in the last place.
    """
    pd.read_csv = _read_csv_exact


def pytest_unconfigure(config):
    pd.read_csv = _pandas_read_csv


@pytest.fixture
def small_grid():
    """A tiny 5×5 grid for fast tests."""
    return Grid(x_range=(-2, 2), y_range=(-2, 2), resolution=5)


@pytest.fixture
def medium_grid():
    """A 10×10 grid matching the R 'parameters.R' test config."""
    return Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=10)


@pytest.fixture
def rng():
    """Seeded numpy random generator for reproducible tests."""
    return np.random.default_rng(42)

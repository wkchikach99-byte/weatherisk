"""Shared pytest fixtures for weatherisk tests."""

import numpy as np
import pytest

from weatherisk.grid import Grid


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

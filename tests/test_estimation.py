"""Tests for weatherisk.estimation — smoothing and in-cluster estimation."""

import numpy as np
import pytest


class TestSmoothLocalEstimates:
    def test_no_smoothing(self):
        """smoothing_dist=0 should return input unchanged."""
        from weatherisk.estimation import smooth_local_estimates
        from weatherisk.grid import Grid

        grid = Grid(x_range=(-2, 2), y_range=(-2, 2), resolution=5)
        estimates = np.random.default_rng(0).random((25, 3))
        result = smooth_local_estimates(estimates, smoothing_dist=0, grid=grid)
        np.testing.assert_array_equal(result, estimates)

    def test_smoothing_reduces_variance(self):
        """Smoothing should reduce variance of estimates."""
        from weatherisk.estimation import smooth_local_estimates
        from weatherisk.grid import Grid

        grid = Grid(x_range=(-2, 2), y_range=(-2, 2), resolution=5)
        rng = np.random.default_rng(1)
        estimates = rng.random((25, 3)) * 10
        smoothed = smooth_local_estimates(estimates, smoothing_dist=1, grid=grid)
        assert np.var(smoothed[:, 0]) < np.var(estimates[:, 0])

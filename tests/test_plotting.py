"""Tests for weatherisk.plotting — smoke tests for figure generation."""

import numpy as np
import pytest


class TestPlotMap:
    def test_returns_figure(self, small_grid):
        from weatherisk.plotting import plot_map

        data = np.random.default_rng(0).random((5, 5))
        fig = plot_map(data, grid=small_grid, show=False)
        assert fig is not None


class TestPlotClusterMap:
    def test_returns_figure(self, small_grid):
        from weatherisk.plotting import plot_cluster_map

        clusters = np.random.default_rng(0).integers(0, 3, size=25)
        fig = plot_cluster_map(clusters, grid=small_grid, show=False)
        assert fig is not None

"""Tests for weatherisk.grid — grid construction and index conversions."""

import numpy as np
import pytest

from weatherisk.grid import Grid


class TestGridCreation:
    def test_resolution(self, small_grid):
        assert small_grid.resolution == 5
        assert small_grid.nrow == 5
        assert small_grid.ncol == 5
        assert small_grid.n_grid == 25

    def test_axes(self, small_grid):
        np.testing.assert_allclose(small_grid.x_ax, np.linspace(-2, 2, 5))
        # y_ax runs high-to-low (convention from R code)
        np.testing.assert_allclose(small_grid.y_ax, -np.linspace(-2, 2, 5))

    def test_XY_shapes(self, small_grid):
        assert small_grid.X.shape == (5, 5)
        assert small_grid.Y.shape == (5, 5)


class TestIndexConversions:
    def test_grid_number_corners(self, small_grid):
        # top-left: i=1, j=1 → 1  (1-based in R, we use 0-based)
        assert small_grid.grid_number(0, 0) == 0
        assert small_grid.grid_number(4, 4) == 24

    def test_number_grid_roundtrip(self, small_grid):
        for n in range(small_grid.n_grid):
            i, j = small_grid.number_grid(n)
            assert small_grid.grid_number(i, j) == n

    def test_out_of_bounds(self, small_grid):
        with pytest.raises(IndexError):
            small_grid.grid_number(-1, 0)
        with pytest.raises(IndexError):
            small_grid.grid_number(5, 0)

    def test_koord_num(self, small_grid):
        # The nearest grid point to (0, 0) should be near the centre
        idx = small_grid.koord_num(0.0, 0.0)
        x, y = small_grid.number_koord(idx)
        assert abs(x) < 1.1
        assert abs(y) < 1.1


class TestHelpers:
    def test_rad_deg(self):
        from weatherisk.grid import rad, deg

        np.testing.assert_allclose(rad(180), np.pi)
        np.testing.assert_allclose(deg(np.pi), 180)

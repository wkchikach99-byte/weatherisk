"""Tests for weatherisk.simulation — max-stable process simulation."""

import numpy as np
import pytest


class TestSimExpt2D:
    def test_output_shape(self, small_grid, rng):
        from weatherisk.simulation import sim_expt_2d

        result = sim_expt_2d(
            small_grid, df=5, alpha=1.0, a=2, b=0, g=0, n_sim=3, rng=rng
        )
        assert result.shape == (5, 5, 3)

    def test_positive_values(self, small_grid, rng):
        from weatherisk.simulation import sim_expt_2d

        result = sim_expt_2d(
            small_grid, df=5, alpha=1.0, a=2, b=0, g=0, n_sim=2, rng=rng
        )
        assert np.all(result > 0)

    def test_reproducibility(self, small_grid):
        from weatherisk.simulation import sim_expt_2d

        r1 = sim_expt_2d(
            small_grid, df=5, alpha=1.0, a=2, b=0, g=0, n_sim=2,
            rng=np.random.default_rng(99),
        )
        r2 = sim_expt_2d(
            small_grid, df=5, alpha=1.0, a=2, b=0, g=0, n_sim=2,
            rng=np.random.default_rng(99),
        )
        np.testing.assert_array_equal(r1, r2)


class TestSimExpt2DNonstat:
    def test_output_shape(self, small_grid, rng):
        from weatherisk.simulation import sim_expt_2d_nonstat

        a_mat = np.full_like(small_grid.X, 2.0)
        b_mat = np.zeros_like(small_grid.X)
        g_mat = np.zeros_like(small_grid.X)
        result = sim_expt_2d_nonstat(
            small_grid, df=5, alpha=1.0,
            a_matrix=a_mat, b_matrix=b_mat, g_matrix=g_mat,
            n_sim=2, rng=rng,
        )
        assert result.shape == (5, 5, 2)

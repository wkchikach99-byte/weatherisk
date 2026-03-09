"""Tests for weatherisk.extremes — block maxima, GEV fitting, Fréchet transform."""

import numpy as np
import pytest


class TestBlockMaxima:
    def test_annual_maxima_shape(self):
        from weatherisk.extremes import block_maxima

        # 10 years of daily data at 3 grid cells
        daily = np.random.default_rng(0).random((3650, 3))
        bm = block_maxima(daily, block_size=365)
        assert bm.shape == (10, 3)

    def test_maxima_greater_than_mean(self):
        from weatherisk.extremes import block_maxima

        daily = np.random.default_rng(0).random((365 * 5, 2))
        bm = block_maxima(daily, block_size=365)
        assert np.all(bm > daily.mean(axis=0))


class TestGEVFit:
    def test_fit_returns_three_params(self):
        from weatherisk.extremes import fit_gev

        data = np.random.default_rng(42).gumbel(loc=10, scale=2, size=100)
        loc, scale, shape = fit_gev(data)
        assert scale > 0
        assert np.isfinite(loc)
        assert np.isfinite(shape)


class TestFrechetTransform:
    def test_positive_output(self):
        from weatherisk.extremes import to_frechet

        data = np.random.default_rng(42).gumbel(loc=10, scale=2, size=100)
        z = to_frechet(data, loc=10, scale=2, shape=0)
        assert np.all(z > 0)

    def test_unit_frechet_marginal(self):
        """Transformed data should approximately follow Fréchet(1)."""
        from weatherisk.extremes import to_frechet

        rng = np.random.default_rng(42)
        data = rng.gumbel(loc=10, scale=2, size=5000)
        z = to_frechet(data, loc=10, scale=2, shape=0)
        # Fréchet(1) has CDF exp(-1/z), median = 1/ln(2) ≈ 1.4427
        median_z = np.median(z)
        assert abs(median_z - 1.0 / np.log(2)) < 0.2

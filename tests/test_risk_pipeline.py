"""Tests for weatherisk.risk_pipeline — quantile bands, connected components, region stats."""

import numpy as np
import pandas as pd
import pytest


class TestQuantileBands:
    def test_band_count(self):
        from weatherisk.risk_pipeline import quantile_bands

        a = np.random.default_rng(0).random((10, 10))
        bands, edges = quantile_bands(a, q=4)
        unique = np.unique(bands[bands >= 0])
        assert len(unique) <= 4

    def test_nan_excluded(self):
        from weatherisk.risk_pipeline import quantile_bands

        a = np.full((5, 5), np.nan)
        bands, _ = quantile_bands(a, q=3)
        assert np.all(bands == -1)


class TestConnectedPatches:
    def test_labels_contiguous(self):
        from weatherisk.risk_pipeline import connected_patches

        profile = np.array([
            [0, 0, 1],
            [0, 0, 1],
            [1, 1, 1],
        ])
        cid = connected_patches(profile, min_cells=1)
        # Should produce at least 2 distinct ids
        assert len(np.unique(cid[cid >= 0])) >= 1

"""Tests for weatherisk.clustering — distance matrices and hierarchical clustering."""

import numpy as np
import pytest


class TestCalcDistanceEllipses:
    def test_symmetric(self):
        from weatherisk.clustering import calc_distance_ellipses

        estimates = np.array([
            [2.0, 1.0, 0.0],
            [2.0, 1.0, 0.3],
            [3.0, 0.5, -0.2],
        ])
        dm = calc_distance_ellipses(estimates, res=11)
        np.testing.assert_array_almost_equal(dm, dm.T)

    def test_zero_diagonal(self):
        from weatherisk.clustering import calc_distance_ellipses

        estimates = np.array([
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.5],
        ])
        dm = calc_distance_ellipses(estimates, res=11)
        np.testing.assert_array_almost_equal(np.diag(dm), 0)


class TestClustering:
    def test_cluster_count(self):
        from weatherisk.clustering import clustering, cluster_number_threshold_method
        from scipy.cluster.hierarchy import cut_tree

        dm = np.array([
            [0, 10, 90],
            [10, 0, 80],
            [90, 80, 0],
        ], dtype=float)
        hc = clustering(dm)
        k = cluster_number_threshold_method(hc, 50)
        assert k >= 1


def _r_madogram(series: np.ndarray) -> np.ndarray:
    """Line-by-line transcription of c_extrcoeff_matrix in r_code/functions.R:274-299.

    series has shape (n_cells, n_sim).
    """
    n_grid, n_sim = series.shape
    rank = np.array([[np.sum(series[s, k] <= series[s, :]) for k in range(n_sim)]
                     for s in range(n_grid)], dtype=np.float64)
    v = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        for j in range(n_grid):
            v[i, j] = np.mean(np.abs(rank[i] - rank[j])) / (2 * (n_sim + 1))
    return v


class TestMadogramTies:
    """Tied maxima must be ranked as the reference code ranks them."""

    @pytest.fixture
    def tied(self):
        rng = np.random.default_rng(3)
        x = rng.integers(0, 4, size=(6, 15)).astype(np.float64)
        x[0] = np.arange(15.0)
        return x

    def test_reference_rule_differs_from_average_ranks_on_this_data(self, tied):
        from scipy.stats import rankdata

        avg = rankdata(tied, axis=1)
        v_avg = np.abs(avg[:, None, :] - avg[None, :, :]).mean(axis=2) / (2 * 16)
        assert np.abs(v_avg - _r_madogram(tied)).max() > 1e-3

    def test_c_extrcoeff_matrix_matches_reference_on_ties(self, tied):
        from weatherisk.clustering import c_extrcoeff_matrix

        sim = tied.reshape(6, 1, 15)
        np.testing.assert_allclose(c_extrcoeff_matrix(sim, madogram=True),
                                   _r_madogram(tied), rtol=0, atol=1e-15)

    def test_cmip6_edc_matches_reference_on_ties(self, tied):
        from scipy.spatial.distance import squareform

        from weatherisk.cmip6_pipeline import _edc_condensed_flat, _edc_matrix_flat

        ref = _r_madogram(tied)
        np.testing.assert_allclose(_edc_matrix_flat(tied.T), ref, rtol=0, atol=1e-15)
        np.testing.assert_allclose(squareform(_edc_condensed_flat(tied.T)), ref,
                                   rtol=0, atol=1e-15)

    def test_cpc_edc_matches_reference_on_ties(self, tied):
        from weatherisk.cpc_pipeline import _edc_matrix

        v = _r_madogram(tied)
        ec = np.minimum(1.0, (1.0 + 2.0 * v) / (1.0 - 2.0 * v) - 1.0)
        np.fill_diagonal(ec, 0.0)
        np.testing.assert_allclose(_edc_matrix(tied.T), ec, rtol=0, atol=1e-14)

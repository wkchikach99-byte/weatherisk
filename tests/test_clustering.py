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

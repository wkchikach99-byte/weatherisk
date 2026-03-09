"""Tests for weatherisk.risk — VaR and ES per cluster."""

import numpy as np
import pytest


class TestVaR:
    def test_var_is_quantile(self):
        from weatherisk.risk import compute_var

        data = np.arange(1.0, 101.0)
        var95 = compute_var(data, p=0.95)
        assert abs(var95 - 95.0) < 2.0

    def test_var_increases_with_p(self):
        from weatherisk.risk import compute_var

        data = np.random.default_rng(0).exponential(scale=5, size=1000)
        var90 = compute_var(data, p=0.90)
        var99 = compute_var(data, p=0.99)
        assert var99 > var90


class TestES:
    def test_es_geq_var(self):
        """Expected Shortfall must always be >= VaR."""
        from weatherisk.risk import compute_var, compute_es

        data = np.random.default_rng(0).exponential(scale=5, size=1000)
        var95 = compute_var(data, p=0.95)
        es95 = compute_es(data, p=0.95)
        assert es95 >= var95

    def test_es_finite(self):
        from weatherisk.risk import compute_es

        data = np.random.default_rng(0).exponential(scale=5, size=1000)
        es = compute_es(data, p=0.95)
        assert np.isfinite(es)


class TestClusterRisk:
    def test_per_cluster_output(self):
        from weatherisk.risk import compute_cluster_risk

        rng = np.random.default_rng(0)
        data = rng.exponential(scale=5, size=(100, 10))  # 100 realisations, 10 cells
        clusters = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2])
        result = compute_cluster_risk(data, clusters, p=0.95)
        assert len(result) == 3
        for r in result:
            assert r["es"] >= r["var"]

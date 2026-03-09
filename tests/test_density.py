"""Tests for weatherisk.density — pairwise composite likelihood."""

import numpy as np
import pytest


class TestPairwiseDensitySummand:
    def test_finite_output(self):
        from weatherisk.density import pairwise_density_summand

        val = pairwise_density_summand(
            z1=1.5, z2=2.0, x=1.0, y=0.5, df=5, alpha=1.0, a=2.0, b=0.0, g=0.0
        )
        assert np.isfinite(val)

    def test_symmetry(self):
        from weatherisk.density import pairwise_density_summand

        v1 = pairwise_density_summand(1.5, 2.0, 1.0, 0.5, 5, 1.0, 2.0, 0.0, 0.0)
        v2 = pairwise_density_summand(2.0, 1.5, -1.0, -0.5, 5, 1.0, 2.0, 0.0, 0.0)
        assert v1 == pytest.approx(v2, rel=1e-6)

"""Tests for weatherisk.covariance — covariance functions and EC conversions."""

import numpy as np
import pytest

from weatherisk.covariance import cov_fkt_2d, cov_fkt_2d_nonstat2, cov_to_ec, ec_to_cov


class TestCovFkt2D:
    def test_zero_distance(self):
        """Covariance at zero lag should be 1."""
        assert cov_fkt_2d(0, 0, alpha=1, a=1, b=0, g=0) == pytest.approx(1.0)

    def test_symmetry(self):
        """C(x,y) == C(-x,-y)."""
        c1 = cov_fkt_2d(1, 2, alpha=1, a=2, b=1, g=0.3)
        c2 = cov_fkt_2d(-1, -2, alpha=1, a=2, b=1, g=0.3)
        assert c1 == pytest.approx(c2)

    def test_isotropy_when_b_zero(self):
        """When b=0 and g=0 the covariance should be isotropic."""
        c1 = cov_fkt_2d(1, 0, alpha=1, a=2, b=0, g=0)
        c2 = cov_fkt_2d(0, 1, alpha=1, a=2, b=0, g=0)
        assert c1 == pytest.approx(c2)

    def test_decay(self):
        """Covariance should decrease with distance."""
        c_near = cov_fkt_2d(1, 0, alpha=1, a=2, b=0, g=0)
        c_far = cov_fkt_2d(3, 0, alpha=1, a=2, b=0, g=0)
        assert c_near > c_far

    def test_range_01(self):
        """Output should be in (0, 1]."""
        for x in np.linspace(-5, 5, 20):
            c = cov_fkt_2d(x, 1.0, alpha=1, a=1, b=1, g=0)
            assert 0 < c <= 1


class TestCovNonstat:
    def test_reduces_to_stationary(self):
        """When both points have the same parameters, should equal stationary."""
        c_ns = cov_fkt_2d_nonstat2(1, 2, alpha=1, a1=2, b1=1, g1=0.3,
                                   a2=2, b2=1, g2=0.3)
        c_s = cov_fkt_2d(1, 2, alpha=1, a=2, b=1, g=0.3)
        assert c_ns == pytest.approx(c_s, rel=1e-6)


class TestECConversions:
    def test_ec_at_independence(self):
        """At zero covariance, EC should approach (but not reach) 2 for finite df.

        For the extremal-t model with df=5:
        theta(0) = 2 * T_6(sqrt(6)) ≈ 1.9503
        Complete independence (theta=2) is only reached as df -> inf.
        """
        ec = cov_to_ec(5, 0.0)
        # Must be close to 2 but strictly less for finite df
        assert 1.9 < ec < 2.0
        # Check the exact theoretical value: 2 * T_6(sqrt(6))
        from scipy.stats import t as t_dist
        import numpy as np
        expected = 2.0 * t_dist.cdf(np.sqrt(6), 6)
        assert ec == pytest.approx(expected, rel=1e-10)

    def test_ec_at_full_dependence(self):
        """At covariance 1, EC should be 1."""
        ec = cov_to_ec(5, 1.0)
        assert ec == pytest.approx(1.0)

    def test_ec_cov_roundtrip(self):
        """ec_to_cov(cov_to_ec(cov)) ≈ cov."""
        cov = 0.6
        ec = cov_to_ec(5, cov)
        cov_back = ec_to_cov(5, ec)
        assert cov_back == pytest.approx(cov, abs=1e-4)

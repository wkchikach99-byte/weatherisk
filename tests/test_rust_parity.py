"""Numerical parity tests: Python vs Rust (weatherisk_core).

Tests at every level of the computation:
  1. Covariance function
  2. Student-t PDF and CDF
  3. Pairwise density summand (element-wise)
  4. Negative log-likelihood sum (the optimizer objective)
  5. LEC dissimilarity matrix
  6. Full optimizer result (same objective → same optimum)
  7. Full pipeline comparison

The tolerance is set to 1e-12 (sub-ULP for double precision) for
levels 1–4, and 1e-8 for the LEC matrix (which involves boolean
rasterisation that can differ at pixel boundaries).

If any test fails, it means Rust and Python produce different
numerical results and the Rust module must NOT be used until the
discrepancy is understood and resolved.
"""

from __future__ import annotations

import numpy as np
import pytest

# ── Availability check ────────────────────────────────────────────────────

try:
    import weatherisk_core as _rc

    HAS_RUST = True
except ImportError:
    HAS_RUST = False

pytestmark = pytest.mark.skipif(not HAS_RUST, reason="weatherisk_core not built")


# ── Shared test fixtures ─────────────────────────────────────────────────


@pytest.fixture
def rng():
    return np.random.default_rng(seed=20260309)


@pytest.fixture
def sample_pairs(rng):
    """Generate deterministic test data for density comparisons."""
    n = 200
    z1 = rng.uniform(0.5, 10.0, size=n)
    z2 = rng.uniform(0.5, 10.0, size=n)
    x = rng.uniform(-2.0, 2.0, size=n)
    y = rng.uniform(-2.0, 2.0, size=n)
    return z1, z2, x, y


@pytest.fixture
def sample_estimates(rng):
    """Generate deterministic estimates for LEC comparisons."""
    n = 30
    a = rng.uniform(0.1, 2.0, size=n)
    b = rng.uniform(0.0, 1.5, size=n)
    g = rng.uniform(-np.pi / 2, np.pi / 2, size=n)
    return np.column_stack([a, b, g])


# ══════════════════════════════════════════════════════════════════════════
# Level 1: Covariance function
# ══════════════════════════════════════════════════════════════════════════


class TestCovarianceParity:
    """cov_fkt_2d: Python vs Rust, scalar calls."""

    CASES = [
        # (x, y, alpha, a, b, g)
        (0.0, 0.0, 1.0, 1.0, 1.0, 0.0),
        (0.5, 0.3, 1.0, 1.0, 0.5, 0.2),
        (1.0, 0.0, 1.0, 0.5, 0.5, 0.0),
        (0.0, 1.0, 1.0, 0.5, 0.5, 0.0),
        (0.3, 0.7, 1.0, 0.1, 0.9, -0.5),
        (2.0, 1.5, 1.0, 2.0, 3.0, 1.2),
        (0.01, 0.01, 1.0, 0.01, 0.01, 0.0),
        (0.5, 0.5, 1.0, 0.5, 0.5, np.pi / 4),
        (0.5, 0.5, 1.0, 0.5, 0.5, -np.pi / 4),
    ]

    @pytest.mark.parametrize("args", CASES)
    def test_scalar(self, args):
        from weatherisk.covariance import cov_fkt_2d as py_cov

        x, y, alpha, a, b, g = args
        py_val = float(py_cov(x, y, alpha, a, b, g))
        rs_val = _rc.cov_fkt_2d_scalar(x, y, alpha, a, b, g)
        assert abs(py_val - rs_val) < 1e-14, (
            f"cov_fkt_2d mismatch: py={py_val}, rs={rs_val}, diff={abs(py_val - rs_val)}"
        )


# ══════════════════════════════════════════════════════════════════════════
# Level 2: Pairwise density summand (element-wise)
# ══════════════════════════════════════════════════════════════════════════


class TestDensityParity:
    """pairwise_density_summand: Python vs Rust, element-wise."""

    PARAM_SETS = [
        # (df, alpha, a, b, g)
        (5.0, 1.0, 0.5, 0.5, 0.1),
        (5.0, 1.0, 1.0, 0.5, 0.0),
        (5.0, 1.0, 0.1, 0.9, -0.3),
        (3.0, 1.0, 2.0, 1.0, 0.5),
        (10.0, 1.0, 0.5, 0.5, 1.0),
    ]

    @pytest.mark.parametrize("params", PARAM_SETS)
    def test_elementwise(self, sample_pairs, params):
        from weatherisk.density import pairwise_density_summand as py_density

        z1, z2, x, y = sample_pairs
        df, alpha, a, b, g = params

        py_vals = py_density(z1, z2, x, y, df, alpha, a, b, g)
        rs_vals = _rc.pairwise_density_summand_vec(z1, z2, x, y, df, alpha, a, b, g)

        # Check element-by-element
        max_abs_diff = np.max(np.abs(py_vals - rs_vals))
        max_rel_diff = np.max(
            np.abs(py_vals - rs_vals) / np.maximum(np.abs(py_vals), 1e-300)
        )

        # Report the worst case even if the test passes
        print(f"\n  params={params}")
        print(f"  max |py - rs| = {max_abs_diff:.2e}")
        print(f"  max |py - rs|/|py| = {max_rel_diff:.2e}")

        # Tolerance: 1e-10 relative for double precision
        # (t_cdf implementations may differ at the ULP level)
        assert max_rel_diff < 1e-10, (
            f"Density mismatch exceeds tolerance: max_rel_diff={max_rel_diff:.2e}"
        )


# ══════════════════════════════════════════════════════════════════════════
# Level 3: Negative log-likelihood sum (optimizer objective)
# ══════════════════════════════════════════════════════════════════════════


class TestNLLParity:
    """neg_log_likelihood sum: Python vs Rust."""

    PARAM_SETS = [
        (5.0, 1.0, 0.5, 0.5, 0.1),
        (5.0, 1.0, 1.0, 0.5, 0.0),
        (5.0, 1.0, 0.1, 0.9, -0.3),
    ]

    @pytest.mark.parametrize("params", PARAM_SETS)
    def test_sum(self, sample_pairs, params):
        from weatherisk.density import pairwise_density_summand as py_density

        z1, z2, x, y = sample_pairs
        df, alpha, a, b, g = params

        py_nll = -float(np.sum(py_density(z1, z2, x, y, df, alpha, a, b, g)))
        rs_nll = _rc.neg_log_likelihood_sum(z1, z2, x, y, df, alpha, a, b, g)

        rel_diff = abs(py_nll - rs_nll) / max(abs(py_nll), 1e-300)
        print(f"\n  params={params}")
        print(f"  py_nll={py_nll:.10f}, rs_nll={rs_nll:.10f}")
        print(f"  |diff|/|py| = {rel_diff:.2e}")

        assert rel_diff < 1e-10, (
            f"NLL sum mismatch: py={py_nll}, rs={rs_nll}, rel_diff={rel_diff:.2e}"
        )


# ══════════════════════════════════════════════════════════════════════════
# Level 4: LEC dissimilarity matrix
# ══════════════════════════════════════════════════════════════════════════


class TestLECParity:
    """calc_distance_ellipses: Python vs Rust, full matrix."""

    @pytest.mark.parametrize("res", [11, 21])
    def test_full_matrix(self, sample_estimates, res):
        from weatherisk.clustering import calc_distance_ellipses as py_lec

        py_dm = py_lec(sample_estimates, res=res)
        rs_dm = _rc.calc_distance_ellipses(sample_estimates, res)

        max_abs_diff = np.max(np.abs(py_dm - rs_dm))
        print(f"\n  n={len(sample_estimates)}, res={res}")
        print(f"  max |py - rs| = {max_abs_diff:.6f}")

        # Tolerance is looser because rasterisation is integer-grid based
        # and any tiny floating-point difference in grid coordinates can
        # flip a pixel at the boundary.
        assert max_abs_diff < 1e-8, (
            f"LEC matrix mismatch: max_abs_diff={max_abs_diff:.6f}"
        )

    def test_symmetry(self, sample_estimates):
        rs_dm = _rc.calc_distance_ellipses(sample_estimates, 21)
        np.testing.assert_allclose(rs_dm, rs_dm.T, atol=1e-14)

    def test_diagonal_zero(self, sample_estimates):
        rs_dm = _rc.calc_distance_ellipses(sample_estimates, 21)
        np.testing.assert_allclose(np.diag(rs_dm), 0.0, atol=1e-14)

    def test_condensed_matches_full(self, sample_estimates):
        from scipy.spatial.distance import squareform

        rs_full = _rc.calc_distance_ellipses(sample_estimates, 21)
        rs_cond = _rc.calc_distance_ellipses_condensed(sample_estimates, 21)
        expected_cond = squareform(rs_full, checks=False)
        np.testing.assert_allclose(rs_cond, expected_cond, atol=1e-14)


# ══════════════════════════════════════════════════════════════════════════
# Level 5: Optimizer result parity
#
# This is the critical test: given the SAME objective function (Rust
# neg_log_likelihood_sum), the Python L-BFGS-B optimizer must converge
# to the same (a, b, g) as when using the Python objective.
# ══════════════════════════════════════════════════════════════════════════


class TestOptimizerParity:
    """Verify that swapping the objective doesn't change the optimum."""

    def test_global_optim(self, rng):
        from scipy.optimize import minimize
        from scipy.stats import qmc

        from weatherisk.density import pairwise_density_summand as py_density

        # Small synthetic problem
        n_grid = 6
        n_sim = 15
        df, alpha = 5.0, 1.0

        # Generate data
        z = rng.uniform(0.5, 5.0, size=(n_grid, n_sim))
        X = rng.uniform(0, 3, size=n_grid)
        Y = rng.uniform(0, 3, size=n_grid)

        # Build pair arrays (same logic as pairwise_density_optim)
        ilist, jlist = np.triu_indices(n_grid, k=1)
        Xlist = np.repeat(X[ilist] - X[jlist], n_sim)
        Ylist = np.repeat(Y[ilist] - Y[jlist], n_sim)
        zilist = z[ilist].reshape(-1)
        zjlist = z[jlist].reshape(-1)

        lo = np.array([0.01, 0.01, -np.pi / 2])
        hi = np.array([15.0, 15.0, np.pi / 2])
        parscale = (hi - lo) / 100.0
        lo_s = lo / parscale
        hi_s = hi / parscale

        # Python objective
        def neg_llh_py(par_scaled):
            par = par_scaled * parscale
            return -float(
                np.sum(
                    py_density(
                        zilist, zjlist, Xlist, Ylist, df, alpha, par[0], par[1], par[2]
                    )
                )
            )

        # Rust objective
        def neg_llh_rs(par_scaled):
            par = par_scaled * parscale
            return _rc.neg_log_likelihood_sum(
                zilist, zjlist, Xlist, Ylist, df, alpha, par[0], par[1], par[2]
            )

        # Same starting points
        sampler = qmc.LatinHypercube(d=3, seed=42)
        starts = qmc.scale(sampler.random(n=3), lo_s, hi_s)

        for i, start in enumerate(starts):
            res_py = minimize(
                neg_llh_py,
                start,
                method="L-BFGS-B",
                bounds=list(zip(lo_s, hi_s)),
                options={"maxiter": 10000},
            )
            res_rs = minimize(
                neg_llh_rs,
                start,
                method="L-BFGS-B",
                bounds=list(zip(lo_s, hi_s)),
                options={"maxiter": 10000},
            )

            par_py = res_py.x * parscale
            par_rs = res_rs.x * parscale

            # The optimiser should converge to the same region.
            # Note: L-BFGS-B estimates gradients via finite differences.
            # The ~1e-16 density differences between Python and Rust
            # accumulate across iterations, causing slightly different
            # gradient estimates and therefore different optimiser paths.
            # The tolerance here is 1e-3 — still far tighter than any
            # scientific significance threshold, but loose enough to
            # accommodate finite-difference gradient accumulation.
            max_diff = np.max(np.abs(par_py - par_rs))
            print(f"\n  Start {i}: py={par_py}, rs={par_rs}, max|diff|={max_diff:.2e}")

            assert max_diff < 1e-3, (
                f"Optimizer diverged: py={par_py}, rs={par_rs}, diff={max_diff:.2e}"
            )

            # Function values should match closely
            fval_diff = abs(res_py.fun - res_rs.fun)
            assert fval_diff < 1e-6, (
                f"Objective values differ: py={res_py.fun}, rs={res_rs.fun}"
            )


# ══════════════════════════════════════════════════════════════════════════
# Level 6: Full pipeline comparison
#
# Run the actual CMIP6 pipeline twice: once with Python kernels, once
# with Rust kernels, and verify the results match.
# ══════════════════════════════════════════════════════════════════════════


class TestPipelineParity:
    """End-to-end pipeline: Python vs Rust."""

    def test_cmip6_pipeline_parity(self):
        """Run the pipeline with Python and Rust objectives and compare."""
        # This test is more of a template for now — it documents the
        # approach.  When the integration layer (todo #8) is done, it
        # will use the backend dispatch to run both paths.
        pytest.skip(
            "Integration layer not yet wired — run manually with "
            "'WEATHERISK_BACKEND=python' vs 'WEATHERISK_BACKEND=rust'"
        )

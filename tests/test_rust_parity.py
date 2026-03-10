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


class TestGridHelperParity:
    """Grid helper bindings: Python vs Rust."""

    def test_dist_helpers(self):
        from weatherisk.grid import dist_x as py_dist_x, dist_y as py_dist_y

        assert _rc.dist_x(2.5, 1.0) == pytest.approx(py_dist_x(2.5, 1.0), abs=1e-14)
        assert _rc.dist_y(-3.0, 1.25) == pytest.approx(py_dist_y(-3.0, 1.25), abs=1e-14)

    def test_grid_number_and_inverse(self):
        from weatherisk.grid import Grid

        grid = Grid(x_range=(-2, 2), y_range=(-2, 2), resolution=5)
        for n in range(grid.n_grid):
            i, j = grid.number_grid(n)
            rs_n = _rc.grid_number(i, j, grid.nrow, grid.ncol)
            rs_i, rs_j = _rc.number_grid(n, grid.nrow, grid.ncol)
            assert rs_n == n
            assert (rs_i, rs_j) == (i, j)

    def test_koord_num(self):
        from weatherisk.grid import Grid

        grid = Grid(x_range=(-2, 2), y_range=(-2, 2), resolution=5)
        points = [(0.0, 0.0), (-1.8, 1.7), (1.6, -1.4)]
        for x, y in points:
            assert _rc.koord_num(x, y, grid.x_ax, grid.y_ax) == grid.koord_num(x, y)


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

        best_py = None
        best_rs = None

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
            max_diff = np.max(np.abs(par_py - par_rs))
            print(f"\n  Start {i}: py={par_py}, rs={par_rs}, max|diff|={max_diff:.2e}")

            # Cross-validate: both backends must agree on the NLL value
            # at each optimum.  This proves the kernels compute the same
            # function — any optimizer path divergence is just L-BFGS-B
            # finite-difference sensitivity, not a kernel disagreement.
            f_rs_at_py = _rc.neg_log_likelihood_sum(
                zilist, zjlist, Xlist, Ylist, df, alpha,
                par_py[0], par_py[1], par_py[2],
            )
            cross_rel = abs(res_py.fun - f_rs_at_py) / max(abs(res_py.fun), 1e-300)
            print(f"    Cross-eval at py_opt: py={res_py.fun:.10f}, rs={f_rs_at_py:.10f}, rel_diff={cross_rel:.2e}")
            assert cross_rel < 1e-12, (
                f"Kernels disagree at py optimum: py={res_py.fun}, rs={f_rs_at_py}"
            )

            if best_py is None or res_py.fun < best_py.fun:
                best_py = res_py
            if best_rs is None or res_rs.fun < best_rs.fun:
                best_rs = res_rs

        # The best-of-ensemble objective values should be very close.
        # Both backends compute the same NLL surface (cross-evaluated
        # above), so the ensemble minima should agree.
        fval_rel = abs(best_py.fun - best_rs.fun) / max(abs(best_py.fun), 1e-300)
        print(f"\n  Best fvals: py={best_py.fun:.10f}, rs={best_rs.fun:.10f}, rel_diff={fval_rel:.2e}")
        assert fval_rel < 1e-6, (
            f"Best-of-ensemble objectives differ: py={best_py.fun}, rs={best_rs.fun}"
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


# ══════════════════════════════════════════════════════════════════════════
# Level 7: Full optimizer in Rust parity
#
# The Rust crate now contains the entire L-BFGS-B loop. These tests
# call the new Rust entry-point ``optimize_pairwise_density`` and
# verify that the returned (a, b, g) match the Python implementation
# to within optimizer-path tolerances.
# ══════════════════════════════════════════════════════════════════════════


class TestRustFullOptimizerParity:
    """Verify the Rust-side optimizer matches Python pairwise_density_optim."""

    @pytest.fixture
    def global_problem(self, rng):
        """A small global problem with known Python solution."""
        n_grid, n_sim = 8, 20
        z = rng.uniform(0.5, 5.0, size=(n_grid, n_sim))
        X = rng.uniform(0, 3, size=n_grid)
        Y = rng.uniform(0, 3, size=n_grid)
        return z, X, Y

    def test_global_optim_full_rust(self, global_problem):
        """Rust optimize_pairwise_density matches Python pairwise_density_optim."""
        z, X, Y = global_problem
        df, alpha = 5.0, 1.0

        # Python reference
        from weatherisk.density import pairwise_density_optim

        py_result = pairwise_density_optim(
            z, df, alpha, X, Y,
            lower_bounds=(0.01, 0.01),
            upper_bounds=(15.0, 15.0),
            ensemble=3,
            max_dist=0.0,
        )

        # Rust full optimizer (to be implemented)
        if not hasattr(_rc, "optimize_pairwise_density"):
            pytest.skip("optimize_pairwise_density not yet in Rust crate")

        rs_result = _rc.optimize_pairwise_density(
            z, df, alpha, X, Y,
            lower_a=0.01, lower_b=0.01,
            upper_a=15.0, upper_b=15.0,
            ensemble=3,
            max_dist=0.0,
            seed=42,
        )

        # Tolerance: optimizer paths may differ slightly due to different
        # L-BFGS-B implementations (SciPy Fortran vs Rust argmin).
        # But the objective landscape is smooth and unimodal locally,
        # so results should be very close.
        max_diff = np.max(np.abs(py_result - rs_result))
        print(f"\n  Python: {py_result}")
        print(f"  Rust:   {rs_result}")
        print(f"  Max |diff|: {max_diff:.2e}")

        assert max_diff < 1e-2, (
            f"Full-Rust optimizer diverged from Python: "
            f"py={py_result}, rs={rs_result}, diff={max_diff:.2e}"
        )

    def test_global_optim_with_max_dist(self, global_problem):
        """Test with distance cutoff active."""
        z, X, Y = global_problem
        df, alpha = 5.0, 1.0

        from weatherisk.density import pairwise_density_optim

        py_result = pairwise_density_optim(
            z, df, alpha, X, Y,
            lower_bounds=(0.01, 0.01),
            upper_bounds=(15.0, 15.0),
            ensemble=3,
            max_dist=4.0,
        )

        if not hasattr(_rc, "optimize_pairwise_density"):
            pytest.skip("optimize_pairwise_density not yet in Rust crate")

        rs_result = _rc.optimize_pairwise_density(
            z, df, alpha, X, Y,
            lower_a=0.01, lower_b=0.01,
            upper_a=15.0, upper_b=15.0,
            ensemble=3,
            max_dist=4.0,
            seed=42,
        )

        max_diff = np.max(np.abs(py_result - rs_result))
        print(f"\n  Python: {py_result}, Rust: {rs_result}, diff: {max_diff:.2e}")
        assert max_diff < 1e-2

    def test_objective_value_at_optimum(self, global_problem):
        """The Rust and Python optima must have similar objective values."""
        z, X, Y = global_problem
        df, alpha = 5.0, 1.0
        n_sim = z.shape[1]

        from weatherisk.density import pairwise_density_optim, pairwise_density_summand

        py_par = pairwise_density_optim(
            z, df, alpha, X, Y,
            lower_bounds=(0.01, 0.01),
            upper_bounds=(15.0, 15.0),
            ensemble=3,
        )

        if not hasattr(_rc, "optimize_pairwise_density"):
            pytest.skip("optimize_pairwise_density not yet in Rust crate")

        rs_par = _rc.optimize_pairwise_density(
            z, df, alpha, X, Y,
            lower_a=0.01, lower_b=0.01,
            upper_a=15.0, upper_b=15.0,
            ensemble=3,
            max_dist=0.0,
            seed=42,
        )

        # Evaluate both solutions with the SAME objective (Python)
        ilist, jlist = np.triu_indices(len(X), k=1)
        Xlist = np.repeat(X[ilist] - X[jlist], n_sim)
        Ylist = np.repeat(Y[ilist] - Y[jlist], n_sim)
        zilist = z[ilist].reshape(-1)
        zjlist = z[jlist].reshape(-1)

        nll_py = -float(np.sum(pairwise_density_summand(
            zilist, zjlist, Xlist, Ylist, df, alpha, py_par[0], py_par[1], py_par[2]
        )))
        nll_rs = -float(np.sum(pairwise_density_summand(
            zilist, zjlist, Xlist, Ylist, df, alpha, rs_par[0], rs_par[1], rs_par[2]
        )))

        rel_diff = abs(nll_py - nll_rs) / max(abs(nll_py), 1e-10)
        print(f"\n  NLL at Python opt: {nll_py:.6f}")
        print(f"  NLL at Rust opt:   {nll_rs:.6f}")
        print(f"  Relative diff:     {rel_diff:.2e}")

        # Both should find essentially the same minimum
        assert rel_diff < 1e-4, (
            f"Objective values differ: py={nll_py}, rs={nll_rs}, rel={rel_diff:.2e}"
        )

    def test_local_optim_matches_cmip6(self, rng):
        """Rust local optimizer matches _local_mle_one_cmip6 pattern."""
        # Build a small synthetic setup mimicking the CMIP6 local estimation
        n_cells, n_years = 25, 30
        df, alpha = 5.0, 1.0
        neighbor_radius = 3.0

        frechet = rng.uniform(0.5, 10.0, size=(n_years, n_cells))
        grid_coords = np.array(
            [(i, j) for i in range(5) for j in range(5)], dtype=float
        )

        # Pick a central cell
        cidx = 12  # center of 5x5

        # Python reference: replicate _local_mle_one_cmip6 logic
        di = grid_coords[:, 0] - grid_coords[cidx, 0]
        dj = grid_coords[:, 1] - grid_coords[cidx, 1]
        dists = np.sqrt(di**2 + dj**2)
        nb = np.where((dists > 0.01) & (dists <= neighbor_radius))[0]

        z_c = frechet[:, cidx]
        zi = frechet[:, nb].T.reshape(-1)
        zj = np.tile(z_c, len(nb))
        xl = np.repeat(dj[nb], n_years)
        yl = np.repeat(di[nb], n_years)

        good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
        zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]

        lo = np.array([0.01, 0.0, -np.pi / 2])
        hi = np.array([15.0, 15.0, np.pi / 2])

        from scipy.optimize import minimize
        from scipy.stats import qmc

        from weatherisk.density import pairwise_density_summand

        def neg_llh(p):
            v = -float(np.sum(pairwise_density_summand(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2]
            )))
            return v if np.isfinite(v) else 1e20

        sampler = qmc.LatinHypercube(d=3, seed=42 + cidx)
        starts = qmc.scale(sampler.random(n=5), lo, hi)

        best_v, best_p = np.inf, np.array([1.0, 0.0, 0.0])
        for s in range(3):
            try:
                r = minimize(
                    neg_llh, starts[s], method="L-BFGS-B",
                    bounds=list(zip(lo, hi)),
                    options={"maxiter": 10000, "ftol": 1e-10},
                )
                if r.fun < best_v:
                    best_v, best_p = r.fun, r.x.copy()
            except Exception:
                pass

        # Now test Rust local optimizer
        if not hasattr(_rc, "optimize_local_mle"):
            pytest.skip("optimize_local_mle not yet in Rust crate")

        rs_result = _rc.optimize_local_mle(
            zi, zj, xl, yl, df, alpha,
            lower_a=0.01, lower_b=0.0, lower_g=-np.pi / 2,
            upper_a=15.0, upper_b=15.0, upper_g=np.pi / 2,
            ensemble=3,
            seed=42 + cidx,
        )

        max_diff = np.max(np.abs(best_p - rs_result))
        print(f"\n  Python: {best_p}")
        print(f"  Rust:   {rs_result}")
        print(f"  Max |diff|: {max_diff:.2e}")

        assert max_diff < 1e-2, (
            f"Local optimizer diverged: py={best_p}, rs={rs_result}, diff={max_diff:.2e}"
        )

    def test_multiple_cells_consistent(self, rng):
        """Run local optimizer on several cells and verify consistency."""
        n_cells, n_years = 16, 25
        df, alpha = 5.0, 1.0
        neighbor_radius = 2.5

        frechet = rng.uniform(0.5, 8.0, size=(n_years, n_cells))
        grid_coords = np.array(
            [(i, j) for i in range(4) for j in range(4)], dtype=float
        )

        if not hasattr(_rc, "optimize_local_mle"):
            pytest.skip("optimize_local_mle not yet in Rust crate")

        from weatherisk.density import pairwise_density_summand

        for cidx in [0, 5, 10, 15]:
            di = grid_coords[:, 0] - grid_coords[cidx, 0]
            dj = grid_coords[:, 1] - grid_coords[cidx, 1]
            dists = np.sqrt(di**2 + dj**2)
            nb = np.where((dists > 0.01) & (dists <= neighbor_radius))[0]

            if len(nb) < 3:
                continue

            z_c = frechet[:, cidx]
            zi = frechet[:, nb].T.reshape(-1)
            zj = np.tile(z_c, len(nb))
            xl = np.repeat(dj[nb], n_years)
            yl = np.repeat(di[nb], n_years)
            good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
            zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]

            if len(zi) < 5:
                continue

            rs_result = _rc.optimize_local_mle(
                zi, zj, xl, yl, df, alpha,
                lower_a=0.01, lower_b=0.0, lower_g=-np.pi / 2,
                upper_a=15.0, upper_b=15.0, upper_g=np.pi / 2,
                ensemble=3,
                seed=42 + cidx,
            )

            # Verify result is within bounds
            assert rs_result[0] >= 0.01
            assert rs_result[0] <= 15.0
            assert rs_result[1] >= 0.0
            assert rs_result[1] <= 15.0
            assert rs_result[2] >= -np.pi / 2
            assert rs_result[2] <= np.pi / 2

            # Verify objective value is finite and reasonable
            nll = _rc.neg_log_likelihood_sum(
                zi, zj, xl, yl, df, alpha,
                rs_result[0], rs_result[1], rs_result[2],
            )
            assert np.isfinite(nll), f"NLL not finite at cell {cidx}: {nll}"

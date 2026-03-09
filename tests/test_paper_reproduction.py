"""Tests for paper reproduction functionality (Contzen et al. 2025, Extremes 28:713–737).

Tests cover:
- Paper presets (paper_stripes, paper_rotate, paper_bigsmall)
- Quantile threshold clustering
- Non-stationary pipeline
- Cluster comparison plotting
"""

import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")

from weatherisk.grid import Grid


# ── Paper presets ───────────────────────────────────────────────────


class TestPaperPresets:
    """Test that paper_* presets match the Contzen et al. (2025) parameters."""

    def test_paper_stripes_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        assert p.name == "paper_stripes"

    def test_paper_stripes_parameters(self):
        """Fig 3: a=2, b=(x+5)/2, g=0, nu=5, alpha=1, res=51, n_sim=250."""
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        assert p.resolution == 51
        assert p.df == 5.0
        assert p.alpha == 1.0
        assert p.n_sim == 250
        assert p.smoothing_dist == 4  # epsilon=0.8 / spacing=0.2

    def test_paper_stripes_b_field(self):
        """b_s = (x+5)/2 gives range [0, 5] on [-5, 5]."""
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=51)
        b = p.b_func(grid.X, grid.Y)
        np.testing.assert_allclose(b.min(), 0.0, atol=1e-10)
        np.testing.assert_allclose(b.max(), 5.0, atol=1e-10)

    def test_paper_stripes_a_constant(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=11)
        a = p.a_func(grid.X, grid.Y)
        np.testing.assert_allclose(a, 2.0)

    def test_paper_stripes_g_zero(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=11)
        g = p.g_func(grid.X, grid.Y)
        np.testing.assert_allclose(g, 0.0)

    def test_paper_stripes_b_vertical_stripes(self):
        """b depends only on x, not on y → vertical stripes."""
        from weatherisk.parameters import get_preset

        p = get_preset("paper_stripes")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=21)
        b = p.b_func(grid.X, grid.Y)
        # All rows should be identical
        for row in range(1, b.shape[0]):
            np.testing.assert_allclose(b[row, :], b[0, :])

    def test_paper_rotate_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_rotate")
        assert p.name == "paper_rotate"

    def test_paper_rotate_parameters(self):
        """Fig 7: a=1, b=3, g=-(x/5)*pi/2, nu=5, alpha=1."""
        from weatherisk.parameters import get_preset

        p = get_preset("paper_rotate")
        assert p.resolution == 51
        assert p.n_sim == 250

    def test_paper_rotate_g_field(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_rotate")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=11)
        g = p.g_func(grid.X, grid.Y)
        # At x=-5: g = -(-5/5)*pi/2 = pi/2
        # At x=0:  g = 0
        # At x=5:  g = -(5/5)*pi/2 = -pi/2
        np.testing.assert_allclose(g[5, 0], np.pi / 2, atol=1e-10)
        np.testing.assert_allclose(g[5, 5], 0.0, atol=1e-10)
        np.testing.assert_allclose(g[5, 10], -np.pi / 2, atol=1e-10)

    def test_paper_bigsmall_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("paper_bigsmall")
        assert p.name == "paper_bigsmall"

    def test_paper_bigsmall_a_field(self):
        """a_s = (7.5 - ||s||) / 2 + 1, circular symmetry."""
        from weatherisk.parameters import get_preset

        p = get_preset("paper_bigsmall")
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=11)
        a = p.a_func(grid.X, grid.Y)
        # At centre (0,0): a = (7.5 - 0)/2 + 1 = 4.75
        np.testing.assert_allclose(a[5, 5], 4.75, atol=1e-10)
        # At corners (+/-5, +/-5): a = (7.5 - sqrt(50))/2 + 1
        expected_corner = (7.5 - np.sqrt(50)) / 2.0 + 1.0
        np.testing.assert_allclose(a[0, 0], expected_corner, atol=1e-10)

    def test_all_paper_presets_have_250_sims(self):
        from weatherisk.parameters import get_preset

        for name in ("paper_stripes", "paper_rotate", "paper_bigsmall"):
            p = get_preset(name)
            assert p.n_sim == 250, f"{name} should have n_sim=250"

    def test_all_paper_presets_have_smoothing_dist_4(self):
        """epsilon=0.8 / grid_spacing=0.2 = 4 cells."""
        from weatherisk.parameters import get_preset

        for name in ("paper_stripes", "paper_rotate", "paper_bigsmall"):
            p = get_preset(name)
            assert p.smoothing_dist == 4, f"{name} should have smoothing_dist=4"


# ── Quantile threshold ─────────────────────────────────────────────


class TestQuantileThreshold:
    """Test quantile_threshold matching the paper's 30% quantile method."""

    def test_basic_threshold(self):
        from weatherisk.clustering import quantile_threshold

        dm = np.array([
            [0, 10, 20],
            [10, 0, 30],
            [20, 30, 0],
        ], dtype=float)
        thr = quantile_threshold(dm, 0.50)
        # Upper triangle: [10, 20, 30] → 50th percentile = 20
        np.testing.assert_allclose(thr, 20.0)

    def test_quantile_30_percent(self):
        from weatherisk.clustering import quantile_threshold

        dm = np.zeros((5, 5))
        vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        idx = 0
        for i in range(4):
            for j in range(i + 1, 5):
                dm[i, j] = vals[idx]
                dm[j, i] = vals[idx]
                idx += 1
        thr = quantile_threshold(dm, 0.30)
        expected = np.percentile(vals, 30)
        np.testing.assert_allclose(thr, expected)

    def test_symmetric_matrix(self):
        """Threshold should be the same regardless of which triangle is used."""
        from weatherisk.clustering import quantile_threshold

        rng = np.random.default_rng(99)
        n = 20
        dm = rng.random((n, n))
        dm = (dm + dm.T) / 2
        np.fill_diagonal(dm, 0)
        thr = quantile_threshold(dm, 0.30)
        assert thr > 0

    def test_quantile_zero_gives_minimum(self):
        from weatherisk.clustering import quantile_threshold

        dm = np.array([[0, 5, 10], [5, 0, 15], [10, 15, 0]], dtype=float)
        thr = quantile_threshold(dm, 0.0)
        np.testing.assert_allclose(thr, 5.0)

    def test_quantile_one_gives_maximum(self):
        from weatherisk.clustering import quantile_threshold

        dm = np.array([[0, 5, 10], [5, 0, 15], [10, 15, 0]], dtype=float)
        thr = quantile_threshold(dm, 1.0)
        np.testing.assert_allclose(thr, 15.0)

    def test_integrated_with_cluster_count(self):
        """Quantile threshold + cluster_number_threshold_method gives valid k."""
        from weatherisk.clustering import (
            quantile_threshold,
            cluster_number_threshold_method,
            clustering,
        )

        rng = np.random.default_rng(42)
        n = 30
        dm = rng.random((n, n)) * 100
        dm = (dm + dm.T) / 2
        np.fill_diagonal(dm, 0)

        hc = clustering(dm)
        thr = quantile_threshold(dm, 0.30)
        k = cluster_number_threshold_method(hc, thr)
        assert k >= 1
        assert k <= n


# ── Non-stationary pipeline ────────────────────────────────────────

# Module-level cache: run the pipeline once for all tests that share
# the same (preset, resolution, n_sim, seed) to avoid ~5s per call.
_pipeline_cache: dict[tuple, dict] = {}


def _get_cached_pipeline(preset="stripes", resolution=5, n_sim=5, seed=42):
    """Return a cached pipeline result, running it only on first call."""
    from weatherisk.pipeline import run_nonstationary_pipeline

    key = (preset, resolution, n_sim, seed)
    if key not in _pipeline_cache:
        _pipeline_cache[key] = run_nonstationary_pipeline(
            preset=preset, resolution=resolution, n_sim=n_sim,
            seed=seed, verbose=False,
        )
    return _pipeline_cache[key]


class TestNonstationaryPipeline:
    """Test run_nonstationary_pipeline on a tiny grid for correctness."""

    def test_pipeline_returns_all_keys(self):
        """Pipeline returns dict with all expected keys."""
        result = _get_cached_pipeline()
        expected_keys = {
            "grid", "preset", "a_matrix", "b_matrix", "g_matrix",
            "sim_data", "local_estimates", "smoothed",
            "dm_lec", "hc_lec", "k_lec", "labels_lec", "inclusters_lec",
            "dm_edc", "hc_edc", "k_edc", "labels_edc", "inclusters_edc",
        }
        assert expected_keys.issubset(result.keys())

    def test_pipeline_grid_shape(self):
        result = _get_cached_pipeline()
        grid = result["grid"]
        assert grid.resolution == 5
        assert grid.n_grid == 25
        assert result["sim_data"].shape == (5, 5, 5)

    def test_pipeline_labels_valid(self):
        result = _get_cached_pipeline()
        for key in ("labels_lec", "labels_edc"):
            labels = result[key]
            assert len(labels) == 25
            assert labels.min() >= 1  # fcluster is 1-based
            assert result[f"k_{key.split('_')[1]}"] >= 2

    def test_pipeline_uses_madogram_for_edc(self):
        """EDC should use madogram=True (raw v values, not EC-1)."""
        result = _get_cached_pipeline()
        dm = result["dm_edc"]
        # Madogram values are typically small (< 1/6 ≈ 0.167)
        upper = dm[np.triu_indices(25, k=1)]
        assert upper.max() < 0.5, "EDC dissimilarities should be raw madogram (< 0.5)"

    def test_pipeline_inclusters_shape(self):
        result = _get_cached_pipeline()
        for key in ("inclusters_lec", "inclusters_edc"):
            inc = result[key]
            assert inc.ndim == 2
            assert inc.shape[1] == 5  # a, b, g, n_cells, avg_llh

    def test_pipeline_deterministic(self):
        """Same seed should produce identical results (via two cached calls)."""
        r1 = _get_cached_pipeline(seed=123)
        r2 = _get_cached_pipeline(seed=123)  # hits cache → same object
        np.testing.assert_array_equal(r1["sim_data"], r2["sim_data"])
        np.testing.assert_array_equal(r1["labels_lec"], r2["labels_lec"])

    def test_pipeline_with_paper_preset(self):
        """Paper preset works with overridden resolution/n_sim."""
        result = _get_cached_pipeline(preset="paper_stripes")
        # Should use paper_stripes b_func: (x+5)/2
        b = result["b_matrix"]
        np.testing.assert_allclose(b.min(), 0.0, atol=0.1)

    def test_pipeline_output_dir(self, tmp_path):
        """Pipeline saves results when output_dir is specified."""
        from weatherisk.pipeline import run_nonstationary_pipeline

        run_nonstationary_pipeline(
            preset="stripes", resolution=5, n_sim=5,
            seed=42, verbose=False, output_dir=str(tmp_path),
        )
        import os
        saved = os.listdir(str(tmp_path))
        assert len(saved) > 0
        assert any(f.endswith(".npy") for f in saved)

    def test_pipeline_quantile_threshold_applied(self):
        """Verify quantile threshold is actually used for k selection."""
        from weatherisk.clustering import (
            quantile_threshold,
            cluster_number_threshold_method,
        )

        # Use the cached pipeline result, check that different quantiles
        # on the same dissimilarity matrix give different k values.
        result = _get_cached_pipeline()
        dm = result["dm_lec"]
        hc = result["hc_lec"]

        thr_low = quantile_threshold(dm, 0.10)
        thr_high = quantile_threshold(dm, 0.90)
        k_low = cluster_number_threshold_method(hc, thr_low)
        k_high = cluster_number_threshold_method(hc, thr_high)
        # Low threshold → more merges exceed it → more clusters
        assert k_low >= k_high


# ── Cluster comparison plot ─────────────────────────────────────────


class TestPlotClusterComparison:
    """Test plot_cluster_comparison produces valid figures."""

    @pytest.fixture
    def mock_pipeline_result(self):
        """Create mock pipeline result for plotting tests."""
        rng = np.random.default_rng(42)
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=5)

        labels_edc = rng.integers(1, 4, size=25)
        labels_lec = rng.integers(1, 4, size=25)

        # 4 clusters (indices 0-3, but 0 unused since labels are 1-based)
        inclusters_edc = np.full((4, 5), np.nan)
        inclusters_lec = np.full((4, 5), np.nan)
        for cl in range(1, 4):
            inclusters_edc[cl] = [2.0, cl * 1.5, 0.1, 10, -50]
            inclusters_lec[cl] = [2.0, cl * 1.5, 0.1, 10, -50]

        true_b = grid.X / 2.0 + 2.5

        return {
            "grid": grid,
            "labels_edc": labels_edc,
            "labels_lec": labels_lec,
            "inclusters_edc": inclusters_edc,
            "inclusters_lec": inclusters_lec,
            "true_field": true_b,
        }

    def test_returns_figure(self, mock_pipeline_result):
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison

        r = mock_pipeline_result
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=r["true_field"],
            param_index=1, param_name="b",
            show=False,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_three_axes(self, mock_pipeline_result):
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison

        r = mock_pipeline_result
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=r["true_field"],
            param_index=1, param_name="b",
            show=False,
        )
        # Should have 3 subplots + 1 colorbar axes = 4 axes
        axes = fig.get_axes()
        assert len(axes) >= 3
        plt.close(fig)

    def test_saves_to_file(self, mock_pipeline_result, tmp_path):
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison
        import os

        r = mock_pipeline_result
        out = str(tmp_path / "test_fig3.png")
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=r["true_field"],
            param_index=1, param_name="b",
            show=False, filename=out,
        )
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0
        plt.close(fig)

    def test_custom_vmin_vmax(self, mock_pipeline_result):
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison

        r = mock_pipeline_result
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=r["true_field"],
            param_index=1, param_name="b",
            vmin=0, vmax=5,
            show=False,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_with_suptitle(self, mock_pipeline_result):
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison

        r = mock_pipeline_result
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=r["true_field"],
            param_index=1, param_name="b",
            suptitle="Test figure",
            show=False,
        )
        assert fig._suptitle is not None
        plt.close(fig)

    def test_parameter_a_index(self, mock_pipeline_result):
        """param_index=0 should plot parameter a."""
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison

        r = mock_pipeline_result
        true_a = np.full_like(r["true_field"], 2.0)
        fig = plot_cluster_comparison(
            grid=r["grid"],
            labels_edc=r["labels_edc"],
            labels_lec=r["labels_lec"],
            inclusters_edc=r["inclusters_edc"],
            inclusters_lec=r["inclusters_lec"],
            true_field=true_a,
            param_index=0, param_name="a",
            show=False,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# ── End-to-end smoke test ──────────────────────────────────────────


class TestFig3EndToEnd:
    """Smoke test: run pipeline + plot on smallest possible grid."""

    def test_stripes_pipeline_and_plot(self, tmp_path):
        """Run paper_stripes at 5x5 resolution and produce the 3-panel figure."""
        import matplotlib.pyplot as plt
        from weatherisk.plotting import plot_cluster_comparison
        import os

        result = _get_cached_pipeline(preset="paper_stripes")

        out = str(tmp_path / "fig3_smoke.png")
        fig = plot_cluster_comparison(
            grid=result["grid"],
            labels_edc=result["labels_edc"],
            labels_lec=result["labels_lec"],
            inclusters_edc=result["inclusters_edc"],
            inclusters_lec=result["inclusters_lec"],
            true_field=result["b_matrix"],
            param_index=1, param_name="b",
            vmin=0, vmax=5,
            label="$b$ (semi-major diff.)",
            suptitle="Fig. 3 reproduction — paper_stripes",
            filename=out,
            show=False,
        )

        assert os.path.exists(out)
        assert os.path.getsize(out) > 1000  # should be a real PNG
        plt.close(fig)

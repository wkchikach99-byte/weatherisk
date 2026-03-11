from pathlib import Path

import numpy as np
import pytest

import weatherisk.cmip6_pipeline as cmip6_pipeline


def _patch_fast_cmip6_pipeline(monkeypatch, *, fail_local_once: bool) -> dict[str, int]:
    """Replace the heavy CMIP6 stages with tiny deterministic stand-ins."""
    state = {
        "load_calls": 0,
        "step2_calls": 0,
        "step3_calls": 0,
    }
    local_failure = {"armed": fail_local_once}

    pr = np.ones((24, 2, 2), dtype=float)
    times = np.arange(
        np.datetime64("2000-01"),
        np.datetime64("2002-01"),
        np.timedelta64(1, "M"),
    )
    lats = np.array([-10.0, 10.0], dtype=float)
    lons = np.array([0.0, 180.0], dtype=float)
    annual_max = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[1.5, 2.5], [3.5, 4.5]],
        ],
        dtype=float,
    )
    years = np.array([2000, 2001], dtype=int)
    frechet = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [1.5, 2.5, 3.5, 4.5],
        ],
        dtype=float,
    )
    valid_idx = np.array([0, 1, 2, 3], dtype=int)

    def fake_load_monthly_precipitation(*args, **kwargs):
        state["load_calls"] += 1
        return pr.copy(), times.copy(), lats.copy(), lons.copy()

    def fake_detrend_grid_fast(pr_in, period=12, *, n_workers=1, verbose=True):
        return pr_in.copy()

    def fake_monthly_annual_maxima(pr_detrended, times_in, verbose=True):
        return annual_max.copy(), years.copy()

    def fake_compute_frechet_global(annual_max_in, n_workers=1, verbose=True):
        state["step2_calls"] += 1
        return frechet.copy(), valid_idx.copy()

    def fake_run_local_estimation_cmip6(frechet_in, grid_coords, cfg, verbose=True):
        state["step3_calls"] += 1
        if local_failure["armed"]:
            local_failure["armed"] = False
            raise RuntimeError("simulated Step 3 crash")
        return np.array(
            [
                [1.0, 0.1, 0.0],
                [1.1, 0.1, 0.0],
                [1.2, 0.2, 0.1],
                [1.3, 0.2, 0.1],
            ],
            dtype=float,
        )

    def fake_smooth_estimates_cmip6(est, grid_coords, cfg, verbose=True):
        return est + np.array([0.0, 0.05, 0.0])

    def fake_run_clustering_cmip6(smoothed, frechet_in, cfg, verbose=True):
        labels_lec = np.array([1, 1, 2, 2], dtype=int)
        labels_edc = np.array([1, 2, 1, 2], dtype=int)
        dm_lec = np.zeros((4, 4), dtype=float)
        dm_edc = np.zeros((4, 4), dtype=float)
        hc = np.zeros((3, 4), dtype=float)
        return {
            "labels_lec": labels_lec,
            "k_lec": 2,
            "hc_lec": hc if cfg.retain_clustering_artifacts else None,
            "dm_lec": dm_lec if cfg.retain_clustering_artifacts else None,
            "labels_edc": labels_edc,
            "k_edc": 2,
            "hc_edc": hc if cfg.retain_clustering_artifacts else None,
            "dm_edc": dm_edc if cfg.retain_clustering_artifacts else None,
        }

    def fake_incluster_reestimate_cmip6(frechet_in, grid_coords, labels, cfg, tag, verbose=True):
        return {
            int(label): np.array([1.0 + 0.1 * int(label), 0.2, 0.0], dtype=float)
            for label in np.unique(labels)
        }

    monkeypatch.setattr(cmip6_pipeline, "load_monthly_precipitation", fake_load_monthly_precipitation)
    monkeypatch.setattr(cmip6_pipeline, "_detrend_grid_fast", fake_detrend_grid_fast)
    monkeypatch.setattr(cmip6_pipeline, "_monthly_annual_maxima", fake_monthly_annual_maxima)
    monkeypatch.setattr(cmip6_pipeline, "_compute_frechet_global", fake_compute_frechet_global)
    monkeypatch.setattr(cmip6_pipeline, "_run_local_estimation_cmip6", fake_run_local_estimation_cmip6)
    monkeypatch.setattr(cmip6_pipeline, "_smooth_estimates_cmip6", fake_smooth_estimates_cmip6)
    monkeypatch.setattr(cmip6_pipeline, "_run_clustering_cmip6", fake_run_clustering_cmip6)
    monkeypatch.setattr(cmip6_pipeline, "_incluster_reestimate_cmip6", fake_incluster_reestimate_cmip6)

    return state


def test_cmip6_pipeline_auto_resumes_and_cleans_checkpoints(tmp_path, monkeypatch):
    state = _patch_fast_cmip6_pipeline(monkeypatch, fail_local_once=True)
    output_dir = tmp_path / "cmip6_output"
    checkpoint_dir = output_dir / "checkpoints"

    cfg = cmip6_pipeline.CMIP6Config(
        data_dir=str(tmp_path / "input_data"),
        output_dir=str(output_dir),
        checkpoint_dir=str(checkpoint_dir),
        n_workers=1,
    )

    with pytest.raises(RuntimeError, match="simulated Step 3 crash"):
        cmip6_pipeline.run_cmip6_pipeline(cfg, verbose=False)

    assert state["load_calls"] == 1
    assert state["step2_calls"] == 1
    assert state["step3_calls"] == 1
    assert checkpoint_dir.exists()
    assert (checkpoint_dir / "step2.npz").exists()
    assert (checkpoint_dir / "manifest.json").exists()

    result = cmip6_pipeline.run_cmip6_pipeline(cfg, verbose=False)

    assert state["load_calls"] == 1
    assert state["step2_calls"] == 1
    assert state["step3_calls"] == 2
    assert not checkpoint_dir.exists()
    assert (output_dir / "pipeline_results.npz").exists()
    assert result["annual_max"].shape == (2, 2, 2)
    np.testing.assert_array_equal(result["years"], np.array([2000, 2001]))
    np.testing.assert_array_equal(result["valid_idx"], np.array([0, 1, 2, 3]))
    assert result["k_lec"] == 2
    assert result["k_edc"] == 2


def test_cmip6_pipeline_discards_clustering_artifacts_by_default(tmp_path, monkeypatch):
    _patch_fast_cmip6_pipeline(monkeypatch, fail_local_once=False)
    cfg = cmip6_pipeline.CMIP6Config(
        data_dir=str(tmp_path / "input_data"),
        output_dir=str(tmp_path / "cmip6_output"),
        n_workers=1,
    )

    result = cmip6_pipeline.run_cmip6_pipeline(cfg, verbose=False)

    assert result["dm_lec"] is None
    assert result["dm_edc"] is None
    assert result["hc_lec"] is None
    assert result["hc_edc"] is None


def test_cmip6_pipeline_retains_clustering_artifacts_when_requested(tmp_path, monkeypatch):
    _patch_fast_cmip6_pipeline(monkeypatch, fail_local_once=False)
    cfg = cmip6_pipeline.CMIP6Config(
        data_dir=str(tmp_path / "input_data"),
        output_dir=str(tmp_path / "cmip6_output"),
        n_workers=1,
        retain_clustering_artifacts=True,
    )

    result = cmip6_pipeline.run_cmip6_pipeline(cfg, verbose=False)

    assert result["dm_lec"].shape == (4, 4)
    assert result["dm_edc"].shape == (4, 4)
    assert result["hc_lec"].shape == (3, 4)
    assert result["hc_edc"].shape == (3, 4)
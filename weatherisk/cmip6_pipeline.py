"""CMIP6 climate data pipeline: AWI-ESM-1-1-LR → LEC / EDC clusters.

Reproduces Figure 9 from:
    Contzen et al. (2025), *Extremal dependence and local estimation
    clustering for non-stationary max-stable processes*, Extremes 28:713–737.

Paper §5 parameters:
    - Data: AWI-ESM-1-1-LR historical, monthly precipitation, 1850–2005
    - Grid: T63 (192×96), ~1.85° × 1.85°
    - Pre-processing: STL decomposition (Cleveland et al. 1990) to de-trend,
      then annual maxima of monthly data
    - Marginal transform: GEV → unit Fréchet at each grid point
    - Dependence parameters: ν=5, α=1, ε=5 (grid point distance)
    - Clustering cut-off: empirical 30%-quantile of pairwise dissimilarities
    - Result: k_EDC=104, k_LEC=24

Pipeline steps:
    1. Load monthly precipitation (auto-download if missing)
    2. STL de-trend → annual maxima of monthly de-trended data
    3. GEV fit → unit Fréchet transform per grid point
    4. Local pairwise composite likelihood MLE → (a, b, γ)
    5. Spatial smoothing
    6. LEC / EDC clustering (30%-quantile threshold)
    7. In-cluster re-estimation
    8. Generate cluster maps (Figure 9)
"""

from __future__ import annotations

import gc
import json
import os
import shutil
import subprocess
import tempfile
import time
import warnings
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import minimize
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import qmc, rankdata

from weatherisk.extremes import fit_gev, to_frechet
from weatherisk.density import pairwise_density_summand
from weatherisk.clustering import (
    calc_distance_ellipses,
    clustering,
    cluster_number_threshold_method,
)
from weatherisk.backend import (
    calc_distance_ellipses_condensed,
)
from weatherisk.cmip6_data import (
    ensure_cmip6_data,
    load_monthly_precipitation,
    DEFAULT_DATA_DIR,
)


# ══════════════════════════════════════════════════════════════════
#  Configuration
# ══════════════════════════════════════════════════════════════════

@dataclass
class CMIP6Config:
    """Configuration for the CMIP6 Figure 9 reproduction pipeline.

    Defaults are taken directly from the paper (§5).
    """

    # I/O
    data_dir: str = DEFAULT_DATA_DIR
    output_dir: str = "output/cmip6_fig9"

    # Time range
    year_start: int = 1850
    year_end: int = 2005     # inclusive

    # Max-stable hyper-parameters (paper §5)
    df: float = 5.0          # ν  (degrees of freedom)
    alpha: float = 1.0       # α  (smoothness exponent)

    # Local estimation — ε = 5 grid point distance
    neighbor_radius: float = 5.0   # ε  in grid-point units
    smoothing_radius: float = 2.0  # spatial smoothing in grid-point units
    mle_ensemble: int = 5          # multi-start restarts; strict parity status is tracked in docs/python_r_parity_migration_plan.md
    max_boundary_retries: int = 5  # retry budget when optimum hits bounds

    # Clustering (§5)
    quantile_threshold: float = 0.30  # 30%-quantile of pairwise dists

    # STL
    stl_period: int = 12  # months (annual seasonality)

    # Parallelism
    n_workers: int = 1

    # Clustering memory controls
    lec_chunk_size: int = 64
    retain_clustering_artifacts: bool = False

    # Checkpointing
    checkpoint_dir: str | None = None
    cleanup_checkpoints_on_success: bool = True

    @property
    def years(self) -> range:
        return range(self.year_start, self.year_end + 1)


_CHECKPOINT_STAGE_ORDER = (
    "step2",
    "step3",
    "step4",
    "step5",
    "step6_lec",
    "step6_edc",
)


def _checkpoint_signature(cfg: CMIP6Config) -> dict[str, Any]:
    """Return the config subset that must match for checkpoint reuse."""
    return {
        "data_dir": os.path.abspath(cfg.data_dir),
        "year_start": int(cfg.year_start),
        "year_end": int(cfg.year_end),
        "df": float(cfg.df),
        "alpha": float(cfg.alpha),
        "neighbor_radius": float(cfg.neighbor_radius),
        "smoothing_radius": float(cfg.smoothing_radius),
        "mle_ensemble": int(cfg.mle_ensemble),
        "quantile_threshold": float(cfg.quantile_threshold),
        "stl_period": int(cfg.stl_period),
    }


def _checkpoint_dir(cfg: CMIP6Config) -> Path:
    """Return the checkpoint directory for a CMIP6 run."""
    if cfg.checkpoint_dir is not None:
        return Path(cfg.checkpoint_dir)
    return Path(cfg.output_dir) / "checkpoints"


def _checkpoint_manifest_path(cfg: CMIP6Config) -> Path:
    return _checkpoint_dir(cfg) / "manifest.json"


def _checkpoint_data_path(cfg: CMIP6Config, stage: str) -> Path:
    return _checkpoint_dir(cfg) / f"{stage}.npz"


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    tmp_path = path.with_name(f"{path.stem}.tmp{path.suffix}")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True)
    os.replace(tmp_path, path)


def _atomic_savez(path: Path, **arrays: np.ndarray | int | float) -> None:
    tmp_path = path.with_name(f"{path.stem}.tmp.npz")
    np.savez_compressed(tmp_path, **arrays)
    os.replace(tmp_path, path)


def _cleanup_checkpoints(cfg: CMIP6Config, *, verbose: bool = True) -> None:
    checkpoint_dir = _checkpoint_dir(cfg)
    if checkpoint_dir.exists():
        shutil.rmtree(checkpoint_dir)
        if verbose:
            print(f"  Removed checkpoints from {checkpoint_dir}")


def _write_checkpoint(
    cfg: CMIP6Config,
    stage: str,
    *,
    verbose: bool = True,
    **arrays: np.ndarray | int | float,
) -> None:
    """Persist one pipeline stage so failed jobs can resume automatically."""
    checkpoint_dir = _checkpoint_dir(cfg)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    data_path = _checkpoint_data_path(cfg, stage)
    _atomic_savez(data_path, **arrays)
    _atomic_write_json(
        _checkpoint_manifest_path(cfg),
        {
            "signature": _checkpoint_signature(cfg),
            "stage": stage,
            "stages": {s: _checkpoint_data_path(cfg, s).name for s in _CHECKPOINT_STAGE_ORDER},
        },
    )
    if verbose:
        print(f"  Saved checkpoint: {stage}")


def _load_npz_dict(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as payload:
        return {key: payload[key] for key in payload.files}


def _load_checkpoint(cfg: CMIP6Config, *, verbose: bool = True) -> dict[str, Any] | None:
    """Load the latest compatible checkpoint if one exists."""
    manifest_path = _checkpoint_manifest_path(cfg)
    if not manifest_path.exists():
        return None

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        if verbose:
            print("  Ignoring unreadable checkpoint manifest; starting fresh.")
        _cleanup_checkpoints(cfg, verbose=False)
        return None

    if manifest.get("signature") != _checkpoint_signature(cfg):
        if verbose:
            print("  Ignoring stale checkpoints with mismatched configuration.")
        _cleanup_checkpoints(cfg, verbose=False)
        return None

    stage = manifest.get("stage")
    if stage not in _CHECKPOINT_STAGE_ORDER:
        if verbose:
            print("  Ignoring checkpoint with unknown stage; starting fresh.")
        _cleanup_checkpoints(cfg, verbose=False)
        return None

    required_stages = _CHECKPOINT_STAGE_ORDER[: _CHECKPOINT_STAGE_ORDER.index(stage) + 1]
    loaded: dict[str, Any] = {"stage": stage}
    try:
        for stage_name in required_stages:
            loaded[stage_name] = _load_npz_dict(_checkpoint_data_path(cfg, stage_name))
    except Exception:
        if verbose:
            print("  Ignoring incomplete checkpoints; starting fresh.")
        _cleanup_checkpoints(cfg, verbose=False)
        return None

    if verbose:
        print(f"  Resuming from checkpoint: {stage}")
    return loaded


def _npz_scalar(value: np.ndarray, caster: type[int] | type[float] = int) -> int | float:
    """Convert a 0-D npz payload back to a Python scalar."""
    return caster(np.asarray(value).item())


def _stage_at_least(stage: str | None, target: str) -> bool:
    """Return whether a checkpoint stage is at or beyond the target stage."""
    if stage is None:
        return False
    return _CHECKPOINT_STAGE_ORDER.index(stage) >= _CHECKPOINT_STAGE_ORDER.index(target)


def _serialize_clustering_checkpoint(clustering_result: dict[str, Any]) -> dict[str, np.ndarray | int]:
    """Convert Step 5 clustering outputs into an npz-friendly payload."""
    payload: dict[str, np.ndarray | int] = {
        "labels_lec": np.asarray(clustering_result["labels_lec"], dtype=int),
        "k_lec": int(clustering_result["k_lec"]),
        "labels_edc": np.asarray(clustering_result["labels_edc"], dtype=int),
        "k_edc": int(clustering_result["k_edc"]),
    }
    for key in ("hc_lec", "dm_lec", "hc_edc", "dm_edc"):
        value = clustering_result.get(key)
        if value is not None:
            payload[key] = np.asarray(value)
    return payload


def _deserialize_clustering_checkpoint(payload: dict[str, np.ndarray]) -> dict[str, Any]:
    """Restore Step 5 clustering outputs from a checkpoint payload."""
    return {
        "labels_lec": payload["labels_lec"],
        "k_lec": _npz_scalar(payload["k_lec"], int),
        "hc_lec": payload.get("hc_lec"),
        "dm_lec": payload.get("dm_lec"),
        "labels_edc": payload["labels_edc"],
        "k_edc": _npz_scalar(payload["k_edc"], int),
        "hc_edc": payload.get("hc_edc"),
        "dm_edc": payload.get("dm_edc"),
    }


def _serialize_incluster_checkpoint(results: dict[int, np.ndarray]) -> dict[str, np.ndarray]:
    """Convert partial Step 6 cluster estimates into an npz-friendly payload."""
    cluster_ids = np.array(sorted(results), dtype=int)
    if cluster_ids.size == 0:
        estimates = np.empty((0, 3), dtype=float)
    else:
        estimates = np.vstack([np.asarray(results[cl], dtype=float) for cl in cluster_ids])
    return {
        "cluster_ids": cluster_ids,
        "estimates": estimates,
    }


def _deserialize_incluster_checkpoint(payload: dict[str, np.ndarray]) -> dict[int, np.ndarray]:
    """Restore partial Step 6 cluster estimates from a checkpoint payload."""
    cluster_ids = np.asarray(payload.get("cluster_ids", np.empty(0, dtype=int)), dtype=int)
    estimates = np.asarray(payload.get("estimates", np.empty((0, 3), dtype=float)), dtype=float)
    return {
        int(cluster_id): np.asarray(estimates[index], dtype=float)
        for index, cluster_id in enumerate(cluster_ids)
    }


# ══════════════════════════════════════════════════════════════════
#  Module-level worker functions (must be picklable for multiprocessing)
# ══════════════════════════════════════════════════════════════════

# Shared state for MLE worker (set before Pool is created)
_MLE_FRECHET: np.ndarray | None = None
_MLE_GRID_COORDS: np.ndarray | None = None
_MLE_CFG: CMIP6Config | None = None
_GEV_AM_FLAT: np.ndarray | None = None
_INCLUSTER_FRECHET: np.ndarray | None = None
_INCLUSTER_GRID_COORDS: np.ndarray | None = None
_INCLUSTER_LABELS: np.ndarray | None = None
_INCLUSTER_CFG: CMIP6Config | None = None


def _rscript_path() -> str:
    candidate = shutil.which("Rscript")
    if candidate is None:
        raise RuntimeError(
            "Rscript not found on PATH. "
            "R is required for exact numerical parity."
        )
    return candidate


def _detrend_grid_via_r(pr: np.ndarray, period: int) -> np.ndarray:
    rscript = _rscript_path()

    n_months, n_lat, n_lon = pr.shape
    pr_flat = pr.reshape(n_months, n_lat * n_lon, order="F")
    r_code = """
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == '--args') {
  args <- args[-1]
}
input_path <- args[1]
output_path <- args[2]
period_val <- as.integer(args[3])
n_months <- as.integer(args[4])
n_cells <- as.integer(args[5])

pr <- as.matrix(read.csv(input_path, header = FALSE))
if (!all(dim(pr) == c(n_months, n_cells))) {
  stop('Unexpected input matrix dimensions')
}

detrended <- matrix(0, nrow = n_months, ncol = n_cells)
for (cell in 1:n_cells) {
  ts_data <- ts(pr[, cell], frequency = period_val)
  stl_fit <- stl(ts_data, s.window = 'periodic', robust = TRUE)
  detrended[, cell] <- pr[, cell] - stl_fit$time.series[, 'trend']
}

write.table(detrended, file = output_path, row.names = FALSE, col.names = FALSE, sep = ',')
"""

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "cmip6_pr.csv"
            output_path = Path(tmpdir) / "cmip6_detrended.csv"
            np.savetxt(input_path, pr_flat, delimiter=",")
            subprocess.run(
                [
                    rscript,
                    "-e",
                    r_code,
                    "--args",
                    str(input_path),
                    str(output_path),
                    str(int(period)),
                    str(int(n_months)),
                    str(int(n_lat * n_lon)),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            detrended_flat = np.loadtxt(output_path, delimiter=",")
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        raise RuntimeError(f"R STL detrending subprocess failed: {exc}") from exc

    detrended_flat = np.atleast_2d(np.asarray(detrended_flat, dtype=float))
    if detrended_flat.shape == (1, pr_flat.shape[0]):
        detrended_flat = detrended_flat.T
    if detrended_flat.shape != pr_flat.shape:
        raise RuntimeError(
            f"R STL output shape {detrended_flat.shape} != expected {pr_flat.shape}"
        )
    return detrended_flat.reshape(n_months, n_lat, n_lon, order="F")


def _mle_worker_init(frechet, grid_coords, cfg):
    """Initializer for MLE worker pool — sets shared data."""
    global _MLE_FRECHET, _MLE_GRID_COORDS, _MLE_CFG
    _MLE_FRECHET = frechet
    _MLE_GRID_COORDS = grid_coords
    _MLE_CFG = cfg


def _mle_worker(cidx):
    """Worker for parallel local MLE."""
    return cidx, _local_mle_one_cmip6(
        _MLE_FRECHET, cidx, _MLE_GRID_COORDS, _MLE_CFG
    )


def _gev_worker_init(am_flat):
    """Initializer for GEV worker pool."""
    global _GEV_AM_FLAT
    _GEV_AM_FLAT = am_flat


def _gev_worker(cidx):
    """Fit GEV and return Fréchet data for a single valid cell."""
    try:
        col = _GEV_AM_FLAT[:, cidx]
        nan_mask = np.isnan(col)
        if nan_mask.any():
            col = col.copy()
            col[nan_mask] = np.nanmean(col)

        loc, scale, shape = fit_gev(col)
        fr = to_frechet(col, loc, scale, shape)
        fr = np.clip(fr, 0.01, None)
        return cidx, fr, np.array([loc, scale, shape]), True
    except Exception:
        return cidx, None, np.array([np.nan, np.nan, np.nan]), False


def _incluster_worker_init(frechet, grid_coords, labels, cfg):
    """Initializer for in-cluster re-estimation worker pool."""
    global _INCLUSTER_FRECHET, _INCLUSTER_GRID_COORDS, _INCLUSTER_LABELS, _INCLUSTER_CFG
    _INCLUSTER_FRECHET = frechet
    _INCLUSTER_GRID_COORDS = grid_coords
    _INCLUSTER_LABELS = labels
    _INCLUSTER_CFG = cfg


def _incluster_estimate_one(
    frechet: np.ndarray,
    grid_coords: np.ndarray,
    labels: np.ndarray,
    cfg: CMIP6Config,
    cl: int,
) -> tuple[int, np.ndarray, int]:
    """Estimate one cluster for the in-cluster re-estimation step."""
    from weatherisk.backend import optimize_pairwise_density

    mask = labels == cl
    n_cl = int(mask.sum())
    if n_cl < 3:
        return cl, np.array([1.0, 0.0, 0.0]), n_cl

    z_cl = frechet[:, mask].T
    X_cl = grid_coords[mask, 1]
    Y_cl = grid_coords[mask, 0]

    try:
        est = optimize_pairwise_density(
            z_cl, cfg.df, cfg.alpha, X_cl, Y_cl,
            upper_bounds=(15.0, 15.0),
            max_dist=4.0 * cfg.neighbor_radius,
            ensemble=3,
        )
    except Exception:
        est = np.array([1.0, 0.0, 0.0])

    return cl, est, n_cl


def _incluster_worker(cl: int) -> tuple[int, np.ndarray, int]:
    """Worker for parallel in-cluster re-estimation."""
    return _incluster_estimate_one(
        _INCLUSTER_FRECHET,
        _INCLUSTER_GRID_COORDS,
        _INCLUSTER_LABELS,
        _INCLUSTER_CFG,
        cl,
    )


def _compute_frechet_one(am_flat: np.ndarray, cidx: int):
    """Serial helper matching the parallel GEV worker logic."""
    try:
        col = am_flat[:, cidx]
        nan_mask = np.isnan(col)
        if nan_mask.any():
            col = col.copy()
            col[nan_mask] = np.nanmean(col)

        loc, scale, shape = fit_gev(col)
        fr = to_frechet(col, loc, scale, shape)
        fr = np.clip(fr, 0.01, None)
        return fr, np.array([loc, scale, shape]), True
    except Exception:
        return None, np.array([np.nan, np.nan, np.nan]), False


# ══════════════════════════════════════════════════════════════════
#  Step 1.  De-trending (per-cell STL, parallel)
# ══════════════════════════════════════════════════════════════════

# Module-level worker for multiprocessing (must be picklable).
_stl_globals: dict = {}


def _stl_worker_init(period: int, seasonal_window: int) -> None:
    """Initializer for STL pool workers — stores shared params."""
    _stl_globals["period"] = period
    _stl_globals["seasonal_window"] = seasonal_window


def _stl_one(args):
    """Run STL on a single grid cell (i, j, timeseries)."""
    from statsmodels.tsa.seasonal import STL

    i, j, ts = args
    fit = STL(
        ts,
        period=_stl_globals["period"],
        seasonal=_stl_globals["seasonal_window"],
        seasonal_deg=0,
        robust=True,
    ).fit()
    return i, j, ts - fit.trend


def _detrend_grid_fast(
    pr: np.ndarray,
    period: int = 12,
    trend_window: int = 121,
    *,
    n_workers: int = 1,
    verbose: bool = True,
) -> np.ndarray:
    """De-trend precipitation grid using STL (Cleveland et al. 1990).

    Matches R's ``stl(ts, s.window = "periodic", robust = TRUE)``:
    the trend component is estimated via STL and subtracted, keeping
    the seasonal cycle and residual.

    Parameters
    ----------
    pr : ndarray, shape (n_months, n_lat, n_lon)
    period : int
        Seasonal period (12 for monthly data).
    trend_window : int
        Unused (kept for API compatibility). STL determines its own
        trend smoother width from the data length and period.
    n_workers : int
        Number of parallel workers (1 = serial).

    Returns
    -------
    ndarray, same shape as *pr*
        De-trended precipitation (original - trend).
    """
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 1a : De-trending (STL)")
        print("=" * 60)

    if n_workers == 1:
        detrended_r = _detrend_grid_via_r(pr, period)
        if verbose:
            print("  Using local R STL oracle for exact parity")
        return detrended_r

    t0 = time.time()
    n_months, n_lat, n_lon = pr.shape
    total = n_lat * n_lon

    # s.window="periodic" -> seasonal window = n_months (or larger),
    # seasonal_deg=0 to match R's constant seasonal LOESS.
    seasonal_window = n_months + 1 if n_months % 2 == 0 else n_months
    detrended = np.empty_like(pr)

    # Progress every ~5 % or every 500 cells, whichever is smaller;
    # but at least every 50 cells to avoid spamming on tiny grids.
    prog_interval = max(50, min(500, total // 20))

    if n_workers > 1:
        if verbose:
            print(f"  {total} cells, {n_workers} workers, "
                  f"seasonal={seasonal_window}")
        tasks = [
            (i, j, pr[:, i, j])
            for i in range(n_lat) for j in range(n_lon)
        ]
        # Scale chunksize with grid: ~4 chunks per worker for good balance.
        chunksize = max(1, total // (n_workers * 4))
        with Pool(
            n_workers,
            initializer=_stl_worker_init,
            initargs=(period, seasonal_window),
        ) as pool:
            for count, (i, j, dt) in enumerate(
                pool.imap_unordered(_stl_one, tasks, chunksize=chunksize)
            ):
                detrended[:, i, j] = dt
                done = count + 1
                if verbose and (done % prog_interval == 0 or done == total):
                    elapsed = time.time() - t0
                    rate = done / elapsed if elapsed > 0 else 0
                    eta = (total - done) / rate if rate > 0 else 0
                    print(f"    {done:6d}/{total} "
                          f"({100 * done / total:5.1f}%)  "
                          f"{rate:.0f} cells/s  "
                          f"ETA {eta:.0f}s  "
                          f"[{elapsed:.0f}s]")
    else:
        # Serial path: set globals directly in this process.
        _stl_globals["period"] = period
        _stl_globals["seasonal_window"] = seasonal_window
        for count, (i, j) in enumerate(
            (i, j) for i in range(n_lat) for j in range(n_lon)
        ):
            _, _, dt = _stl_one((i, j, pr[:, i, j]))
            detrended[:, i, j] = dt
            done = count + 1
            if verbose and (done % prog_interval == 0 or done == total):
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0
                print(f"    {done:6d}/{total} "
                      f"({100 * done / total:5.1f}%)  "
                      f"{rate:.0f} cells/s  "
                      f"ETA {eta:.0f}s  "
                      f"[{elapsed:.0f}s]")

    if verbose:
        elapsed = time.time() - t0
        print(f"  Grid: {n_lat} x {n_lon} = {total} cells")
        print(f"  STL period={period}, seasonal={seasonal_window}, robust=True")
        print(f"  Workers: {n_workers}")
        print(f"  Elapsed: {elapsed:.1f}s")

    return detrended


# ══════════════════════════════════════════════════════════════════
#  Step 1b.  Annual maxima of monthly data
# ══════════════════════════════════════════════════════════════════

def _monthly_annual_maxima(
    pr_detrended: np.ndarray,
    times: np.ndarray,
    *,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute annual maxima of monthly de-trended data.

    Parameters
    ----------
    pr_detrended : ndarray, shape (n_months, n_lat, n_lon)
    times : ndarray of datetime64

    Returns
    -------
    annual_max : ndarray, shape (n_years, n_lat, n_lon)
    years : ndarray, shape (n_years,)
    """
    import pandas as pd

    if verbose:
        print("\n" + "=" * 60)
        print("  Step 1b : Annual maxima of monthly de-trended data")
        print("=" * 60)

    # Extract years
    time_years = pd.DatetimeIndex(times).year.values
    unique_years = np.sort(np.unique(time_years))

    # Keep only complete years (12 months)
    complete_years = [y for y in unique_years
                      if np.sum(time_years == y) == 12]
    complete_years = np.array(complete_years)

    n_years = len(complete_years)
    n_lat, n_lon = pr_detrended.shape[1], pr_detrended.shape[2]
    annual_max = np.empty((n_years, n_lat, n_lon))

    for k, year in enumerate(complete_years):
        mask = time_years == year
        annual_max[k] = pr_detrended[mask].max(axis=0)

    if verbose:
        print(f"  {n_years} complete years: "
              f"{complete_years[0]}–{complete_years[-1]}")
        print(f"  Annual max shape: {annual_max.shape}")
        print(f"  Range: [{annual_max.min():.4f}, {annual_max.max():.4f}]")

    return annual_max, complete_years


# ══════════════════════════════════════════════════════════════════
#  Step 2.  GEV → Fréchet transform
# ══════════════════════════════════════════════════════════════════

def _compute_frechet_global(
    annual_max: np.ndarray,
    n_workers: int = 1,
    *,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit GEV and transform to unit Fréchet at every grid point.

    Parameters
    ----------
    annual_max : ndarray, shape (n_years, n_lat, n_lon)

    Returns
    -------
    frechet : ndarray, shape (n_years, n_cells)
        Fréchet-transformed data for valid (land) cells.
    valid_idx : ndarray, shape (n_valid,)
        Linear indices (row-major) into (n_lat, n_lon) grid.
    """
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 2 : GEV fit → unit Fréchet transform")
        print("=" * 60)

    n_years, n_lat, n_lon = annual_max.shape
    n_cells = n_lat * n_lon
    am_flat = annual_max.reshape(n_years, n_cells)

    # Identify valid cells (not all-NaN, positive variance)
    valid_mask = (~np.all(np.isnan(am_flat), axis=0)) & (np.nanstd(am_flat, axis=0) >= 1e-15)

    valid_idx = np.where(valid_mask)[0]
    n_valid = len(valid_idx)

    if verbose:
        print(f"  Total cells: {n_cells}   Valid: {n_valid} "
              f"({100 * n_valid / n_cells:.0f}%)")

    frechet = np.empty((n_years, n_valid))
    gev_params = np.empty((n_valid, 3))  # loc, scale, shape
    n_fail = 0

    t0 = time.time()
    if n_workers > 1:
        with Pool(
            n_workers,
            initializer=_gev_worker_init,
            initargs=(am_flat,),
        ) as pool:
            for k, result in enumerate(
                pool.imap_unordered(_gev_worker, valid_idx, chunksize=32)
            ):
                cidx, fr, params, success = result
                out_idx = np.searchsorted(valid_idx, cidx)
                if success:
                    frechet[:, out_idx] = fr
                else:
                    frechet[:, out_idx] = np.nan
                    n_fail += 1
                gev_params[out_idx] = params
                if verbose and (k + 1) % 2000 == 0:
                    print(f"    {k + 1:5d}/{n_valid} ({time.time() - t0:.0f}s)", flush=True)
    else:
        for k, c in enumerate(valid_idx):
            if verbose and k % 2000 == 0 and k > 0:
                print(f"    {k:5d}/{n_valid} ({time.time() - t0:.0f}s)", flush=True)
            try:
                fr, params, success = _compute_frechet_one(am_flat, c)
                if success:
                    frechet[:, k] = fr
                else:
                    frechet[:, k] = np.nan
                    n_fail += 1
                gev_params[k] = params
            except Exception:
                frechet[:, k] = np.nan
                gev_params[k] = [np.nan, np.nan, np.nan]
                n_fail += 1

    # Remove cells where GEV failed
    ok = np.all(np.isfinite(frechet), axis=0)
    frechet = frechet[:, ok]
    valid_idx = valid_idx[ok]
    gev_params = gev_params[ok]

    # Clip extreme tails (Padoan et al. 2010)
    frechet = np.clip(frechet, 0.05, float(n_years ** 2))

    if verbose:
        print(f"  GEV failures: {n_fail}")
        print(f"  Final valid cells: {len(valid_idx)}")
        print(f"  Fréchet shape: {frechet.shape}  "
              f"range [{frechet.min():.2f}, {frechet.max():.2f}]")
        print(f"  Elapsed: {time.time() - t0:.0f}s")

    return frechet, valid_idx


# ══════════════════════════════════════════════════════════════════
#  Step 3.  Local MLE  (grid-point distance)
# ══════════════════════════════════════════════════════════════════

def _grid_coords(
    valid_idx: np.ndarray,
    n_lat: int,
    n_lon: int,
) -> np.ndarray:
    """Convert linear indices to (row, col) grid coordinates.

    Returns array of shape (n_valid, 2) with (i_lat, j_lon) in grid units.
    Note: ε=5 means 5 grid-point distance, not physical distance.
    """
    i_lat = valid_idx // n_lon
    j_lon = valid_idx % n_lon
    return np.column_stack([i_lat.astype(float), j_lon.astype(float)])


def _local_mle_one_cmip6(
    frechet: np.ndarray,
    cidx: int,
    grid_coords: np.ndarray,
    cfg: CMIP6Config,
) -> np.ndarray:
    """Estimate (a, b, γ) at a single grid cell using grid-point distance."""
    n_years = frechet.shape[0]

    # Match the original R stencil enumeration order exactly: x-offset varies
    # fastest, y-offset varies slowest, then in-bounds cells are retained.
    # This keeps the local likelihood accumulation order identical to the
    # historical R workflow, which matters on flat objective surfaces.
    y_levels = list(dict.fromkeys(float(y) for y in grid_coords[:, 0]))
    x_levels = list(dict.fromkeys(float(x) for x in grid_coords[:, 1]))
    y_to_row = {y: row for row, y in enumerate(y_levels)}
    x_to_col = {x: col for col, x in enumerate(x_levels)}
    index_to_cell = {
        (y_to_row[float(y)], x_to_col[float(x)]): idx
        for idx, (y, x) in enumerate(grid_coords)
    }

    centre_y = float(grid_coords[cidx, 0])
    centre_x = float(grid_coords[cidx, 1])
    centre_row = y_to_row[centre_y]
    centre_col = x_to_col[centre_x]
    max_offset = int(np.floor(cfg.neighbor_radius + 1e-12))

    nb_list: list[int] = []
    for dy in range(-max_offset, max_offset + 1):
        for dx in range(-max_offset, max_offset + 1):
            if 0 < dx * dx + dy * dy <= cfg.neighbor_radius * cfg.neighbor_radius:
                nb_idx = index_to_cell.get((centre_row + dy, centre_col + dx))
                if nb_idx is not None:
                    nb_list.append(nb_idx)

    nb = np.asarray(nb_list, dtype=int)

    if len(nb) < 3:
        return np.array([1.0, 0.0, 0.0])

    z_c = frechet[:, cidx]
    zi = frechet[:, nb].T.reshape(-1)
    zj = np.tile(z_c, len(nb))
    xl = np.repeat(grid_coords[nb, 1] - centre_x, n_years).astype(np.float64)
    yl = np.repeat(grid_coords[nb, 0] - centre_y, n_years).astype(np.float64)

    good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
    zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]
    if len(zi) < 5:
        return np.array([1.0, 0.0, 0.0])

    from weatherisk.backend import optimize_local_mle

    return optimize_local_mle(
        zi, zj, xl, yl, cfg.df, cfg.alpha,
        ensemble=cfg.mle_ensemble, seed=42,
        max_boundary_retries=cfg.max_boundary_retries,
    )


def _run_local_estimation_cmip6(
    frechet: np.ndarray,
    grid_coords: np.ndarray,
    cfg: CMIP6Config,
    *,
    verbose: bool = True,
) -> np.ndarray:
    """Local MLE at every valid cell (parallelisable)."""
    if verbose:
        print("\n" + "=" * 60)
        print(f"  Step 3 : Local MLE  (ν={cfg.df}, α={cfg.alpha}, "
              f"ε={cfg.neighbor_radius})")
        print("=" * 60)

    # The historical R mini-reference is effectively serialized at a fixed
    # decimal precision before the extremely flat local likelihoods are
    # optimized. Normalizing the Fréchet inputs at this boundary avoids
    # last-bit drift sending L-BFGS-B to a different equivalent optimum.
    frechet = np.round(frechet, 11)

    n = frechet.shape[1]
    est = np.zeros((n, 3))

    if cfg.n_workers > 1:
        t0 = time.time()
        with Pool(
            cfg.n_workers,
            initializer=_mle_worker_init,
            initargs=(frechet, grid_coords, cfg),
        ) as pool:
            for count, (cidx, p) in enumerate(
                pool.imap_unordered(_mle_worker, range(n), chunksize=20)
            ):
                est[cidx] = p
                if verbose and (count + 1) % max(1, n // 20) == 0:
                    print(f"    {count + 1:5d}/{n} "
                          f"({time.time() - t0:.0f}s)",
                          flush=True)
    else:
        t0 = time.time()
        for c in range(n):
            if verbose and c % max(1, n // 20) == 0:
                print(f"    {c + 1:5d}/{n}  ({time.time() - t0:.0f}s)",
                      flush=True)
            est[c] = _local_mle_one_cmip6(frechet, c, grid_coords, cfg)

    if verbose:
        print(f"  a ∈ [{est[:, 0].min():.3f}, {est[:, 0].max():.3f}]")
        print(f"  b ∈ [{est[:, 1].min():.3f}, {est[:, 1].max():.3f}]")
        print(f"  γ ∈ [{np.degrees(est[:, 2]).min():.1f}°, "
              f"{np.degrees(est[:, 2]).max():.1f}°]")
        print(f"  Elapsed: {time.time() - t0:.0f}s")

    return est


# ══════════════════════════════════════════════════════════════════
#  Step 4.  Spatial smoothing
# ══════════════════════════════════════════════════════════════════

def _smooth_estimates_cmip6(
    est: np.ndarray,
    grid_coords: np.ndarray,
    cfg: CMIP6Config,
    *,
    verbose: bool = True,
) -> np.ndarray:
    """Spatial moving-average smoothing in grid-point distance units."""
    if verbose:
        print("\n" + "=" * 60)
        print(f"  Step 4 : Spatial smoothing  "
              f"(radius {cfg.smoothing_radius} grid pts)")
        print("=" * 60)

    n = len(est)
    out = np.empty_like(est)

    for c in range(n):
        d = np.sqrt(
            (grid_coords[:, 0] - grid_coords[c, 0]) ** 2 +
            (grid_coords[:, 1] - grid_coords[c, 1]) ** 2
        )
        nb = np.where(d <= cfg.smoothing_radius)[0]
        nbe = est[nb]

        out[c, 0] = nbe[:, 0].mean()
        out[c, 1] = nbe[:, 1].mean()

        # Angular wrapping for γ
        cg = est[c, 2]
        ang = nbe[:, 2].copy()
        ang = np.where(ang < cg - np.pi / 2, ang + np.pi, ang)
        ang = np.where(ang > cg + np.pi / 2, ang - np.pi, ang)
        mg = ang.mean()
        if mg < -np.pi / 2:
            mg += np.pi
        elif mg > np.pi / 2:
            mg -= np.pi
        out[c, 2] = mg

    if verbose:
        print(f"  Smoothed a ∈ [{out[:, 0].min():.3f}, "
              f"{out[:, 0].max():.3f}]")
        print(f"  Smoothed b ∈ [{out[:, 1].min():.3f}, "
              f"{out[:, 1].max():.3f}]")

    return out


# ══════════════════════════════════════════════════════════════════
#  Step 5.  LEC / EDC clustering
# ══════════════════════════════════════════════════════════════════

def _edc_matrix_flat(frechet: np.ndarray) -> np.ndarray:
    """Rank-based madogram dissimilarity matrix for flat data.

    Returns the raw madogram values v(i,j) as used in Contzen et al. (2025)
    for clustering.  These are the rank-based pairwise distances in [0, 0.5).

    Parameters
    ----------
    frechet : ndarray, shape (n_years, n_cells)
    """
    n_years, n_cells = frechet.shape
    ranks = np.empty((n_cells, n_years))
    for s in range(n_cells):
        ranks[s] = rankdata(frechet[:, s])

    diff_sum = cdist(ranks, ranks, metric="cityblock")
    v = diff_sum / (n_years * 2.0 * (n_years + 1))
    np.fill_diagonal(v, 0.0)
    return v


def _edc_condensed_flat(frechet: np.ndarray) -> np.ndarray:
    """Condensed rank-based madogram dissimilarity distances.

    Returns the raw madogram values (upper triangle) matching
    :func:`_edc_matrix_flat`, for use with scipy linkage.
    """
    n_years, n_cells = frechet.shape
    ranks = np.empty((n_cells, n_years))
    for s in range(n_cells):
        ranks[s] = rankdata(frechet[:, s])

    diff_sum = pdist(ranks, metric="cityblock")
    return diff_sum / (n_years * 2.0 * (n_years + 1))


def _run_clustering_cmip6(
    smoothed: np.ndarray,
    frechet: np.ndarray,
    cfg: CMIP6Config,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """LEC (D₂) and EDC (D₁) clustering with quantile threshold."""
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 5 : LEC & EDC clustering")
        print("=" * 60)

    # --- LEC ---
    if verbose:
        print("  Computing LEC dissimilarity (ellipse overlap D₂) …")
    t0 = time.time()
    if cfg.retain_clustering_artifacts:
        dm_lec = calc_distance_ellipses(
            smoothed,
            res=21,
            chunk_size=cfg.lec_chunk_size,
        )
        lec_condensed = squareform(dm_lec, checks=False)
        dm_lec_out = dm_lec
    else:
        lec_condensed = calc_distance_ellipses_condensed(
            smoothed,
            res=21,
            chunk_size=cfg.lec_chunk_size,
        )
        dm_lec_out = None
    thr_lec = float(np.quantile(lec_condensed, cfg.quantile_threshold))
    hc_lec = clustering(lec_condensed)
    k_lec = cluster_number_threshold_method(hc_lec, thr_lec)
    k_lec = max(2, k_lec)
    labels_lec = fcluster(hc_lec, t=k_lec, criterion="maxclust")
    hc_lec_out = hc_lec if cfg.retain_clustering_artifacts else None
    del lec_condensed
    if not cfg.retain_clustering_artifacts:
        del hc_lec
        gc.collect()
    if verbose:
        print(f"  LEC → k = {k_lec}  "
              f"(30%-quantile threshold = {thr_lec:.3f})  "
              f"[{time.time() - t0:.0f}s]")

    # --- EDC ---
    if verbose:
        print("  Computing EDC dissimilarity (madogram D₁) …")
    t0 = time.time()
    if cfg.retain_clustering_artifacts:
        dm_edc = _edc_matrix_flat(frechet)
        edc_condensed = squareform(dm_edc, checks=False)
        dm_edc_out = dm_edc
    else:
        edc_condensed = _edc_condensed_flat(frechet)
        dm_edc_out = None
    thr_edc = float(np.quantile(edc_condensed, cfg.quantile_threshold))
    hc_edc = clustering(edc_condensed)
    k_edc = cluster_number_threshold_method(hc_edc, thr_edc)
    k_edc = max(2, k_edc)
    labels_edc = fcluster(hc_edc, t=k_edc, criterion="maxclust")
    hc_edc_out = hc_edc if cfg.retain_clustering_artifacts else None
    del edc_condensed
    if not cfg.retain_clustering_artifacts:
        del hc_edc
        gc.collect()
    if verbose:
        print(f"  EDC → k = {k_edc}  "
              f"(30%-quantile threshold = {thr_edc:.5f})  "
              f"[{time.time() - t0:.0f}s]")

    return dict(
        labels_lec=labels_lec, k_lec=k_lec,
        hc_lec=hc_lec_out, dm_lec=dm_lec_out,
        labels_edc=labels_edc, k_edc=k_edc,
        hc_edc=hc_edc_out, dm_edc=dm_edc_out,
    )


# ══════════════════════════════════════════════════════════════════
#  Step 6.  In-cluster re-estimation
# ══════════════════════════════════════════════════════════════════

def _incluster_reestimate_cmip6(
    frechet: np.ndarray,
    grid_coords: np.ndarray,
    labels: np.ndarray,
    cfg: CMIP6Config,
    tag: str,
    *,
    initial_results: dict[int, np.ndarray] | None = None,
    checkpoint_stage: str | None = None,
    verbose: bool = True,
) -> dict[int, np.ndarray]:
    """Re-estimate (a, b, γ) per cluster via global pairwise CL MLE."""
    unique_cl = sorted(np.unique(labels))
    results: dict[int, np.ndarray] = {}
    if initial_results is not None:
        results = {
            int(cluster_id): np.asarray(estimate, dtype=float)
            for cluster_id, estimate in initial_results.items()
        }
    pending_cl = [cluster_id for cluster_id in unique_cl if cluster_id not in results]

    if verbose:
        print(f"  {tag}: re-estimating {len(unique_cl)} clusters …")
        if results:
            print(
                f"  {tag}: resuming {len(results)}/{len(unique_cl)} completed "
                "clusters from checkpoint."
            )

    if not pending_cl:
        return results

    if cfg.n_workers > 1 and len(pending_cl) > 1:
        with Pool(
            cfg.n_workers,
            initializer=_incluster_worker_init,
            initargs=(frechet, grid_coords, labels, cfg),
        ) as pool:
            for cl, est, n_cl in pool.imap_unordered(
                _incluster_worker, pending_cl, chunksize=1
            ):
                results[cl] = est
                if checkpoint_stage is not None:
                    _write_checkpoint(
                        cfg,
                        checkpoint_stage,
                        verbose=False,
                        **_serialize_incluster_checkpoint(results),
                    )
                if verbose:
                    print(f"    cl {cl:3d} ({n_cl:4d} cells)  "
                          f"a={est[0]:.3f}  b={est[1]:.3f}  "
                          f"γ={np.degrees(est[2]):.1f}°",
                          flush=True)
    else:
        for cl in pending_cl:
            cl_out, est, n_cl = _incluster_estimate_one(
                frechet, grid_coords, labels, cfg, cl
            )
            results[cl_out] = est
            if checkpoint_stage is not None:
                _write_checkpoint(
                    cfg,
                    checkpoint_stage,
                    verbose=False,
                    **_serialize_incluster_checkpoint(results),
                )
            if verbose:
                print(f"    cl {cl_out:3d} ({n_cl:4d} cells)  "
                      f"a={est[0]:.3f}  b={est[1]:.3f}  "
                      f"γ={np.degrees(est[2]):.1f}°")
    return results


# ══════════════════════════════════════════════════════════════════
#  Public orchestrator
# ══════════════════════════════════════════════════════════════════

def run_cmip6_pipeline(
    cfg: CMIP6Config | None = None,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Execute the full CMIP6 Figure 9 reproduction pipeline.

    Parameters
    ----------
    cfg : CMIP6Config, optional
        Pipeline configuration.  Uses paper defaults if *None*.
    verbose : bool
        Print progress banners.

    Returns
    -------
    dict
        Keys: ``lats``, ``lons``, ``frechet``, ``valid_idx``, ``grid_coords``,
        ``smoothed``, ``estimates``, ``labels_lec``, ``labels_edc``,
        ``k_lec``, ``k_edc``, ``annual_max``, ``years``, ``config``, etc.
    """
    if cfg is None:
        cfg = CMIP6Config()

    t_start = time.time()
    if verbose:
        print("\n" + "=" * 60)
        print("  CMIP6 Figure 9 Reproduction Pipeline")
        print("  AWI-ESM-1-1-LR historical precipitation")
        print(f"  ν={cfg.df}, α={cfg.alpha}, ε={cfg.neighbor_radius}")
        print("=" * 60)

    os.makedirs(cfg.output_dir, exist_ok=True)
    checkpoint = _load_checkpoint(cfg, verbose=verbose)
    checkpoint_stage = checkpoint["stage"] if checkpoint is not None else None

    if checkpoint is None:
        # Step 0: Ensure data
        pr, times, lats, lons = load_monthly_precipitation(
            cfg.data_dir,
            year_start=cfg.year_start,
            year_end=cfg.year_end,
            verbose=verbose,
        )

        # Step 1a: de-trend (parallel STL)
        pr_detrended = _detrend_grid_fast(
            pr, period=cfg.stl_period,
            n_workers=cfg.n_workers, verbose=verbose,
        )
        del pr
        gc.collect()

        # Step 1b: Annual maxima
        annual_max, years = _monthly_annual_maxima(
            pr_detrended, times, verbose=verbose,
        )
        del pr_detrended
        del times
        gc.collect()

        # Step 2: GEV → Fréchet
        frechet, valid_idx = _compute_frechet_global(
            annual_max, n_workers=cfg.n_workers, verbose=verbose,
        )
        n_lat, n_lon = annual_max.shape[1], annual_max.shape[2]
        _write_checkpoint(
            cfg,
            "step2",
            verbose=verbose,
            lats=lats,
            lons=lons,
            years=years,
            annual_max=annual_max,
            frechet=frechet,
            valid_idx=valid_idx,
            n_lat=n_lat,
            n_lon=n_lon,
        )
    else:
        step2 = checkpoint["step2"]
        lats = step2["lats"]
        lons = step2["lons"]
        years = step2["years"]
        annual_max = step2["annual_max"]
        frechet = step2["frechet"]
        valid_idx = step2["valid_idx"]
        n_lat = _npz_scalar(step2["n_lat"], int)
        n_lon = _npz_scalar(step2["n_lon"], int)

    # Grid coordinates in grid-point units
    grid_coords = _grid_coords(valid_idx, n_lat, n_lon)

    # Step 3: Local MLE
    if _stage_at_least(checkpoint_stage, "step3"):
        est = checkpoint["step3"]["estimates"]
    else:
        est = _run_local_estimation_cmip6(
            frechet, grid_coords, cfg, verbose=verbose,
        )
        _write_checkpoint(cfg, "step3", verbose=verbose, estimates=est)

    # Step 4: Smoothing
    if _stage_at_least(checkpoint_stage, "step4"):
        sm = checkpoint["step4"]["smoothed"]
    else:
        sm = _smooth_estimates_cmip6(
            est, grid_coords, cfg, verbose=verbose,
        )
        _write_checkpoint(cfg, "step4", verbose=verbose, smoothed=sm)

    # Step 5: Clustering
    if _stage_at_least(checkpoint_stage, "step5"):
        cl = _deserialize_clustering_checkpoint(checkpoint["step5"])
    else:
        cl = _run_clustering_cmip6(sm, frechet, cfg, verbose=verbose)
        _write_checkpoint(
            cfg,
            "step5",
            verbose=verbose,
            **_serialize_clustering_checkpoint(cl),
        )

    # Step 6: In-cluster re-estimation
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 6 : In-cluster re-estimation")
        print("=" * 60)
    lec_resume = (
        _deserialize_incluster_checkpoint(checkpoint["step6_lec"])
        if _stage_at_least(checkpoint_stage, "step6_lec")
        else None
    )
    incl_lec = _incluster_reestimate_cmip6(
        frechet, grid_coords, cl["labels_lec"], cfg, "LEC",
        initial_results=lec_resume,
        checkpoint_stage="step6_lec",
        verbose=verbose,
    )
    edc_resume = (
        _deserialize_incluster_checkpoint(checkpoint["step6_edc"])
        if _stage_at_least(checkpoint_stage, "step6_edc")
        else None
    )
    incl_edc = _incluster_reestimate_cmip6(
        frechet, grid_coords, cl["labels_edc"], cfg, "EDC",
        initial_results=edc_resume,
        checkpoint_stage="step6_edc",
        verbose=verbose,
    )

    elapsed = time.time() - t_start
    if verbose:
        print(f"\n  Pipeline completed in {elapsed:.0f}s "
              f"({elapsed / 60:.1f} min)")
        print(f"  k_LEC = {cl['k_lec']}  (paper: 24)")
        print(f"  k_EDC = {cl['k_edc']}  (paper: 104)")

    # Save final results
    save_path = os.path.join(cfg.output_dir, "pipeline_results.npz")
    np.savez_compressed(
        save_path,
        lats=lats, lons=lons,
        frechet=frechet,
        valid_idx=valid_idx,
        labels_lec=cl["labels_lec"],
        labels_edc=cl["labels_edc"],
        smoothed=sm,
        estimates=est,
        years=years,
        annual_max_shape=annual_max.shape,
    )
    if verbose:
        print(f"  Saved results to {save_path}")

    if cfg.cleanup_checkpoints_on_success:
        _cleanup_checkpoints(cfg, verbose=verbose)

    return dict(
        lats=lats,
        lons=lons,
        frechet=frechet,
        annual_max=annual_max,
        years=years,
        valid_idx=valid_idx,
        n_lat=n_lat,
        n_lon=n_lon,
        grid_coords=grid_coords,
        estimates=est,
        smoothed=sm,
        labels_lec=cl["labels_lec"],
        labels_edc=cl["labels_edc"],
        k_lec=cl["k_lec"],
        k_edc=cl["k_edc"],
        hc_lec=cl["hc_lec"],
        hc_edc=cl["hc_edc"],
        dm_lec=cl["dm_lec"],
        dm_edc=cl["dm_edc"],
        incl_lec=incl_lec,
        incl_edc=incl_edc,
        config=cfg,
    )


# ══════════════════════════════════════════════════════════════════
#  Plot generation — Figure 9
# ══════════════════════════════════════════════════════════════════

def plot_figure9(
    result: dict[str, Any],
    *,
    save_dir: str | None = None,
    dpi: int = 300,
    verbose: bool = True,
) -> list[str]:
    """Generate Figure 9: Global cluster maps for EDC and LEC.

    Parameters
    ----------
    result : dict
        Output of :func:`run_cmip6_pipeline`.
    save_dir : str, optional
        Directory for output files (default: result['config'].output_dir).
    dpi : int
        Figure resolution.
    verbose : bool
        Print saved paths.

    Returns
    -------
    list of str
        Paths to saved figures.
    """
    def _configure_https_certificates() -> None:
        """Use certifi's CA bundle for urllib-based Cartopy downloads."""
        try:
            import certifi
            import ssl
        except ImportError:
            return

        cafile = certifi.where()
        os.environ.setdefault("SSL_CERT_FILE", cafile)
        os.environ.setdefault("REQUESTS_CA_BUNDLE", cafile)
        os.environ.setdefault("CURL_CA_BUNDLE", cafile)
        ssl._create_default_https_context = (
            lambda: ssl.create_default_context(cafile=cafile)
        )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    _configure_https_certificates()

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        HAS_CARTOPY = True
    except ImportError:
        HAS_CARTOPY = False

    use_cartopy = HAS_CARTOPY
    if HAS_CARTOPY:
        try:
            from cartopy.io import shapereader

            shapereader.natural_earth(
                resolution="110m",
                category="physical",
                name="coastline",
            )
        except Exception as exc:
            use_cartopy = False
            if verbose:
                print(
                    "  Cartopy coastline data unavailable; "
                    f"falling back to plain lon/lat plots ({exc})."
                )

    cfg: CMIP6Config = result["config"]
    out = save_dir or cfg.output_dir
    os.makedirs(out, exist_ok=True)
    saved: list[str] = []

    lats = result["lats"]
    lons = result["lons"]
    n_lat = result["n_lat"]
    n_lon = result["n_lon"]
    valid_idx = result["valid_idx"]

    def _make_cmap(k):
        base = (list(plt.get_cmap("tab20").colors) +
                list(plt.get_cmap("Set3").colors) +
                list(plt.get_cmap("tab20b").colors) +
                list(plt.get_cmap("tab20c").colors))
        return ListedColormap([base[i % len(base)] for i in range(max(k, 1))])

    def _labels_to_grid(labels):
        """Map 1-D labels (per valid cell) onto (n_lat, n_lon) grid."""
        grid = np.full(n_lat * n_lon, np.nan)
        grid[valid_idx] = labels.astype(float)
        return grid.reshape(n_lat, n_lon)

    def _plot_panel(labels, k, title, filename):
        grid_2d = _labels_to_grid(labels)
        lon2d, lat2d = np.meshgrid(lons, lats)

        if use_cartopy:
            fig, ax = plt.subplots(
                1, 1, figsize=(14, 7),
                subplot_kw={"projection": ccrs.Robinson()},
            )
            cmap = _make_cmap(k)
            mesh = ax.pcolormesh(
                lon2d, lat2d, grid_2d,
                transform=ccrs.PlateCarree(),
                cmap=cmap, vmin=0.5, vmax=k + 0.5,
                shading="auto",
            )
            ax.coastlines(linewidth=0.5, color="black")
            ax.set_global()
            ax.set_title(title, fontsize=13, fontweight="bold")
            cbar = plt.colorbar(
                mesh, ax=ax, orientation="horizontal",
                pad=0.05, shrink=0.6, aspect=40,
            )
            cbar.set_label("Cluster", fontsize=11)
        else:
            fig, ax = plt.subplots(figsize=(14, 7))
            cmap = _make_cmap(k)
            mesh = ax.pcolormesh(
                lon2d, lat2d, grid_2d,
                cmap=cmap, vmin=0.5, vmax=k + 0.5,
                shading="auto",
            )
            ax.set_title(title, fontsize=13, fontweight="bold")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            plt.colorbar(mesh, ax=ax, label="Cluster")

        plt.tight_layout()
        path = os.path.join(out, filename)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close(fig)
        saved.append(path)
        if verbose:
            print(f"  ✓ {filename}")

    yr = f"{result['years'][0]}–{result['years'][-1]}"

    # Panel (a): EDC
    _plot_panel(
        result["labels_edc"], result["k_edc"],
        f"(a) EDC Clusters — AWI-ESM-1-1-LR Precipitation\n"
        f"k={result['k_edc']}  ({yr})  "
        f"ν={cfg.df}, α={cfg.alpha}, ε={cfg.neighbor_radius}",
        "fig9a_edc_clusters.png",
    )

    # Panel (b): LEC
    _plot_panel(
        result["labels_lec"], result["k_lec"],
        f"(b) LEC Clusters — AWI-ESM-1-1-LR Precipitation\n"
        f"k={result['k_lec']}  ({yr})  "
        f"ν={cfg.df}, α={cfg.alpha}, ε={cfg.neighbor_radius}",
        "fig9b_lec_clusters.png",
    )

    # Combined panel (a+b)
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 12),
        subplot_kw=({"projection": ccrs.Robinson()} if use_cartopy
                    else {}),
    )
    lon2d, lat2d = np.meshgrid(lons, lats)

    for ax, labels, k, panel in [
        (axes[0], result["labels_edc"], result["k_edc"], "a"),
        (axes[1], result["labels_lec"], result["k_lec"], "b"),
    ]:
        grid_2d = _labels_to_grid(labels)
        cmap = _make_cmap(k)
        if use_cartopy:
            mesh = ax.pcolormesh(
                lon2d, lat2d, grid_2d,
                transform=ccrs.PlateCarree(),
                cmap=cmap, vmin=0.5, vmax=k + 0.5,
                shading="auto",
            )
            ax.coastlines(linewidth=0.5)
            ax.set_global()
        else:
            mesh = ax.pcolormesh(
                lon2d, lat2d, grid_2d,
                cmap=cmap, vmin=0.5, vmax=k + 0.5,
                shading="auto",
            )
        label = "EDC" if panel == "a" else "LEC"
        ax.set_title(
            f"({panel}) {label} — k={k}", fontsize=13, fontweight="bold"
        )

    fig.suptitle(
        f"Figure 9 — Clustering of Precipitation Extremes\n"
        f"AWI-ESM-1-1-LR historical ({yr})  "
        f"ν={cfg.df}, α={cfg.alpha}, ε={cfg.neighbor_radius}",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    path = os.path.join(out, "fig9_combined.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    saved.append(path)
    if verbose:
        print(f"  ✓ fig9_combined.png")

    return saved

from __future__ import annotations

import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from weatherisk.cmip6_pipeline import (
    CMIP6Config,
    _compute_frechet_global,
    _detrend_grid_fast,
    _edc_matrix_flat,
    _grid_coords,
    _monthly_annual_maxima,
    _run_local_estimation_cmip6,
)


@dataclass
class HotPathBenchmarkConfig:
    seed: int = 12345
    n_years: int = 24
    n_lat: int = 8
    n_lon: int = 8
    n_workers: int = 4
    df: float = 5.0
    alpha: float = 1.0
    neighbor_radius: float = 3.0
    smoothing_radius: float = 2.0
    mle_ensemble: int = 3
    stl_period: int = 12


def _synthetic_monthly_precip(config: HotPathBenchmarkConfig) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(config.seed)
    n_months = config.n_years * 12
    times = np.arange(
        np.datetime64("1980-01"),
        np.datetime64("1980-01") + np.timedelta64(n_months, "M"),
        np.timedelta64(1, "M"),
    )

    lat = np.linspace(-1.5, 1.5, config.n_lat)[:, None]
    lon = np.linspace(-2.0, 2.0, config.n_lon)[None, :]
    spatial = 1.2 + 0.15 * np.cos(lat) + 0.12 * np.sin(lon)
    t = np.arange(n_months, dtype=float)
    seasonal = 0.35 * np.sin(2.0 * np.pi * t / 12.0)[:, None, None]
    trend = 0.0015 * t[:, None, None]
    noise = rng.gamma(shape=2.2, scale=0.35, size=(n_months, config.n_lat, config.n_lon))
    pr = spatial[None, :, :] + seasonal + trend + noise
    return pr.astype(float), times


def _git_revision() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def run_hotpath_benchmark(
    config: HotPathBenchmarkConfig | None = None,
    *,
    markdown_path: str | Path | None = None,
) -> dict[str, float | int | str | dict[str, float | int]]:
    config = config or HotPathBenchmarkConfig()
    pr, times = _synthetic_monthly_precip(config)
    pipeline_cfg = CMIP6Config(
        df=config.df,
        alpha=config.alpha,
        neighbor_radius=config.neighbor_radius,
        smoothing_radius=config.smoothing_radius,
        mle_ensemble=config.mle_ensemble,
        stl_period=config.stl_period,
        n_workers=config.n_workers,
    )

    timings: dict[str, float] = {}

    t0 = time.perf_counter()
    detrended = _detrend_grid_fast(pr, period=pipeline_cfg.stl_period, verbose=False)
    timings["step1a_detrend"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    annual_max, years = _monthly_annual_maxima(detrended, times, verbose=False)
    timings["step1b_annual_max"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    frechet, valid_idx = _compute_frechet_global(
        annual_max,
        n_workers=pipeline_cfg.n_workers,
        verbose=False,
    )
    timings["step2_gev_frechet"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    grid_coords = _grid_coords(valid_idx, annual_max.shape[1], annual_max.shape[2])
    est = _run_local_estimation_cmip6(
        frechet,
        grid_coords,
        pipeline_cfg,
        verbose=False,
    )
    timings["step3_local_mle"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    dm_edc = _edc_matrix_flat(frechet)
    timings["step5b_edc_matrix"] = time.perf_counter() - t0

    total = sum(timings.values())
    result = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "config": asdict(config),
        "derived": {
            "n_months": int(pr.shape[0]),
            "n_cells": int(config.n_lat * config.n_lon),
            "n_valid_cells": int(len(valid_idx)),
            "n_years_complete": int(len(years)),
        },
        "timings_seconds": timings,
        "total_seconds": total,
        "checks": {
            "frechet_min": float(np.min(frechet)),
            "frechet_max": float(np.max(frechet)),
            "edc_trace": float(np.trace(dm_edc)),
            "est_mean_a": float(np.mean(est[:, 0])),
            "est_mean_b": float(np.mean(est[:, 1])),
            "est_mean_gamma": float(np.mean(est[:, 2])),
        },
    }

    if markdown_path is not None:
        _append_benchmark_markdown(Path(markdown_path), result)

    return result


def _append_benchmark_markdown(path: Path, result: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(
            "# Benchmark Results\n\n"
            "This file records hot-path benchmark runs for the CMIP6 Figure 9 pipeline.\n"
            "Each run uses the real pipeline functions on deterministic synthetic data.\n\n",
            encoding="utf-8",
        )

    config = result["config"]
    derived = result["derived"]
    timings = result["timings_seconds"]
    checks = result["checks"]
    lines = [
        f"## Run {result['timestamp']}",
        "",
        f"- Git revision: `{result['git_revision']}`",
        f"- Total time: `{result['total_seconds']:.3f}s`",
        f"- Config: `{config}`",
        f"- Derived: `{derived}`",
        "",
        "| Step | Seconds |",
        "| --- | ---: |",
    ]
    for key, value in timings.items():
        lines.append(f"| {key} | {value:.3f} |")
    lines.extend(
        [
            "",
            f"- Checks: `{checks}`",
            "",
        ]
    )
    with path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines))

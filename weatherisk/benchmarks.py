from __future__ import annotations

import subprocess
import time
import os
import threading
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import weatherisk.cmip6_pipeline as cmip6_pipeline
from weatherisk.cmip6_pipeline import CMIP6Config, run_cmip6_pipeline


class _PeakMemorySampler:
    """Sample peak RSS across the current process tree during a benchmark."""

    def __init__(self, *, interval_seconds: float = 0.05):
        try:
            import psutil
        except ImportError as exc:
            raise RuntimeError(
                "psutil is required for benchmark memory profiling. "
                "Install weatherisk[dev] or pip install psutil."
            ) from exc

        self._psutil = psutil
        self._root = psutil.Process(os.getpid())
        self._interval_seconds = interval_seconds
        self._stop_event = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self.max_rss_bytes = 0

    @property
    def interval_seconds(self) -> float:
        return self._interval_seconds

    def start(self) -> None:
        self.sample_once()
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        self._thread.join(timeout=max(1.0, 5.0 * self._interval_seconds))
        self.sample_once()

    def sample_once(self) -> int:
        total_rss = 0
        procs = [self._root]
        try:
            procs.extend(self._root.children(recursive=True))
        except self._psutil.Error:
            pass

        for proc in procs:
            try:
                total_rss += proc.memory_info().rss
            except self._psutil.Error:
                continue

        if total_rss > self.max_rss_bytes:
            self.max_rss_bytes = total_rss
        return total_rss

    def _run(self) -> None:
        while not self._stop_event.wait(self._interval_seconds):
            self.sample_once()


@contextmanager
def _profile_peak_memory(*, interval_seconds: float = 0.05):
    sampler = _PeakMemorySampler(interval_seconds=interval_seconds)
    sampler.start()
    try:
        yield sampler
    finally:
        sampler.stop()


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


@contextmanager
def _patched_pipeline_environment(
    pr: np.ndarray,
    times: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
):
    """Patch the pipeline module so run_cmip6_pipeline uses synthetic input."""
    originals = {
        "load_monthly_precipitation": cmip6_pipeline.load_monthly_precipitation,
        "ensure_cmip6_data": cmip6_pipeline.ensure_cmip6_data,
    }

    def _fake_load_monthly_precipitation(*args, **kwargs):
        return pr.copy(), times.copy(), lats.copy(), lons.copy()

    def _fake_ensure_cmip6_data(*args, **kwargs):
        return "synthetic-benchmark"

    cmip6_pipeline.load_monthly_precipitation = _fake_load_monthly_precipitation
    cmip6_pipeline.ensure_cmip6_data = _fake_ensure_cmip6_data
    try:
        yield
    finally:
        cmip6_pipeline.load_monthly_precipitation = originals["load_monthly_precipitation"]
        cmip6_pipeline.ensure_cmip6_data = originals["ensure_cmip6_data"]


@contextmanager
def _timed_pipeline_steps(timings: dict[str, float]):
    """Wrap full-pipeline step functions so run_cmip6_pipeline records timings."""
    step_names = [
        "_detrend_grid_fast",
        "_monthly_annual_maxima",
        "_compute_frechet_global",
        "_run_local_estimation_cmip6",
        "_smooth_estimates_cmip6",
        "_run_clustering_cmip6",
        "_incluster_reestimate_cmip6",
    ]
    originals = {name: getattr(cmip6_pipeline, name) for name in step_names}

    def make_wrapper(name, func):
        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            out = func(*args, **kwargs)
            timings[name] = timings.get(name, 0.0) + (time.perf_counter() - t0)
            return out

        return wrapper

    for name, func in originals.items():
        setattr(cmip6_pipeline, name, make_wrapper(name, func))

    try:
        yield
    finally:
        for name, func in originals.items():
            setattr(cmip6_pipeline, name, func)


def run_hotpath_benchmark(
    config: HotPathBenchmarkConfig | None = None,
    *,
    markdown_path: str | Path | None = None,
) -> dict[str, float | int | str | dict[str, float | int | str]]:
    config = config or HotPathBenchmarkConfig()
    pr, times = _synthetic_monthly_precip(config)
    lats = np.linspace(-60.0, 60.0, config.n_lat)
    lons = np.linspace(0.0, 360.0 - 360.0 / config.n_lon, config.n_lon)
    pipeline_cfg = CMIP6Config(
        df=config.df,
        alpha=config.alpha,
        neighbor_radius=config.neighbor_radius,
        smoothing_radius=config.smoothing_radius,
        mle_ensemble=config.mle_ensemble,
        stl_period=config.stl_period,
        n_workers=config.n_workers,
        retain_clustering_artifacts=False,
    )

    timings: dict[str, float] = {}

    with TemporaryDirectory(prefix="weatherisk-benchmark-") as tmpdir:
        pipeline_cfg.output_dir = tmpdir
        with _profile_peak_memory() as memory_sampler:
            t0 = time.perf_counter()
            with _patched_pipeline_environment(pr, times, lats, lons), _timed_pipeline_steps(timings):
                result = run_cmip6_pipeline(pipeline_cfg, verbose=False)
            total = time.perf_counter() - t0

    result = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "entrypoint": "weatherisk.cmip6_pipeline.run_cmip6_pipeline",
        "config": asdict(config),
        "derived": {
            "n_months": int(pr.shape[0]),
            "n_cells": int(config.n_lat * config.n_lon),
            "n_valid_cells": int(len(result["valid_idx"])),
            "n_years_complete": int(len(result["years"])),
            "k_lec": int(result["k_lec"]),
            "k_edc": int(result["k_edc"]),
        },
        "timings_seconds": timings,
        "total_seconds": total,
        "memory_profile": {
            "metric": "peak_process_tree_rss",
            "sample_interval_seconds": memory_sampler.interval_seconds,
            "max_rss_bytes": int(memory_sampler.max_rss_bytes),
            "max_rss_gib": float(memory_sampler.max_rss_bytes / (1024 ** 3)),
        },
        "checks": {
            "frechet_min": float(np.min(result["frechet"])),
            "frechet_max": float(np.max(result["frechet"])),
            "labels_edc_sum": int(np.sum(result["labels_edc"])),
            "est_mean_a": float(np.mean(result["estimates"][:, 0])),
            "est_mean_b": float(np.mean(result["estimates"][:, 1])),
            "est_mean_gamma": float(np.mean(result["estimates"][:, 2])),
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
            "This file records benchmark runs for the CMIP6 Figure 9 pipeline.\n"
            "Each run uses the same `run_cmip6_pipeline` orchestration path used by the HPC script,\n"
            "with deterministic synthetic data injected only at the data-loading boundary.\n\n",
            encoding="utf-8",
        )

    config = result["config"]
    derived = result["derived"]
    timings = result["timings_seconds"]
    memory_profile = result["memory_profile"]
    checks = result["checks"]
    lines = [
        f"## Run {result['timestamp']}",
        "",
        f"- Git revision: `{result['git_revision']}`",
        f"- Entrypoint: `{result['entrypoint']}`",
        f"- Total time: `{result['total_seconds']:.3f}s`",
        f"- Config: `{config}`",
        f"- Derived: `{derived}`",
        f"- Max memory: `{memory_profile['max_rss_gib']:.3f} GiB` "
        f"(`{memory_profile['max_rss_bytes']}` bytes; "
        f"{memory_profile['metric']}, Δt={memory_profile['sample_interval_seconds']}s)",
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

from __future__ import annotations

import io
import importlib.util
import json
import os
import subprocess
import threading
import time
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import weatherisk.cmip6_pipeline as cmip6_pipeline
from weatherisk.cmip6_pipeline import CMIP6Config


_REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_reproduce_fig9_module():
    module_path = _REPO_ROOT / "scripts" / "reproduce_fig9.py"
    spec = importlib.util.spec_from_file_location("weatherisk_reproduce_fig9", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load benchmark entrypoint from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
class Figure9BenchmarkConfig:
    seed: int = 12345
    n_years: int = 12
    n_lat: int = 5
    n_lon: int = 5
    n_workers: int = 4
    year_start: int = 1980
    dpi: int = 300
    generate_plots: bool = True
    backend: str = "python"
    suppress_script_output: bool = True


@dataclass
class DecisionBenchmarkConfig:
    benchmark: Figure9BenchmarkConfig = field(
        default_factory=lambda: Figure9BenchmarkConfig(
            n_years=12,
            n_lat=5,
            n_lon=5,
            n_workers=4,
        )
    )
    warmup_runs: int = 0
    measured_runs: int = 1


HotPathBenchmarkConfig = Figure9BenchmarkConfig


def _synthetic_monthly_precip(config: Figure9BenchmarkConfig) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(config.seed)
    n_months = config.n_years * 12
    times = np.arange(
        np.datetime64(f"{config.year_start:04d}-01"),
        np.datetime64(f"{config.year_start:04d}-01") + np.timedelta64(n_months, "M"),
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
def _forced_backend(backend: str):
    old_backend = os.environ.get("WEATHERISK_BACKEND")
    os.environ["WEATHERISK_BACKEND"] = backend
    try:
        yield
    finally:
        if old_backend is None:
            os.environ.pop("WEATHERISK_BACKEND", None)
        else:
            os.environ["WEATHERISK_BACKEND"] = old_backend


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
        "plot_figure9",
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


def run_figure9_benchmark(
    config: Figure9BenchmarkConfig | None = None,
    *,
    markdown_path: str | Path | None = None,
) -> dict[str, float | int | str | dict[str, float | int | str]]:
    config = config or Figure9BenchmarkConfig()
    pr, times = _synthetic_monthly_precip(config)
    lats = np.linspace(-60.0, 60.0, config.n_lat)
    lons = np.linspace(0.0, 360.0 - 360.0 / config.n_lon, config.n_lon)
    timings: dict[str, float] = {}

    with TemporaryDirectory(prefix="weatherisk-benchmark-") as tmpdir:
        with _profile_peak_memory() as memory_sampler:
            t0 = time.perf_counter()
            with _forced_backend(config.backend), _patched_pipeline_environment(pr, times, lats, lons), _timed_pipeline_steps(timings):
                reproduce_fig9 = _load_reproduce_fig9_module()

                argv = [
                    "--data-dir", "synthetic-benchmark",
                    "--workers", str(config.n_workers),
                    "--output-dir", tmpdir,
                    "--year-start", str(config.year_start),
                    "--year-end", str(config.year_start + config.n_years - 1),
                    "--dpi", str(config.dpi),
                ]
                if not config.generate_plots:
                    argv.append("--skip-plots")

                if config.suppress_script_output:
                    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                        run_out = reproduce_fig9.main(argv)
                else:
                    run_out = reproduce_fig9.main(argv)
            total = time.perf_counter() - t0

        result = run_out["result"]
        saved_figures = run_out["saved_figures"]

    result = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "script_entrypoint": "scripts.reproduce_fig9.main",
        "pipeline_entrypoint": "weatherisk.cmip6_pipeline.run_cmip6_pipeline",
        "config": asdict(config),
        "derived": {
            "n_months": int(pr.shape[0]),
            "n_cells": int(config.n_lat * config.n_lon),
            "n_valid_cells": int(len(result["valid_idx"])),
            "n_years_complete": int(len(result["years"])),
            "k_lec": int(result["k_lec"]),
            "k_edc": int(result["k_edc"]),
            "saved_figure_count": int(len(saved_figures)),
        },
        "invocation": {
            "backend": config.backend,
            "plots_enabled": bool(config.generate_plots),
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


run_hotpath_benchmark = run_figure9_benchmark


def run_decision_benchmark(
    config: DecisionBenchmarkConfig | None = None,
    *,
    markdown_path: str | Path | None = None,
) -> dict[str, object]:
    config = config or DecisionBenchmarkConfig()

    warmups = []
    for _ in range(config.warmup_runs):
        warmups.append(run_figure9_benchmark(config.benchmark, markdown_path=None))

    measured = []
    for _ in range(config.measured_runs):
        measured.append(run_figure9_benchmark(config.benchmark, markdown_path=None))

    total_seconds = np.array([run["total_seconds"] for run in measured], dtype=float)
    peak_rss_bytes = np.array(
        [run["memory_profile"]["max_rss_bytes"] for run in measured], dtype=np.int64
    )
    step_names = list(measured[0]["timings_seconds"].keys())
    step_stats: dict[str, dict[str, float]] = {}
    for step in step_names:
        values = np.array([run["timings_seconds"][step] for run in measured], dtype=float)
        step_stats[step] = {
            "mean_seconds": float(values.mean()),
            "min_seconds": float(values.min()),
            "max_seconds": float(values.max()),
            "std_seconds": float(values.std(ddof=0)),
        }

    checks_reference = measured[0]["checks"]
    checks_stable = all(run["checks"] == checks_reference for run in measured[1:])

    result: dict[str, object] = {
        "timestamp": measured[-1]["timestamp"],
        "git_revision": measured[-1]["git_revision"],
        "script_entrypoint": measured[-1]["script_entrypoint"],
        "pipeline_entrypoint": measured[-1]["pipeline_entrypoint"],
        "benchmark_method": {
            "case": "decision",
            "warmup_runs": config.warmup_runs,
            "measured_runs": config.measured_runs,
            "summary": "report mean/min/max/std over measured runs; warmups excluded",
        },
        "config": asdict(config.benchmark),
        "derived": measured[-1]["derived"],
        "warmup_runs": warmups,
        "measured_runs": measured,
        "summary": {
            "total_seconds": {
                "mean_seconds": float(total_seconds.mean()),
                "min_seconds": float(total_seconds.min()),
                "max_seconds": float(total_seconds.max()),
                "std_seconds": float(total_seconds.std(ddof=0)),
            },
            "memory_profile": {
                "metric": measured[-1]["memory_profile"]["metric"],
                "sample_interval_seconds": measured[-1]["memory_profile"]["sample_interval_seconds"],
                "mean_max_rss_bytes": float(peak_rss_bytes.mean()),
                "max_max_rss_bytes": int(peak_rss_bytes.max()),
                "mean_max_rss_gib": float(peak_rss_bytes.mean() / (1024 ** 3)),
                "max_max_rss_gib": float(peak_rss_bytes.max() / (1024 ** 3)),
            },
            "step_timings": step_stats,
            "checks_stable": checks_stable,
            "checks_reference": checks_reference,
        },
    }

    if markdown_path is not None:
        _append_decision_benchmark_markdown(Path(markdown_path), result)

    return result


def _append_benchmark_markdown(path: Path, result: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(
            "# Benchmark Results\n\n"
            "This file is the canonical benchmark log for the CMIP6 Figure 9 path.\n"
            "Each benchmark runs the real `scripts/reproduce_fig9.py` entrypoint used by the SLURM job,\n"
            "with deterministic synthetic data injected only at the data-loading boundary so the reduced case\n"
            "finishes quickly while preserving the production call chain.\n\n"
            "## Benchmark Workflow\n\n"
            "- Canonical runner: `python scripts/run_fig9_benchmark.py`\n"
            "- Script entrypoint under test: `scripts.reproduce_fig9.main`\n"
            "- Pipeline entrypoint under test: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`\n"
            "- Default benchmark shape: reduced synthetic Figure 9 case intended to finish in about a minute or less\n"
            "- Measured metrics: total wall time, per-step timings, and peak process-tree RSS\n\n"
            "## Key Learnings\n\n"
            "- A real script entrypoint is required for multiprocessing benchmarks; heredoc and stdin-driven entrypoints are not reliable on macOS spawn mode.\n"
            "- CMIP6 local MLE remains the dominant hotspot; pair-array assembly matters less than optimizer cost once the production objective is in play.\n"
            "- Memory wins in clustering matter most on larger grids; small synthetic cases can hide them behind interpreter and worker overhead.\n\n",
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
        f"- Script entrypoint: `{result['script_entrypoint']}`",
        f"- Pipeline entrypoint: `{result['pipeline_entrypoint']}`",
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


def _append_decision_benchmark_markdown(path: Path, result: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(
            "# Benchmark Results\n\n"
            "This file is the canonical benchmark log for the CMIP6 Figure 9 path.\n"
            "Each benchmark runs the real `scripts/reproduce_fig9.py` entrypoint used by the SLURM job,\n"
            "with deterministic synthetic data injected only at the data-loading boundary so the reduced case\n"
            "finishes quickly while preserving the production call chain.\n\n",
            encoding="utf-8",
        )

    summary = result["summary"]
    total = summary["total_seconds"]
    memory = summary["memory_profile"]
    step_timings = summary["step_timings"]
    lines = [
        f"## Decision Benchmark {result['timestamp']}",
        "",
        f"- Git revision: `{result['git_revision']}`",
        f"- Script entrypoint: `{result['script_entrypoint']}`",
        f"- Pipeline entrypoint: `{result['pipeline_entrypoint']}`",
        f"- Method: `{result['benchmark_method']['warmup_runs']} warmup + {result['benchmark_method']['measured_runs']} measured runs` (warmups excluded from summary)",
        f"- Benchmark case: `reduced Figure 9`",
        f"- Config: `{result['config']}`",
        f"- Derived: `{result['derived']}`",
        f"- Total time summary: mean `{total['mean_seconds']:.3f}s`, min `{total['min_seconds']:.3f}s`, max `{total['max_seconds']:.3f}s`, std `{total['std_seconds']:.3f}s`",
        f"- Peak memory summary: mean `{memory['mean_max_rss_gib']:.3f} GiB`, max `{memory['max_max_rss_gib']:.3f} GiB` (`{memory['metric']}`, Δt={memory['sample_interval_seconds']}s)",
        f"- Checks stable across measured runs: `{summary['checks_stable']}`",
        "",
        "| Step | Mean (s) | Min (s) | Max (s) | Std (s) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for step, stats in step_timings.items():
        lines.append(
            f"| {step} | {stats['mean_seconds']:.3f} | {stats['min_seconds']:.3f} | {stats['max_seconds']:.3f} | {stats['std_seconds']:.3f} |"
        )
    lines.extend(
        [
            "",
            f"- Reference checks: `{summary['checks_reference']}`",
            "",
        ]
    )
    with path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines))

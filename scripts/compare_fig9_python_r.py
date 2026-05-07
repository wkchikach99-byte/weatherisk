#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
import tempfile
import time
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np

import weatherisk.cmip6_pipeline as cm
from weatherisk.backend import calc_distance_ellipses_condensed
from weatherisk.benchmarks import (
    Figure9BenchmarkConfig,
    _patched_pipeline_environment,
    _synthetic_monthly_precip,
    run_figure9_benchmark,
)
from weatherisk.cmip6_pipeline import CMIP6Config


def _write_monthly_input(config: Figure9BenchmarkConfig, out_dir: Path) -> None:
    pr, _times = _synthetic_monthly_precip(config)
    n_months, n_lat, n_lon = pr.shape

    monthly_pr_path = out_dir / "monthly_pr.csv"
    with monthly_pr_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([f"cell_{idx}" for idx in range(n_lat * n_lon)])
        flat = pr.reshape(n_months, n_lat * n_lon, order="F")
        writer.writerows(flat)

    scalar_params_path = out_dir / "scalar_params.csv"
    with scalar_params_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "value"])
        writer.writerows(
            [
                ["resolution", n_lat],
                ["n_years", config.n_years],
                ["n_months", n_months],
                ["df", 5.0],
                ["alpha", 1.0],
                ["neighbor_radius", 5.0],
                ["smoothing_radius", 2.0],
                ["mle_ensemble", 5],
            ]
        )


def _python_thresholds(config: Figure9BenchmarkConfig) -> dict[str, float | int]:
    pr, times = _synthetic_monthly_precip(config)
    lats = np.linspace(-60.0, 60.0, config.n_lat)
    lons = np.linspace(0.0, 360.0 - 360.0 / config.n_lon, config.n_lon)
    pipe_cfg = CMIP6Config(
        output_dir="output/_tmp_py_compare",
        year_start=config.year_start,
        year_end=config.year_start + config.n_years - 1,
        n_workers=config.n_workers,
        retain_clustering_artifacts=False,
    )
    with _patched_pipeline_environment(pr, times, lats, lons):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            result = cm.run_cmip6_pipeline(pipe_cfg, verbose=False)

    lec_condensed = calc_distance_ellipses_condensed(
        result["smoothed"],
        res=21,
        chunk_size=pipe_cfg.lec_chunk_size,
    )
    edc_condensed = cm._edc_condensed_flat(result["frechet"])
    return {
        "k_lec": int(result["k_lec"]),
        "k_edc": int(result["k_edc"]),
        "q30_lec": float(np.quantile(lec_condensed, pipe_cfg.quantile_threshold)),
        "q30_edc": float(np.quantile(edc_condensed, pipe_cfg.quantile_threshold)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare Python and R Figure 9-style k values on the same benchmark input.")
    parser.add_argument("--years", type=int, default=12)
    parser.add_argument("--lat", type=int, default=6)
    parser.add_argument("--lon", type=int, default=6)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--skip-python-benchmark", action="store_true")
    args = parser.parse_args()

    config = Figure9BenchmarkConfig(
        n_years=args.years,
        n_lat=args.lat,
        n_lon=args.lon,
        n_workers=args.workers,
        generate_plots=False,
    )

    python_result = None
    python_thresholds = _python_thresholds(config)
    if not args.skip_python_benchmark:
        python_result = run_figure9_benchmark(config, markdown_path=None)

    with tempfile.TemporaryDirectory(prefix="weatherisk-r-compare-") as tmpdir:
        tmp_path = Path(tmpdir)
        input_dir = tmp_path / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        _write_monthly_input(config, input_dir)

        output_json = tmp_path / "r_result.json"
        t0_r = time.perf_counter()
        subprocess.run(
            [
                "/usr/local/bin/Rscript",
                "scripts/replay_fig9_benchmark_in_r.R",
                str(input_dir),
                str(output_json),
            ],
            check=True,
            cwd=Path(__file__).resolve().parents[1],
        )
        r_total_seconds = time.perf_counter() - t0_r
        r_result = json.loads(output_json.read_text(encoding="utf-8"))

    runtime_compare = {
        "python_total_seconds": None if python_result is None else float(python_result["total_seconds"]),
        "r_total_seconds": float(r_total_seconds),
        "r_div_python": None,
        "python_div_r": None,
    }
    if python_result is not None and python_result["total_seconds"] > 0:
        runtime_compare["r_div_python"] = float(r_total_seconds / python_result["total_seconds"])
        runtime_compare["python_div_r"] = float(python_result["total_seconds"] / r_total_seconds)

    summary = {
        "config": {
            "n_years": config.n_years,
            "n_lat": config.n_lat,
            "n_lon": config.n_lon,
            "n_workers": config.n_workers,
        },
        "python": None if python_result is None else {
            "total_seconds": python_result["total_seconds"],
            **python_thresholds,
        },
        "python_thresholds": python_thresholds if python_result is None else None,
        "r": {
            "total_seconds": float(r_total_seconds),
            **r_result,
        },
        "runtime_compare": runtime_compare,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
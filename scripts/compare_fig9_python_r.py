#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path

from weatherisk.benchmarks import Figure9BenchmarkConfig, _synthetic_monthly_precip, run_figure9_benchmark


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
    if not args.skip_python_benchmark:
        python_result = run_figure9_benchmark(config, markdown_path=None)

    with tempfile.TemporaryDirectory(prefix="weatherisk-r-compare-") as tmpdir:
        tmp_path = Path(tmpdir)
        input_dir = tmp_path / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        _write_monthly_input(config, input_dir)

        output_json = tmp_path / "r_result.json"
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
        r_result = json.loads(output_json.read_text(encoding="utf-8"))

    summary = {
        "config": {
            "n_years": config.n_years,
            "n_lat": config.n_lat,
            "n_lon": config.n_lon,
            "n_workers": config.n_workers,
        },
        "python": None if python_result is None else {
            "k_lec": python_result["derived"]["k_lec"],
            "k_edc": python_result["derived"]["k_edc"],
            "total_seconds": python_result["total_seconds"],
        },
        "r": r_result,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
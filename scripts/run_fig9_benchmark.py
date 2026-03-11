#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from weatherisk.benchmarks import (
    DecisionBenchmarkConfig,
    Figure9BenchmarkConfig,
    run_decision_benchmark,
    run_figure9_benchmark,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the canonical reduced Figure 9 benchmark."
    )
    parser.add_argument("--years", type=int, default=12)
    parser.add_argument("--lat", type=int, default=5)
    parser.add_argument("--lon", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--year-start", type=int, default=1980)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--backend",
        choices=["python"],
        default="python",
        help="Backend to use during the benchmark. Only the Python backend is active.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation if you only want the pipeline timing.",
    )
    parser.add_argument("--warmup-runs", type=int, default=0)
    parser.add_argument("--measured-runs", type=int, default=1)
    parser.add_argument("--markdown-path", default="docs/benchmark_results.md")
    args = parser.parse_args()

    benchmark = Figure9BenchmarkConfig(
        n_years=args.years,
        n_lat=args.lat,
        n_lon=args.lon,
        n_workers=args.workers,
        year_start=args.year_start,
        dpi=args.dpi,
        generate_plots=not args.skip_plots,
        backend=args.backend,
    )

    if args.warmup_runs == 0 and args.measured_runs == 1:
        result = run_figure9_benchmark(benchmark, markdown_path=args.markdown_path)
        print(f"Total: {result['total_seconds']:.3f}s")
        print(
            f"Peak RSS: {result['memory_profile']['max_rss_gib']:.3f} GiB "
            f"({result['memory_profile']['max_rss_bytes']} bytes)"
        )
        print(
            f"k_LEC={result['derived']['k_lec']}  "
            f"k_EDC={result['derived']['k_edc']}  "
            f"figures={result['derived']['saved_figure_count']}"
        )
        print(json.dumps(result, indent=2))
        return

    result = run_decision_benchmark(
        DecisionBenchmarkConfig(
            benchmark=benchmark,
            warmup_runs=args.warmup_runs,
            measured_runs=args.measured_runs,
        ),
        markdown_path=args.markdown_path,
    )
    summary = result["summary"]
    total = summary["total_seconds"]
    mem = summary["memory_profile"]
    print(
        f"Total: {total['mean_seconds']:.3f}s "
        f"(min={total['min_seconds']:.3f}, max={total['max_seconds']:.3f}, std={total['std_seconds']:.3f})"
    )
    print(
        f"Peak RSS: {mem['mean_max_rss_gib']:.3f} GiB "
        f"(max={mem['max_max_rss_gib']:.3f})"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
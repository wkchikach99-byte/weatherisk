#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from weatherisk.benchmarks import (
    Figure9BenchmarkConfig,
    run_figure9_benchmark,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the canonical reduced Figure 9 benchmark."
    )
    parser.add_argument("--years", type=int, default=12)
    parser.add_argument("--lat", type=int, default=6)
    parser.add_argument("--lon", type=int, default=6)
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

    result = run_figure9_benchmark(benchmark, markdown_path=args.markdown_path)
    mem = result["memory_profile"]["max_rss"]
    print(
        f"Total: {result['total_seconds']:.3f}s"
    )
    print(
        f"Peak RSS: {mem['bytes']} bytes "
        f"({mem['kib']:.1f} KiB, {mem['mib']:.3f} MiB, {mem['gib']:.3f} GiB)"
    )
    print(
        f"k_LEC={result['derived']['k_lec']}  "
        f"k_EDC={result['derived']['k_edc']}  "
        f"figures={result['derived']['saved_figure_count']}"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
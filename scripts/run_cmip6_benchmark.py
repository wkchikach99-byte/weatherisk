#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from weatherisk.benchmarks import HotPathBenchmarkConfig, run_hotpath_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CMIP6 pipeline benchmark harness")
    parser.add_argument("--years", type=int, default=36)
    parser.add_argument("--lat", type=int, default=12)
    parser.add_argument("--lon", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--neighbor-radius", type=float, default=3.0)
    parser.add_argument("--smoothing-radius", type=float, default=2.0)
    parser.add_argument("--mle-ensemble", type=int, default=3)
    parser.add_argument("--markdown-path", default="docs/benchmark_results.md")
    args = parser.parse_args()

    result = run_hotpath_benchmark(
        HotPathBenchmarkConfig(
            n_years=args.years,
            n_lat=args.lat,
            n_lon=args.lon,
            n_workers=args.workers,
            neighbor_radius=args.neighbor_radius,
            smoothing_radius=args.smoothing_radius,
            mle_ensemble=args.mle_ensemble,
        ),
        markdown_path=args.markdown_path,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
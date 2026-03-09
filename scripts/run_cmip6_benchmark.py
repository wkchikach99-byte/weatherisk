#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from weatherisk.benchmarks import (
    DecisionBenchmarkConfig,
    HotPathBenchmarkConfig,
    run_decision_benchmark,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CMIP6 pipeline decision benchmark harness")
    parser.add_argument("--years", type=int, default=48)
    parser.add_argument("--lat", type=int, default=16)
    parser.add_argument("--lon", type=int, default=16)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--neighbor-radius", type=float, default=3.0)
    parser.add_argument("--smoothing-radius", type=float, default=2.0)
    parser.add_argument("--mle-ensemble", type=int, default=3)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=5)
    parser.add_argument("--markdown-path", default="docs/benchmark_results.md")
    args = parser.parse_args()

    result = run_decision_benchmark(
        DecisionBenchmarkConfig(
            benchmark=HotPathBenchmarkConfig(
                n_years=args.years,
                n_lat=args.lat,
                n_lon=args.lon,
                n_workers=args.workers,
                neighbor_radius=args.neighbor_radius,
                smoothing_radius=args.smoothing_radius,
                mle_ensemble=args.mle_ensemble,
            ),
            warmup_runs=args.warmup_runs,
            measured_runs=args.measured_runs,
        ),
        markdown_path=args.markdown_path,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
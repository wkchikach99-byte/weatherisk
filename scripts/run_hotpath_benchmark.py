from __future__ import annotations

import argparse
import json

from weatherisk.benchmarks import HotPathBenchmarkConfig, run_hotpath_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CMIP6 hot-path benchmark.")
    parser.add_argument(
        "--markdown-path",
        default="docs/benchmark_results.md",
        help="Markdown file where benchmark runs are appended.",
    )
    parser.add_argument("--years", type=int, default=24)
    parser.add_argument("--lat", type=int, default=8)
    parser.add_argument("--lon", type=int, default=8)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    config = HotPathBenchmarkConfig(
        n_years=args.years,
        n_lat=args.lat,
        n_lon=args.lon,
        n_workers=args.workers,
    )
    result = run_hotpath_benchmark(config, markdown_path=args.markdown_path)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
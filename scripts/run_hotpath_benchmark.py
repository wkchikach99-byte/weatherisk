from __future__ import annotations

import argparse
import json
import os

from weatherisk.benchmarks import DecisionBenchmarkConfig, HotPathBenchmarkConfig, run_decision_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CMIP6 decision benchmark.")
    parser.add_argument(
        "--markdown-path",
        default=None,
        help="Markdown file where benchmark runs are appended (default: none).",
    )
    parser.add_argument("--years", type=int, default=48)
    parser.add_argument("--lat", type=int, default=16)
    parser.add_argument("--lon", type=int, default=16)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--backend",
        choices=["auto", "python", "rust"],
        default="auto",
        help="Force a specific backend (sets WEATHERISK_BACKEND).",
    )
    args = parser.parse_args()

    if args.backend != "auto":
        os.environ["WEATHERISK_BACKEND"] = args.backend

    config = DecisionBenchmarkConfig(
        benchmark=HotPathBenchmarkConfig(
            n_years=args.years,
            n_lat=args.lat,
            n_lon=args.lon,
            n_workers=args.workers,
        ),
    )

    from weatherisk.backend import _USE_RUST
    backend_label = "rust" if _USE_RUST else "python"
    print(f"Backend: {backend_label}")
    print(f"Protocol: {config.warmup_runs} warmup + {config.measured_runs} measured runs")
    print(f"Config: {args.years}yr {args.lat}x{args.lon} {args.workers} workers")
    print()

    result = run_decision_benchmark(config, markdown_path=args.markdown_path)

    summary = result["summary"]
    total = summary["total_seconds"]
    mem = summary["memory_profile"]
    steps = summary["step_timings"]
    checks_stable = summary["checks_stable"]

    print(f"Total: {total['mean_seconds']:.3f}s (min={total['min_seconds']:.3f}, max={total['max_seconds']:.3f}, std={total['std_seconds']:.3f})")
    print(f"Peak RSS: {mem['mean_max_rss_gib']:.3f} GiB (max={mem['max_max_rss_gib']:.3f})")
    print(f"Checks stable: {checks_stable}")
    print()
    print(f"{'Step':<35} {'Mean':>8} {'Min':>8} {'Max':>8} {'Std':>8}")
    print("-" * 75)
    for step, stats in steps.items():
        print(f"{step:<35} {stats['mean_seconds']:>8.3f} {stats['min_seconds']:>8.3f} {stats['max_seconds']:>8.3f} {stats['std_seconds']:>8.3f}")

    print()
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
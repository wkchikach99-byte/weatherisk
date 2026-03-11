#!/usr/bin/env python3
"""Reproduce Figure 9 from Contzen et al. (2025, Extremes 28:713–737).

Generates the global LEC and EDC cluster maps applied to monthly
precipitation data from AWI-ESM-1-1-LR (historical, 1850–2005).

Paper parameters:
    ν = 5, α = 1, ε = 5 (grid point distance)
    30%-quantile threshold → k_EDC = 104, k_LEC = 24

Usage:
    # Quick test (downloads data if needed, serial):
    python scripts/reproduce_fig9.py

    # Full run with 16 parallel workers:
    python scripts/reproduce_fig9.py --workers 16

    # Use pre-existing data directory:
    python scripts/reproduce_fig9.py --data-dir /pool/data/.../pr/gn/

    # Custom output:
    python scripts/reproduce_fig9.py --output-dir output/fig9 --dpi 600
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(
        description="Reproduce Figure 9 from the Extremes paper"
    )
    parser.add_argument(
        "--data-dir", default=None,
        help="Directory containing AWI-ESM-1-1-LR pr NetCDF files. "
             "Auto-discovered/downloaded if not provided.",
    )
    parser.add_argument(
        "--output-dir", default="output/cmip6_fig9",
        help="Output directory for figures and intermediate data.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of parallel workers (default: 1; use 16+ on HPC).",
    )
    parser.add_argument(
        "--year-start", type=int, default=1850,
        help="First year (inclusive).",
    )
    parser.add_argument(
        "--year-end", type=int, default=2005,
        help="Last year (inclusive).",
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="Figure DPI.",
    )
    parser.add_argument(
        "--no-download", action="store_true",
        help="Don't attempt to download data if missing.",
    )
    parser.add_argument(
        "--skip-plots", action="store_true",
        help="Run pipeline only, skip figure generation.",
    )
    args = parser.parse_args(argv)

    from weatherisk.cmip6_pipeline import CMIP6Config, run_cmip6_pipeline, plot_figure9
    from weatherisk.cmip6_data import DEFAULT_DATA_DIR

    cfg = CMIP6Config(
        data_dir=args.data_dir or DEFAULT_DATA_DIR,
        output_dir=args.output_dir,
        year_start=args.year_start,
        year_end=args.year_end,
        n_workers=args.workers,
    )

    os.makedirs(cfg.output_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print("  Reproducing Figure 9 — Contzen et al. (2025, Extremes)")
    print("  AWI-ESM-1-1-LR historical precipitation")
    print(f"  ν={cfg.df}, α={cfg.alpha}, ε={cfg.neighbor_radius}")
    print(f"  Period: {cfg.year_start}–{cfg.year_end}")
    print(f"  Workers: {cfg.n_workers}")
    from weatherisk.backend import _USE_RUST
    print(f"  Backend: {'Rust (weatherisk_core)' if _USE_RUST else 'Python (pure NumPy/SciPy)'}")
    print(f"  Output: {cfg.output_dir}")
    print("=" * 60)

    # Run the pipeline
    result = run_cmip6_pipeline(cfg, verbose=True)

    print(f"\n  ━━━ Results ━━━")
    print(f"  k_LEC = {result['k_lec']}  (paper: 24)")
    print(f"  k_EDC = {result['k_edc']}  (paper: 104)")

    # Generate plots
    saved: list[str] = []
    if not args.skip_plots:
        print("\n  Generating Figure 9 …")
        saved = plot_figure9(result, dpi=args.dpi, verbose=True)
        print(f"\n  {len(saved)} figure(s) saved to {cfg.output_dir}")
    else:
        print("\n  Plot generation skipped (--skip-plots).")

    print("\nDone.")
    return {
        "result": result,
        "saved_figures": saved,
        "config": cfg,
    }


if __name__ == "__main__":
    main()

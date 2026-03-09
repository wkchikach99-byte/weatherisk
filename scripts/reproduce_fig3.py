#!/usr/bin/env python3
"""Reproduce Fig. 3 from Contzen et al. (2025, Extremes 28:713-737).

Usage:
    # Sanity check (11x11, 50 obs, ~2 min):
    python scripts/reproduce_fig3.py --resolution 11 --n_sim 50

    # Full resolution with 8 parallel workers:
    python scripts/reproduce_fig3.py --workers 8
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.pipeline import run_nonstationary_pipeline
from weatherisk.plotting import plot_cluster_comparison


def main():
    parser = argparse.ArgumentParser(description="Reproduce Fig. 3")
    parser.add_argument("--resolution", type=int, default=None,
                        help="Grid resolution (default: preset value, 51)")
    parser.add_argument("--n_sim", type=int, default=None,
                        help="Number of simulations (default: preset value, 250)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default="docs/figures/fig3_paper_stripes.png")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel workers for local MLE (default: 1)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    result = run_nonstationary_pipeline(
        "paper_stripes",
        resolution=args.resolution,
        n_sim=args.n_sim,
        seed=args.seed,
        n_workers=args.workers,
    )

    plot_cluster_comparison(
        result["grid"], result["labels_edc"], result["labels_lec"],
        result["inclusters_edc"], result["inclusters_lec"],
        result["b_matrix"],
        param_index=1, param_name="b",
        vmin=0, vmax=5,
        label="$b$",
        suptitle=f"Fig. 3 — EDC vs LEC (k_edc={result['k_edc']}, k_lec={result['k_lec']})",
        filename=args.output,
    )

    print(f"\nk_lec = {result['k_lec']}")
    print(f"k_edc = {result['k_edc']}")
    print("LEC in-cluster estimates (a, b, gamma):")
    for cl in range(result["inclusters_lec"].shape[0]):
        row = result["inclusters_lec"][cl]
        if row[3] > 0:
            print(f"  cluster {cl}: a={row[0]:.2f}, b={row[1]:.2f}, g={row[2]:.3f}, n={int(row[3])}")
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()

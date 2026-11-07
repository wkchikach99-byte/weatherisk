#!/usr/bin/env python3
"""Recovery metrics for the stripes experiment (paper Figure 2).

The true range field b(x) = (x + 5) / 2 is continuous, so there is no true
partition and the adjusted Rand index is undefined. Two metrics replace it:

* between-cluster share of the variance of true b over cells, 1 - SS_within / SS_total;
* RMSE over cells of the in-cluster estimate b_hat(cluster(cell)) against true b(cell),
  restricted to cells whose cluster has an in-cluster estimate (at least 5 cells).

Usage:
    python scripts/stripes_metrics.py OUTPUT_DIR [OUTPUT_DIR ...]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.grid import Grid
from weatherisk.parameters import get_preset


def stripe_metrics(labels: np.ndarray, inclusters: np.ndarray, b_true: np.ndarray) -> dict:
    total = np.sum((b_true - b_true.mean()) ** 2)
    within = sum(
        np.sum((b_true[labels == c] - b_true[labels == c].mean()) ** 2)
        for c in np.unique(labels)
    )
    b_hat = inclusters[labels, 1]
    covered = np.isfinite(b_hat)
    rmse = float(np.sqrt(np.mean((b_hat[covered] - b_true[covered]) ** 2)))
    return {
        "k": int(np.unique(labels).size),
        "between_share": float(1.0 - within / total),
        "rmse_b": rmse,
        "covered_cells": int(covered.sum()),
        "cells": int(labels.size),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dirs", nargs="+")
    args = parser.parse_args()

    preset = get_preset("paper_stripes")
    print("dir\tmethod\tk\tbetween_share\trmse_b\tcovered_cells\tcells")
    for d in args.output_dirs:
        labels_lec = np.load(os.path.join(d, "labels_lec.npy"))
        res = int(round(np.sqrt(labels_lec.size)))
        grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=res)
        b_true = np.asarray(preset.b_func(grid.X_flat, grid.Y_flat), dtype=float)
        for method in ("edc", "lec"):
            labels = np.load(os.path.join(d, f"labels_{method}.npy")).astype(int)
            inclusters = np.load(os.path.join(d, f"inclusters_{method}.npy"))
            m = stripe_metrics(labels, inclusters, b_true)
            print(f"{d}\t{method}\t{m['k']}\t{m['between_share']:.4f}\t{m['rmse_b']:.4f}"
                  f"\t{m['covered_cells']}\t{m['cells']}")


if __name__ == "__main__":
    main()

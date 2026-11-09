#!/usr/bin/env python3
"""Two-by-two Figure 9 panel: EDC and LEC clusters with and without STL detrending.

Usage:
    python scripts/plot_fig9_comparison.py NODETREND.npz DETREND.npz OUTPUT.png

Each input is a ``pipeline_results.npz`` written by ``scripts/reproduce_fig9.py``.
"""

import argparse

import cartopy.crs as ccrs
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402


def cluster_cmap(k: int) -> ListedColormap:
    base = (list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
            + list(plt.get_cmap("tab20b").colors) + list(plt.get_cmap("tab20c").colors))
    return ListedColormap([base[i % len(base)] for i in range(k)])


def label_grid(res, key: str) -> tuple[np.ndarray, int]:
    labels = res[key]
    k = int(np.unique(labels).size)
    grid = np.full(res["lats"].size * res["lons"].size, np.nan)
    grid[res["valid_idx"]] = np.unique(labels, return_inverse=True)[1] + 1
    return grid.reshape(res["lats"].size, res["lons"].size), k


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("nodetrend")
    ap.add_argument("detrend")
    ap.add_argument("output")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    rows = [("no detrending", np.load(args.nodetrend)), ("STL detrending", np.load(args.detrend))]
    cols = [("EDC", "labels_edc"), ("LEC", "labels_lec")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.4),
                             subplot_kw={"projection": ccrs.Robinson()})
    for r, (row_name, res) in enumerate(rows):
        lon2d, lat2d = np.meshgrid(res["lons"], res["lats"])
        for c, (method, key) in enumerate(cols):
            grid, k = label_grid(res, key)
            ax = axes[r, c]
            ax.pcolormesh(lon2d, lat2d, grid, transform=ccrs.PlateCarree(),
                          cmap=cluster_cmap(k), vmin=0.5, vmax=k + 0.5, shading="auto")
            ax.coastlines(linewidth=0.4, color="black")
            ax.set_global()
            letter = "abcd"[2 * r + c]
            ax.set_title(f"({letter}) {method}, {row_name}: $k = {k}$", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")


if __name__ == "__main__":
    main()

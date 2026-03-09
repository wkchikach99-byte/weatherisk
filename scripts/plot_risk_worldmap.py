#!/usr/bin/env python3
"""Plot risk metrics (ES, exposure, clusters) on world maps with coastlines.

Uses the weatherisk risk_pipeline module to load and process the data,
then plots with Cartopy projections — matching the style of the
Risk_Analysis paper (Figures 1–2: global Expected Shortfall maps).

Outputs are saved to docs/figures/.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, BoundaryNorm, ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ── paths ────────────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

CSV_PATH = os.path.join(ROOT, "data", "risk_map_grid.csv")
OUT_DIR  = os.path.join(ROOT, "docs", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ── settings ─────────────────────────────────────────────────────────
N_BANDS      = 6
SMOOTH_SIGMA = 0.8
MIN_PATCH    = 30
DPI          = 300

# ── load data ────────────────────────────────────────────────────────
from weatherisk.risk_pipeline import (
    load_and_grid,
    smooth_field,
    quantile_bands,
    connected_patches,
    merge_tiny_regions,
    remap_ids_to_sequential,
    compute_cluster_stats,
)


def load_data():
    d = load_and_grid(CSV_PATH)
    df = d["df"]
    lats, lons = d["lats"], d["lons"]
    ES, VaR, EXP = d["ES"], d["VaR"], d["EXP"]
    land_mask = d["land_mask"]
    lon_grid, lat_grid = d["lon_grid"], d["lat_grid"]
    return df, lats, lons, ES, VaR, EXP, land_mask, lon_grid, lat_grid


def run_clustering(ES, land_mask, lon_grid, lat_grid):
    ES_s = smooth_field(ES, SMOOTH_SIGMA, land_mask)
    bands, _ = quantile_bands(ES_s, N_BANDS)
    cluster_id = connected_patches(bands, MIN_PATCH)
    cluster_id = merge_tiny_regions(cluster_id, lon_grid, lat_grid)
    cluster_id, K = remap_ids_to_sequential(cluster_id)
    return ES_s, cluster_id, K


# ── helper: world map axis ─────────────────────────────────────────
def make_map_ax(fig, pos=111, projection=None):
    """Create a GeoAxes with coastlines and gridlines."""
    if projection is None:
        projection = ccrs.Robinson()
    ax = fig.add_subplot(pos, projection=projection)
    ax.set_global()
    ax.coastlines(linewidth=0.4, color="0.3")
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor="0.5")
    gl = ax.gridlines(draw_labels=False, linewidth=0.15,
                      color="grey", alpha=0.5)
    return ax


# ════════════════════════════════════════════════════════════════════
#  PLOT 1 — Expected Shortfall (raw) on a world map
# ════════════════════════════════════════════════════════════════════
def plot_es_worldmap(lons, lats, ES, out_path):
    """Global map of Expected Shortfall (95 %) — log scale,
    matching Risk_Analysis paper Figure 1."""
    es_plot = ES.copy().astype(float)
    es_plot[es_plot <= 0] = np.nan

    fig = plt.figure(figsize=(14, 7))
    ax = make_map_ax(fig)

    vmin = np.nanpercentile(es_plot[np.isfinite(es_plot)], 1)
    vmax = np.nanpercentile(es_plot[np.isfinite(es_plot)], 99)
    if vmin <= 0:
        vmin = 1.0

    im = ax.pcolormesh(
        lons, lats, es_plot,
        transform=ccrs.PlateCarree(),
        cmap="YlOrRd",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        rasterized=True,
    )
    cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                      pad=0.04, shrink=0.7, aspect=35)
    cb.set_label("Expected Shortfall (ES$_{95}$)", fontsize=12)
    ax.set_title("Global Expected Shortfall (95 %)", fontsize=14, pad=12)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  PLOT 2 — Exposure (population) on a world map
# ════════════════════════════════════════════════════════════════════
def plot_exposure_worldmap(lons, lats, EXP, out_path):
    """Global map of exposure (population per cell)."""
    exp_plot = EXP.copy().astype(float)
    exp_plot[exp_plot <= 0] = np.nan

    fig = plt.figure(figsize=(14, 7))
    ax = make_map_ax(fig)

    vmin = np.nanpercentile(exp_plot[np.isfinite(exp_plot)], 5)
    vmax = np.nanpercentile(exp_plot[np.isfinite(exp_plot)], 99)
    if vmin <= 0:
        vmin = 1.0

    im = ax.pcolormesh(
        lons, lats, exp_plot,
        transform=ccrs.PlateCarree(),
        cmap="viridis",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        rasterized=True,
    )
    cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                      pad=0.04, shrink=0.7, aspect=35)
    cb.set_label("Population (exposure)", fontsize=12)
    ax.set_title("Global Exposure (Population per Grid Cell)", fontsize=14, pad=12)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  PLOT 3 — ES risk clusters on a world map
# ════════════════════════════════════════════════════════════════════
def plot_clusters_worldmap(lons, lats, cluster_id, K, out_path):
    """Global map of ES-based risk clusters (contiguous regions)."""
    cl_plot = cluster_id.astype(float)
    cl_plot[cluster_id < 0] = np.nan

    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    cmap = ListedColormap([base[i % len(base)] for i in range(K)])

    fig = plt.figure(figsize=(14, 7))
    ax = make_map_ax(fig)

    im = ax.pcolormesh(
        lons, lats, cl_plot,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=-0.5, vmax=K - 0.5,
        shading="auto",
        rasterized=True,
    )
    cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                      pad=0.04, shrink=0.7, aspect=35,
                      ticks=range(0, K, max(1, K // 10)))
    cb.set_label("Cluster ID", fontsize=12)
    ax.set_title(f"ES$_{{95}}$ Risk Clusters ({K} regions)", fontsize=14, pad=12)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  PLOT 4 — Choropleth: population-weighted ES per cluster
# ════════════════════════════════════════════════════════════════════
def plot_es_choropleth_worldmap(lons, lats, cluster_id, stats, out_path):
    """World-map choropleth coloured by pop-weighted mean ES per cluster."""
    metric_dict = dict(zip(
        stats["cluster"].astype(int),
        stats["ES_popw_mean"].astype(float),
    ))
    filled = np.full(cluster_id.shape, np.nan, dtype=float)
    for cid, val in metric_dict.items():
        filled[cluster_id == cid] = val

    v = filled[np.isfinite(filled)]

    fig = plt.figure(figsize=(14, 7))
    ax = make_map_ax(fig)

    if v.size > 0:
        vmin = max(np.nanpercentile(v, 2), 1e-6)
        vmax = np.nanpercentile(v, 98)
        norm = LogNorm(vmin=vmin, vmax=vmax) if vmax / vmin > 20 else None
        im = ax.pcolormesh(
            lons, lats, filled,
            transform=ccrs.PlateCarree(),
            cmap="inferno",
            norm=norm,
            shading="auto",
            rasterized=True,
        )
    else:
        im = ax.pcolormesh(
            lons, lats, filled,
            transform=ccrs.PlateCarree(),
            cmap="inferno",
            shading="auto",
            rasterized=True,
        )

    cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                      pad=0.04, shrink=0.7, aspect=35)
    cb.set_label("Population-weighted mean ES$_{95}$", fontsize=12)
    ax.set_title("Regional Expected Shortfall (Population-Weighted)",
                 fontsize=14, pad=12)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  PLOT 5 — Top-20 regions bar chart
# ════════════════════════════════════════════════════════════════════
def plot_top_regions_bar(stats, out_path):
    """Horizontal bar chart of the 20 highest-ES clusters."""
    top = stats.sort_values("ES_sum", ascending=False).head(20)

    fig, ax = plt.subplots(figsize=(9, 7))
    bars = ax.barh(
        top["cluster"].astype(str),
        top["ES_sum"] / 1e6,
        color=plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(top))),
        edgecolor="0.3", linewidth=0.4,
    )
    ax.invert_yaxis()
    ax.set_xlabel("Total ES$_{95}$ (millions)", fontsize=12)
    ax.set_ylabel("Cluster ID", fontsize=12)
    ax.set_title("Top-20 Risk Regions by Total Expected Shortfall",
                 fontsize=13, pad=10)
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  PLOT 6 — 2×2 summary panel (ES, Exposure, Clusters, Choropleth)
# ════════════════════════════════════════════════════════════════════
def plot_summary_panel(lons, lats, ES, EXP, cluster_id, K, stats, out_path):
    """Four-panel summary figure for the paper."""
    fig = plt.figure(figsize=(18, 12))
    proj = ccrs.Robinson()

    # --- Panel A: Raw ES ---
    ax1 = fig.add_subplot(2, 2, 1, projection=proj)
    ax1.set_global(); ax1.coastlines(lw=0.4, color="0.3")
    es_plot = ES.copy().astype(float); es_plot[es_plot <= 0] = np.nan
    v = es_plot[np.isfinite(es_plot)]
    vmin, vmax = max(np.nanpercentile(v, 1), 1), np.nanpercentile(v, 99)
    im1 = ax1.pcolormesh(lons, lats, es_plot, transform=ccrs.PlateCarree(),
                         cmap="YlOrRd", norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading="auto", rasterized=True)
    fig.colorbar(im1, ax=ax1, orientation="horizontal", pad=0.05, shrink=0.65)
    ax1.set_title("(a) Expected Shortfall (ES$_{95}$)", fontsize=12)

    # --- Panel B: Exposure ---
    ax2 = fig.add_subplot(2, 2, 2, projection=proj)
    ax2.set_global(); ax2.coastlines(lw=0.4, color="0.3")
    exp_plot = EXP.copy().astype(float); exp_plot[exp_plot <= 0] = np.nan
    v2 = exp_plot[np.isfinite(exp_plot)]
    vmin2, vmax2 = max(np.nanpercentile(v2, 5), 1), np.nanpercentile(v2, 99)
    im2 = ax2.pcolormesh(lons, lats, exp_plot, transform=ccrs.PlateCarree(),
                         cmap="viridis", norm=LogNorm(vmin=vmin2, vmax=vmax2),
                         shading="auto", rasterized=True)
    fig.colorbar(im2, ax=ax2, orientation="horizontal", pad=0.05, shrink=0.65)
    ax2.set_title("(b) Exposure (Population)", fontsize=12)

    # --- Panel C: Clusters ---
    ax3 = fig.add_subplot(2, 2, 3, projection=proj)
    ax3.set_global(); ax3.coastlines(lw=0.4, color="0.3")
    cl_plot = cluster_id.astype(float); cl_plot[cluster_id < 0] = np.nan
    cmap_cl = ListedColormap(
        [list(plt.get_cmap("tab20").colors)[i % 20] for i in range(K)]
    )
    im3 = ax3.pcolormesh(lons, lats, cl_plot, transform=ccrs.PlateCarree(),
                         cmap=cmap_cl, vmin=-0.5, vmax=K-0.5,
                         shading="auto", rasterized=True)
    fig.colorbar(im3, ax=ax3, orientation="horizontal", pad=0.05, shrink=0.65)
    ax3.set_title(f"(c) ES$_{{95}}$ Risk Clusters ({K} regions)", fontsize=12)

    # --- Panel D: Choropleth ---
    ax4 = fig.add_subplot(2, 2, 4, projection=proj)
    ax4.set_global(); ax4.coastlines(lw=0.4, color="0.3")
    metric_dict = dict(zip(stats["cluster"].astype(int),
                           stats["ES_popw_mean"].astype(float)))
    filled = np.full(cluster_id.shape, np.nan, dtype=float)
    for cid, val in metric_dict.items():
        filled[cluster_id == cid] = val
    vf = filled[np.isfinite(filled)]
    vmin4 = max(np.nanpercentile(vf, 2), 1e-6) if vf.size else 1
    vmax4 = np.nanpercentile(vf, 98) if vf.size else 10
    norm4 = LogNorm(vmin=vmin4, vmax=vmax4) if vmax4/vmin4 > 20 else None
    im4 = ax4.pcolormesh(lons, lats, filled, transform=ccrs.PlateCarree(),
                         cmap="inferno", norm=norm4,
                         shading="auto", rasterized=True)
    fig.colorbar(im4, ax=ax4, orientation="horizontal", pad=0.05, shrink=0.65)
    ax4.set_title("(d) Pop-Weighted Mean ES$_{95}$ per Region", fontsize=12)

    fig.suptitle("Climate Risk Analysis — Expected Shortfall Metrics",
                 fontsize=15, y=0.98)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {out_path}")


# ════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════
def main():
    print("Loading data …")
    df, lats, lons, ES, VaR, EXP, land_mask, lon_grid, lat_grid = load_data()
    print(f"  Grid: {len(lats)}×{len(lons)} = {len(lats)*len(lons):,} cells")
    print(f"  Land cells (exposure > 0): {(EXP > 0).sum():,}")

    print("Running ES clustering …")
    ES_s, cluster_id, K = run_clustering(ES, land_mask, lon_grid, lat_grid)
    print(f"  → {K} contiguous risk regions")

    print("Computing per-cluster statistics …")
    out_df = pd.DataFrame({
        "lat": lat_grid.ravel(),
        "lon": lon_grid.ravel(),
        "cluster": cluster_id.ravel(),
    })
    stats = compute_cluster_stats(out_df, df)
    print(f"  Top-5 clusters by total ES:")
    top5 = stats.sort_values("ES_sum", ascending=False).head(5)
    for _, r in top5.iterrows():
        print(f"    Cluster {int(r['cluster']):3d}: "
              f"ES_sum={r['ES_sum']/1e6:10.2f}M  "
              f"cells={int(r['n_cells']):6d}  "
              f"pop={r['exposure_sum']/1e6:8.2f}M")

    print("\nGenerating world-map plots …")
    plot_es_worldmap(
        lons, lats, ES,
        os.path.join(OUT_DIR, "risk_es_worldmap.pdf"))
    plot_exposure_worldmap(
        lons, lats, EXP,
        os.path.join(OUT_DIR, "risk_exposure_worldmap.pdf"))
    plot_clusters_worldmap(
        lons, lats, cluster_id, K,
        os.path.join(OUT_DIR, "risk_clusters_worldmap.pdf"))
    plot_es_choropleth_worldmap(
        lons, lats, cluster_id, stats,
        os.path.join(OUT_DIR, "risk_es_choropleth_worldmap.pdf"))
    plot_top_regions_bar(
        stats,
        os.path.join(OUT_DIR, "risk_top20_bar.pdf"))
    plot_summary_panel(
        lons, lats, ES, EXP, cluster_id, K, stats,
        os.path.join(OUT_DIR, "risk_summary_panel.pdf"))

    print("\n✔ All plots saved to docs/figures/")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate clear, well-annotated LEC cluster risk maps.

Focuses on making the per-cluster risk metrics from the CPC pipeline
visible and interpretable for a first-time reader.

Produces 4 annotated maps + 1 summary panel:
  1. LEC clusters (k=26) with cluster-ID labels
  2. Per-cluster ES₉₅ of spatial block maximum (mm/day)
  3. Per-cell ES₉₅ — marginal hazard (mm/day)
  4. Per-cluster spatial coherence (co-exceedance rate)
  5. Summary panel combining all four

Usage:
    python scripts/cluster_risk_maps.py
"""

import os, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import run_cpc_pipeline, PipelineConfig
from weatherisk.risk import compute_var, compute_es

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

OUT_DIR = "docs/figures/cluster_risk"
os.makedirs(OUT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
#  1. Run pipeline
# ═══════════════════════════════════════════════════════════════

print("=" * 65)
print("  LEC CLUSTER RISK MAPS")
print("=" * 65)

t0 = time.time()

cfg = PipelineConfig(
    gdp_path="data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc",
    gdp_year=2015,
)
print("\n[1/5] Running CPC pipeline (this takes ~4 min) ...")
result = run_cpc_pipeline(cfg, verbose=False)

bm = result["bm"]              # (20, n_cells) mm/day annual block maxima
frechet = result["frechet"]     # (20, n_cells) Fréchet-transformed
labels = result["labels_lec"]   # (n_cells,) cluster labels 1..k
k = result["k_lec"]
lat = result["lat"]             # (n_cells,) geographic latitudes
lon = result["lon"]             # (n_cells,) geographic longitudes
gdp = result["gdp_per_cell"]
land_idx = result["land_idx"]
lats_1d = result["lats_1d"]
lons_1d = result["lons_1d"]
n_lat = result["n_lat"]
n_lon = result["n_lon"]

n_years, n_cells = bm.shape
clusters = sorted(np.unique(labels))
n_clusters = len(clusters)

print(f"    {n_years} years (2000–2019), {n_cells} land cells, "
      f"k_LEC = {n_clusters} clusters")
print(f"    Block maxima: [{bm.min():.1f}, {bm.max():.1f}] mm/day")
print(f"    Fréchet:      [{frechet.min():.2f}, {frechet.max():.2f}]")


# ═══════════════════════════════════════════════════════════════
#  2. Compute risk metrics per cluster AND per cell
# ═══════════════════════════════════════════════════════════════

print("\n[2/5] Computing risk metrics ...")

# --- Per-cell ES₉₅ in mm/day ---
cell_es_mm = np.array([compute_es(bm[:, c], 0.95) for c in range(n_cells)])
cell_q90 = np.quantile(bm, 0.90, axis=0)  # 90th pct per cell

# --- Per-cluster metrics ---
cluster_info = []
for cl in clusters:
    mask = labels == cl
    idx = np.where(mask)[0]
    n_cl = len(idx)

    bm_cl = bm[:, mask]          # (20, n_cl) mm/day
    fr_cl = frechet[:, mask]     # (20, n_cl) Fréchet

    # Spatial block maximum per year (the "loss" for this cluster)
    L_mm = bm_cl.max(axis=1)     # (20,) in mm/day
    L_fr = fr_cl.max(axis=1)     # (20,) in Fréchet units

    # ES₉₅ of spatial maximum
    es_mm = compute_es(L_mm, 0.95)
    var_mm = compute_var(L_mm, 0.95)
    es_fr = compute_es(L_fr, 0.95)

    # Mean per-cell ES within cluster
    mean_cell_es = cell_es_mm[mask].mean()

    # Co-exceedance: in years with any cell > q90, fraction exceeding
    exceed = bm_cl >= cell_q90[mask][None, :]
    any_exceed = exceed.any(axis=1)
    if any_exceed.sum() > 0:
        co_exceedance = exceed[any_exceed].mean(axis=1).mean()
    else:
        co_exceedance = 0.0

    # Spearman ρ̄
    if n_cl >= 2:
        from scipy.stats import spearmanr
        rho_mat, _ = spearmanr(bm_cl)
        if n_cl == 2:
            mean_rho = float(rho_mat)
        else:
            tri = np.triu_indices(n_cl, k=1)
            mean_rho = rho_mat[tri].mean()
    else:
        mean_rho = 1.0

    # Worst year (year of spatial max)
    worst_year_idx = np.argmax(L_mm)
    worst_year = 2000 + worst_year_idx
    worst_year_mm = L_mm[worst_year_idx]

    # Cluster centroid (for labelling)
    centroid_lat = lat[mask].mean()
    centroid_lon = lon[mask].mean()

    # GDP
    gdp_total = gdp[mask].sum() if gdp is not None else 0.0

    cluster_info.append(dict(
        cluster=cl, n_cells=n_cl,
        es_mm=es_mm, var_mm=var_mm, es_fr=es_fr,
        mean_cell_es=mean_cell_es,
        co_exceedance=co_exceedance,
        mean_rho=mean_rho,
        worst_year=worst_year, worst_mm=worst_year_mm,
        centroid_lat=centroid_lat, centroid_lon=centroid_lon,
        gdp_total=gdp_total,
    ))

# Print summary table
print(f"\n    {'Cl':>3s} {'N':>4s} {'ES₉₅mm':>8s} {'VaR₉₅':>7s} "
      f"{'CellES':>7s} {'ρ̄':>6s} {'CoExc':>6s} {'Worst':>5s} {'Peak':>7s}")
print("    " + "-" * 65)
for s in sorted(cluster_info, key=lambda x: x["es_mm"], reverse=True):
    print(f"    {s['cluster']:3d} {s['n_cells']:4d} "
          f"{s['es_mm']:8.1f} {s['var_mm']:7.1f} "
          f"{s['mean_cell_es']:7.1f} {s['mean_rho']:6.3f} "
          f"{s['co_exceedance']:6.3f} {s['worst_year']:5d} "
          f"{s['worst_mm']:7.1f}")


# ═══════════════════════════════════════════════════════════════
#  3. Helper functions for plotting
# ═══════════════════════════════════════════════════════════════

print("\n[3/5] Preparing maps ...")

extent = [cfg.lon_range[0] - 2, cfg.lon_range[1] + 2,
          cfg.lat_range[0] - 2, cfg.lat_range[1] + 2]


def to_grid(values):
    """Map cell-level values to the full (n_lat, n_lon) grid."""
    grid = np.full(n_lat * n_lon, np.nan)
    grid[land_idx] = values
    return grid.reshape(n_lat, n_lon)


def make_edges(centres):
    c = np.asarray(centres, dtype=float)
    d = np.diff(c)
    e = np.empty(len(c) + 1)
    e[1:-1] = c[:-1] + d / 2
    e[0] = c[0] - d[0] / 2
    e[-1] = c[-1] + d[-1] / 2
    return e


lat_e = make_edges(lats_1d)
lon_e = make_edges(lons_1d)


def setup_ax(fig, pos, title):
    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
    ax.add_feature(cfeature.LAND, facecolor="#f5f0e8", alpha=0.15)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.7, color="#444444")
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--",
                   color="#888888")
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, alpha=0.4,
                      color="#999999")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {"size": 8}
    gl.ylabel_style = {"size": 8}
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    return ax


def add_cluster_labels(ax, cluster_info, field="cluster", fmt="{:.0f}",
                       fontsize=7, color="black"):
    """Add centroid labels on each cluster."""
    text_effect = [pe.withStroke(linewidth=2.5, foreground="white")]
    for s in cluster_info:
        val = s[field]
        txt = fmt.format(val) if isinstance(val, float) else str(val)
        ax.text(s["centroid_lon"], s["centroid_lat"], txt,
                fontsize=fontsize, ha="center", va="center",
                fontweight="bold", color=color,
                transform=ccrs.PlateCarree(),
                path_effects=text_effect)


# ═══════════════════════════════════════════════════════════════
#  4. Individual maps
# ═══════════════════════════════════════════════════════════════

print("\n[4/5] Generating individual maps ...")

# ── MAP 1: LEC clusters with IDs ────────────────────────────

fig = plt.figure(figsize=(12, 7.5))
ax = setup_ax(fig, 111,
    f"Map 1 — LEC Dependence Clusters (k = {n_clusters})\n"
    "Cells grouped by similarity of estimated anisotropy ellipses\n"
    "Contzen et al. (2025) method · CPC precipitation 2000–2019")

base_colors = (list(plt.get_cmap("tab20").colors)
              + list(plt.get_cmap("Set3").colors))
cmap_cl = ListedColormap([base_colors[i % len(base_colors)]
                          for i in range(k)])
norm_cl = BoundaryNorm(np.arange(0.5, k + 1.5), cmap_cl.N)

grid = to_grid(labels.astype(float))
masked = np.ma.masked_invalid(grid)
ax.pcolormesh(lon_e, lat_e, masked, cmap=cmap_cl, norm=norm_cl,
              transform=ccrs.PlateCarree(), shading="flat")

# Label each cluster with its ID and cell count
text_effect = [pe.withStroke(linewidth=2.5, foreground="white")]
for s in cluster_info:
    ax.text(s["centroid_lon"], s["centroid_lat"],
            f"{s['cluster']}\n({s['n_cells']})",
            fontsize=6.5, ha="center", va="center",
            fontweight="bold", color="black",
            transform=ccrs.PlateCarree(), path_effects=text_effect)

# Add a note
ax.text(0.01, -0.06, 
        "Numbers: cluster ID (cell count). "
        "Spatial contiguity is not enforced — clusters are defined by\n"
        "ellipse-shape similarity. Domain: 30–65°N, 5–55°E, ~2° resolution.",
        fontsize=8, transform=ax.transAxes, color="#555555",
        style="italic")

p = os.path.join(OUT_DIR, "map1_lec_clusters.png")
fig.savefig(p, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── MAP 2: Per-cluster ES₉₅ of spatial maximum (mm/day) ─────

fig = plt.figure(figsize=(12, 7.5))
ax = setup_ax(fig, 111,
    "Map 2 — Cluster-Level Tail Risk: Expected Shortfall\n"
    "ES₉₅ of the annual spatial block maximum $L_t = \\max_{s \\in A} X_t(s)$\n"
    "Computed per LEC cluster from 20 annual maxima (2000–2019), mm/day")

# Build per-cell array from cluster-level ES
es_cluster_per_cell = np.zeros(n_cells)
for s in cluster_info:
    es_cluster_per_cell[labels == s["cluster"]] = s["es_mm"]

grid = to_grid(es_cluster_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7, pad=0.02)
cb.set_label("Expected Shortfall ES₉₅ (mm/day)", fontsize=10)

# Label each cluster with its ES value and worst year
text_effect = [pe.withStroke(linewidth=2.5, foreground="white")]
for s in cluster_info:
    ax.text(s["centroid_lon"], s["centroid_lat"],
            f"{s['es_mm']:.0f}\n({s['worst_year']})",
            fontsize=6, ha="center", va="center",
            fontweight="bold", color="black",
            transform=ccrs.PlateCarree(), path_effects=text_effect)

ax.text(0.01, -0.06,
        "Numbers: ES₉₅ in mm/day (worst year in the 2000–2019 record).\n"
        "ES₉₅ = mean of annual spatial maxima exceeding the 95th percentile. "
        "Higher values → more intense cluster-wide extremes.",
        fontsize=8, transform=ax.transAxes, color="#555555", style="italic")

p = os.path.join(OUT_DIR, "map2_cluster_es.png")
fig.savefig(p, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── MAP 3: Per-cell ES₉₅ (marginal hazard, mm/day) ──────────

fig = plt.figure(figsize=(12, 7.5))
ax = setup_ax(fig, 111,
    "Map 3 — Per-Cell Marginal Hazard: Expected Shortfall\n"
    "ES₉₅ of each cell's own 20 annual block maxima (2000–2019)\n"
    "Independent of cluster assignment — purely local tail intensity")

grid = to_grid(cell_es_mm)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7, pad=0.02)
cb.set_label("Per-cell ES₉₅ (mm/day)", fontsize=10)

ax.text(0.01, -0.06,
        "Each cell independently: ES₉₅ = mean of annual maxima above "
        "the 95th percentile (top 1–2 years out of 20).\n"
        "This is the marginal (single-site) tail risk — it does not "
        "capture spatial co-occurrence of extremes.",
        fontsize=8, transform=ax.transAxes, color="#555555", style="italic")

p = os.path.join(OUT_DIR, "map3_cell_es.png")
fig.savefig(p, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── MAP 4: Co-exceedance rate per cluster ────────────────────

fig = plt.figure(figsize=(12, 7.5))
ax = setup_ax(fig, 111,
    "Map 4 — Spatial Coherence: Co-Exceedance Rate\n"
    "During years with at least one tail event (any cell > its 90th pct),\n"
    "what fraction of cells in the cluster simultaneously exceed?")

coexc_per_cell = np.zeros(n_cells)
for s in cluster_info:
    coexc_per_cell[labels == s["cluster"]] = s["co_exceedance"]

grid = to_grid(coexc_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="OrRd",
                     vmin=0.05, vmax=0.70,
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7, pad=0.02)
cb.set_label("Co-exceedance rate (fraction of cluster)", fontsize=10)

# Label with rate
text_effect = [pe.withStroke(linewidth=2.5, foreground="white")]
for s in cluster_info:
    ax.text(s["centroid_lon"], s["centroid_lat"],
            f"{s['co_exceedance']:.0%}",
            fontsize=6.5, ha="center", va="center",
            fontweight="bold", color="black",
            transform=ccrs.PlateCarree(), path_effects=text_effect)

ax.text(0.01, -0.06,
        "Higher values mean extremes are spatially widespread within "
        "the cluster (coherent hazard).\n"
        "Lower values mean extremes tend to be localised to a few "
        "cells (incoherent / diversified).",
        fontsize=8, transform=ax.transAxes, color="#555555", style="italic")

p = os.path.join(OUT_DIR, "map4_coexceedance.png")
fig.savefig(p, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── MAP 5: Mean cell ES₉₅ per cluster (size-bias-free) ──────

fig = plt.figure(figsize=(12, 7.5))
ax = setup_ax(fig, 111,
    "Map 5 — Cluster Hazard Intensity (Size-Bias-Free)\n"
    "Mean per-cell ES₉₅ within each cluster — removes mechanical\n"
    "inflation from taking spatial max over many cells")

mean_es_per_cell = np.zeros(n_cells)
for s in cluster_info:
    mean_es_per_cell[labels == s["cluster"]] = s["mean_cell_es"]

grid = to_grid(mean_es_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7, pad=0.02)
cb.set_label("Mean per-cell ES₉₅ (mm/day)", fontsize=10)

text_effect = [pe.withStroke(linewidth=2.5, foreground="white")]
for s in cluster_info:
    ax.text(s["centroid_lon"], s["centroid_lat"],
            f"{s['mean_cell_es']:.0f}",
            fontsize=6.5, ha="center", va="center",
            fontweight="bold", color="black",
            transform=ccrs.PlateCarree(), path_effects=text_effect)

ax.text(0.01, -0.06,
        "For each cluster: average the per-cell ES₉₅ over all cells "
        "in the cluster.\nThis measures the typical local tail severity "
        "without inflating large clusters.",
        fontsize=8, transform=ax.transAxes, color="#555555", style="italic")

p = os.path.join(OUT_DIR, "map5_mean_cell_es.png")
fig.savefig(p, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ═══════════════════════════════════════════════════════════════
#  6. Combined 3×2 summary panel
# ═══════════════════════════════════════════════════════════════

print("\n[5/6] Generating summary panel ...")

fig = plt.figure(figsize=(22, 20))
fig.suptitle(
    "LEC Cluster Risk Analysis — CPC Precipitation Extremes (2000–2019)\n"
    f"Europe / Middle East · {n_cells} land cells · {n_clusters} LEC clusters "
    f"· 20 annual block maxima",
    fontsize=15, fontweight="bold", y=0.98)

# --- Panel (a): LEC clusters ---
ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.PlateCarree())
ax1.set_extent(extent, crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax1.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(labels.astype(float))
masked = np.ma.masked_invalid(grid)
ax1.pcolormesh(lon_e, lat_e, masked, cmap=cmap_cl, norm=norm_cl,
               transform=ccrs.PlateCarree(), shading="flat")
te = [pe.withStroke(linewidth=2, foreground="white")]
for s in cluster_info:
    ax1.text(s["centroid_lon"], s["centroid_lat"],
             f"{s['cluster']}", fontsize=5.5, ha="center", va="center",
             fontweight="bold", transform=ccrs.PlateCarree(),
             path_effects=te)
ax1.set_title(f"(a) LEC Clusters (k = {n_clusters})", fontsize=12,
              fontweight="bold")

# --- Panel (b): Cluster-level ES₉₅ (spatial max) ---
ax2 = fig.add_subplot(3, 2, 2, projection=ccrs.PlateCarree())
ax2.set_extent(extent, crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax2.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(es_cluster_per_cell)
masked = np.ma.masked_invalid(grid)
m2 = ax2.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb2 = plt.colorbar(m2, ax=ax2, shrink=0.6, pad=0.02)
cb2.set_label("ES₉₅ (mm/day)")
for s in cluster_info:
    ax2.text(s["centroid_lon"], s["centroid_lat"],
             f"{s['es_mm']:.0f}", fontsize=5.5, ha="center", va="center",
             fontweight="bold", transform=ccrs.PlateCarree(),
             path_effects=te)
ax2.set_title("(b) Cluster ES₉₅ — Spatial Max (mm/day)\n"
              "[Biased by cluster size]",
              fontsize=11, fontweight="bold")

# --- Panel (c): Mean cell ES per cluster (size-bias-free) ---
ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
ax3.set_extent(extent, crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax3.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(mean_es_per_cell)
masked = np.ma.masked_invalid(grid)
m3 = ax3.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb3 = plt.colorbar(m3, ax=ax3, shrink=0.6, pad=0.02)
cb3.set_label("Mean cell ES₉₅ (mm/day)")
for s in cluster_info:
    ax3.text(s["centroid_lon"], s["centroid_lat"],
             f"{s['mean_cell_es']:.0f}", fontsize=5.5,
             ha="center", va="center", fontweight="bold",
             transform=ccrs.PlateCarree(), path_effects=te)
ax3.set_title("(c) Mean Cell ES₉₅ per Cluster (mm/day)\n"
              "[Size-bias-free hazard intensity]",
              fontsize=11, fontweight="bold")

# --- Panel (d): Per-cell ES₉₅ ---
ax4 = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree())
ax4.set_extent(extent, crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax4.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(cell_es_mm)
masked = np.ma.masked_invalid(grid)
m4 = ax4.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
cb4 = plt.colorbar(m4, ax=ax4, shrink=0.6, pad=0.02)
cb4.set_label("ES₉₅ (mm/day)")
ax4.set_title("(d) Per-Cell ES₉₅ — Marginal Hazard (mm/day)",
              fontsize=11, fontweight="bold")

# --- Panel (e): Co-exceedance ---
ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree())
ax5.set_extent(extent, crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax5.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax5.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(coexc_per_cell)
masked = np.ma.masked_invalid(grid)
m5 = ax5.pcolormesh(lon_e, lat_e, masked, cmap="OrRd",
                     vmin=0.05, vmax=0.70,
                     transform=ccrs.PlateCarree(), shading="flat")
cb5 = plt.colorbar(m5, ax=ax5, shrink=0.6, pad=0.02)
cb5.set_label("Co-exceedance rate")
for s in cluster_info:
    ax5.text(s["centroid_lon"], s["centroid_lat"],
             f"{s['co_exceedance']:.0%}", fontsize=5.5,
             ha="center", va="center", fontweight="bold",
             transform=ccrs.PlateCarree(), path_effects=te)
ax5.set_title("(e) Co-Exceedance Rate per Cluster",
              fontsize=11, fontweight="bold")

# --- Panel (f): Spearman ρ̄ ---
rho_per_cell = np.zeros(n_cells)
for s in cluster_info:
    rho_per_cell[labels == s["cluster"]] = s["mean_rho"]

ax6 = fig.add_subplot(3, 2, 6, projection=ccrs.PlateCarree())
ax6.set_extent(extent, crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.OCEAN, facecolor="#ddeeff", alpha=0.4)
ax6.add_feature(cfeature.COASTLINE, linewidth=0.6)
ax6.add_feature(cfeature.BORDERS, linewidth=0.25, linestyle="--")
grid = to_grid(rho_per_cell)
masked = np.ma.masked_invalid(grid)
m6 = ax6.pcolormesh(lon_e, lat_e, masked, cmap="RdYlBu_r",
                     vmin=-0.5, vmax=0.5,
                     transform=ccrs.PlateCarree(), shading="flat")
cb6 = plt.colorbar(m6, ax=ax6, shrink=0.6, pad=0.02)
cb6.set_label("Spearman ρ̄")
for s in cluster_info:
    if s["n_cells"] >= 2:
        ax6.text(s["centroid_lon"], s["centroid_lat"],
                 f"{s['mean_rho']:.2f}", fontsize=5.5,
                 ha="center", va="center", fontweight="bold",
                 transform=ccrs.PlateCarree(), path_effects=te)
ax6.set_title("(f) Spatial Coherence — Mean Spearman ρ̄",
              fontsize=11, fontweight="bold")

plt.tight_layout(rect=[0, 0.04, 1, 0.95])

fig.text(0.5, 0.008,
    "Data: NOAA CPC daily precipitation (2000–2019), 30–65°N, 5–55°E, "
    "coarsened to ~2°, 384 land cells.\n"
    "Method: annual block maxima → GEV → Fréchet → local MLE(a,b,γ) → "
    "LEC clustering (30% quantile threshold) → risk metrics.\n"
    "ES₉₅ = mean of annual maxima exceeding the 95th percentile. "
    "Co-exceedance = fraction of cluster cells simultaneously above "
    "their 90th percentile. ρ̄ = mean pairwise Spearman correlation.",
    ha="center", fontsize=8, style="italic", color="#555555")

p = os.path.join(OUT_DIR, "summary_6panel.png")
fig.savefig(p, dpi=250, bbox_inches="tight")
plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ═══════════════════════════════════════════════════════════════
#  7. Comprehensive analytical table
# ═══════════════════════════════════════════════════════════════

print("\n[6/6] Writing comprehensive analytical report ...")

sorted_info = sorted(cluster_info, key=lambda x: x["mean_cell_es"],
                     reverse=True)

def _geo_desc(lat, lon):
    """Rough geographic label from centroid coordinates."""
    labels_geo = []
    if lat > 55: labels_geo.append("Northern Europe")
    elif lat > 45: labels_geo.append("Central Europe")
    elif lat > 35: labels_geo.append("Southern/Mediterranean")
    else: labels_geo.append("North Africa/Middle East")
    if lon > 40: labels_geo.append("East")
    elif lon > 25: labels_geo.append("Central-East")
    elif lon > 15: labels_geo.append("Central")
    else: labels_geo.append("West")
    return ", ".join(labels_geo)

lines = []
lines.append("=" * 85)
lines.append("  LEC CLUSTER RISK ANALYSIS — COMPREHENSIVE REPORT")
lines.append("=" * 85)
lines.append("")
lines.append("DATA SOURCE")
lines.append("  Dataset     : NOAA CPC Global Unified Gauge-Based Analysis")
lines.append("  Variable    : Daily precipitation (mm/day)")
lines.append("  Period      : 2000–2019 (20 years)")
lines.append("  Domain      : 30°N–65°N, 5°E–55°E (Europe / Middle East)")
lines.append("  Resolution  : 0.5° coarsened ×4 → ~2° grid")
lines.append(f"  Land cells  : {n_cells}")
lines.append(f"  Block maxima: {n_years} annual maxima per cell "
             f"(block size = 365 days)")
lines.append(f"  BM range    : [{bm.min():.1f}, {bm.max():.1f}] mm/day")
lines.append("")
lines.append("CLUSTERING")
lines.append("  Method      : LEC (Localised Ellipse-shape Clustering)")
lines.append("  Reference   : Contzen, Dickhaus & Lohmeyer (2025, Extremes)")
lines.append("  Parameters  : ν=5, α=1.0, ε=3.0 (normalised units)")
lines.append("  Threshold   : 30th percentile of pairwise LEC dissimilarities")
lines.append(f"  Result      : k_LEC = {n_clusters} clusters")
lines.append("")
lines.append("-" * 85)
lines.append("")

# ── Main summary table ──
lines.append("TABLE 1: CLUSTER RISK SUMMARY")
lines.append("Sorted by mean per-cell ES₉₅ (descending) — "
             "size-bias-free hazard ranking")
lines.append("")
hdr = (f"{'Cl':>3s}  {'N':>4s}  {'MeanES':>7s}  {'MaxES':>6s}  "
       f"{'SpatES':>7s}  {'VaR₉₅':>7s}  {'ρ̄':>6s}  {'CoExc':>6s}  "
       f"{'Worst':>5s}  {'Peak':>7s}  {'Region'}")
lines.append(hdr)
lines.append("-" * 95)

for s in sorted_info:
    region = _geo_desc(s["centroid_lat"], s["centroid_lon"])
    lines.append(
        f"{s['cluster']:3d}  {s['n_cells']:4d}  "
        f"{s['mean_cell_es']:7.1f}  {cell_es_mm[labels == s['cluster']].max():6.1f}  "
        f"{s['es_mm']:7.1f}  {s['var_mm']:7.1f}  "
        f"{s['mean_rho']:6.3f}  {s['co_exceedance']:6.3f}  "
        f"{s['worst_year']:5d}  {s['worst_mm']:7.1f}  {region}"
    )

lines.append("")
lines.append("")

# ── Column definitions ──
lines.append("COLUMN DEFINITIONS")
lines.append("-" * 85)
lines.append("")
lines.append("  Cl      Cluster ID assigned by the LEC algorithm (arbitrary "
             "label, not ranked).")
lines.append("")
lines.append("  N       Number of ~2° grid cells in this cluster. The LEC "
             "algorithm groups cells")
lines.append("          by anisotropy ellipse similarity. Spatial contiguity "
             "is NOT enforced.")
lines.append("")
lines.append("  MeanES  Mean per-cell ES₉₅ (mm/day). For each cell in the "
             "cluster, take its 20")
lines.append("          annual block maxima, compute ES₉₅ = mean of values "
             "≥ 95th percentile,")
lines.append("          then average across all cells in the cluster.")
lines.append("          → This is the SIZE-BIAS-FREE hazard intensity. "
             "Best metric for comparing")
lines.append("          clusters of different sizes.")
lines.append("")
lines.append("  MaxES   Maximum per-cell ES₉₅ (mm/day) within the cluster. "
             "The single most")
lines.append("          extreme location.")
lines.append("")
lines.append("  SpatES  ES₉₅ of the spatial block maximum (mm/day). "
             "For each year t, define")
lines.append("          L_t = max over all cells s in cluster of X_t(s). "
             "Then ES₉₅ = mean of")
lines.append("          L_t values exceeding the 95th percentile of "
             "{L_1, ..., L_20}.")
lines.append("          ⚠ CAUTION: Mechanically inflated by cluster size. "
             "A cluster with 155 cells")
lines.append("          will always have a larger spatial max than one with "
             "5 cells, even if the")
lines.append("          individual cells are equally extreme. "
             "Use MeanES for fair comparison.")
lines.append("")
lines.append("  VaR₉₅  Value at Risk at 95% of the spatial block maximum. "
             "The 95th percentile")
lines.append("          of {L_1, ..., L_20}. With 20 years, this is "
             "approximately the 19th")
lines.append("          largest annual spatial max.")
lines.append("")
lines.append("  ρ̄       Mean pairwise Spearman rank correlation of "
             "annual block maxima time")
lines.append("          series within the cluster. Measures whether "
             "cells tend to be extreme in")
lines.append("          the same years (ρ̄ > 0 → simultaneous, "
             "ρ̄ ≈ 0 → independent,")
lines.append("          ρ̄ < 0 → anti-correlated).")
lines.append("          → Single-cell clusters trivially have ρ̄ = 1.000 "
             "(only self-correlation).")
lines.append("")
lines.append("  CoExc   Co-exceedance rate. In years where at least one "
             "cell exceeds its own")
lines.append("          90th percentile, what fraction of all cluster "
             "cells simultaneously exceed")
lines.append("          their own 90th percentile?")
lines.append("          Example: CoExc = 0.100 means 10% of cells "
             "co-exceedance on average.")
lines.append("          → Single-cell clusters trivially have "
             "CoExc = 1.000.")
lines.append("")
lines.append("  Worst   The calendar year (2000–2019) in which the "
             "cluster's spatial block")
lines.append("          maximum was highest — the 'worst event year' "
             "for this cluster.")
lines.append("")
lines.append("  Peak    Spatial block maximum in the worst year (mm/day). "
             "The single highest")
lines.append("          block maximum observed across all cells in "
             "the cluster in any year.")
lines.append("")
lines.append("  Region  Approximate geographic label based on "
             "cluster centroid coordinates.")
lines.append("")
lines.append("")

# ── Per-cluster detailed cards ──
lines.append("=" * 85)
lines.append("  DETAILED PER-CLUSTER ANALYSIS CARDS")
lines.append("=" * 85)
lines.append("")

for i, s in enumerate(sorted_info):
    cl = s["cluster"]
    mask = labels == cl
    idx = np.where(mask)[0]
    bm_cl = bm[:, mask]
    L_mm = bm_cl.max(axis=1)  # spatial max per year

    region = _geo_desc(s["centroid_lat"], s["centroid_lon"])

    lines.append(f"  CLUSTER {cl}  |  {s['n_cells']} cells  |  "
                 f"{region}")
    lines.append(f"  Centroid: ({s['centroid_lat']:.1f}°N, "
                 f"{s['centroid_lon']:.1f}°E)")
    lines.append(f"  " + "-" * 60)
    lines.append(f"")
    lines.append(f"  Hazard metrics:")
    lines.append(f"    Mean per-cell ES₉₅ = {s['mean_cell_es']:.1f} mm/day")
    lines.append(f"    Max per-cell ES₉₅  = "
                 f"{cell_es_mm[mask].max():.1f} mm/day")
    lines.append(f"    Min per-cell ES₉₅  = "
                 f"{cell_es_mm[mask].min():.1f} mm/day")
    lines.append(f"    Std per-cell ES₉₅  = "
                 f"{cell_es_mm[mask].std():.1f} mm/day")
    lines.append(f"    Spatial max ES₉₅   = {s['es_mm']:.1f} mm/day  "
                 f"(⚠ size-biased)")
    lines.append(f"    Spatial max VaR₉₅  = {s['var_mm']:.1f} mm/day")
    lines.append(f"")
    lines.append(f"  Spatial coherence:")
    if s["n_cells"] >= 2:
        lines.append(f"    Mean Spearman ρ̄  = {s['mean_rho']:.3f}")
        lines.append(f"    Co-exceedance    = "
                     f"{s['co_exceedance']:.1%}")
        if s["mean_rho"] > 0.1:
            lines.append(f"    → Moderately coherent: extremes tend "
                         f"to co-occur across the cluster")
        elif s["mean_rho"] > -0.05:
            lines.append(f"    → Weakly coherent / near-independent: "
                         f"extremes are mostly localised")
        else:
            lines.append(f"    → Anti-correlated: when one part is "
                         f"wet, another tends to be dry")
    else:
        lines.append(f"    Single-cell cluster — coherence metrics "
                     f"not applicable")
    lines.append(f"")

    # Year-by-year spatial max
    lines.append(f"  Year-by-year spatial block maximum L_t (mm/day):")
    for t in range(n_years):
        year = 2000 + t
        val = L_mm[t]
        marker = " ◀ VaR" if val >= s["var_mm"] else ""
        marker = " ◀◀ WORST" if t == np.argmax(L_mm) else marker
        lines.append(f"    {year}: {val:7.1f}{marker}")

    if s["gdp_total"] > 0:
        def _fmt_gdp(v):
            if v >= 1e12: return f"${v/1e12:.1f}T"
            if v >= 1e9: return f"${v/1e9:.1f}B"
            if v >= 1e6: return f"${v/1e6:.0f}M"
            return f"${v:,.0f}"
        lines.append(f"")
        lines.append(f"  Economic exposure:")
        lines.append(f"    Total GDP (PPP 2015) = "
                     f"{_fmt_gdp(s['gdp_total'])}")
        lines.append(f"    Mean GDP per cell    = "
                     f"{_fmt_gdp(s['gdp_total'] / s['n_cells'])}")

    lines.append(f"")
    lines.append(f"")

# ── Sanity checks ──
lines.append("=" * 85)
lines.append("  SANITY CHECKS & DIAGNOSTICS")
lines.append("=" * 85)
lines.append("")

# Size bias diagnostic
multi = [s for s in cluster_info if s["n_cells"] >= 2]
sizes = np.array([s["n_cells"] for s in multi])
spat_es = np.array([s["es_mm"] for s in multi])
mean_es = np.array([s["mean_cell_es"] for s in multi])
from scipy.stats import spearmanr as _sp
rho_size_spat, p_size_spat = _sp(sizes, spat_es)
rho_size_mean, p_size_mean = _sp(sizes, mean_es)

lines.append("  1. SIZE BIAS CHECK")
lines.append(f"     Spearman correlation (cluster size vs spatial-max ES₉₅): "
             f"ρ = {rho_size_spat:.3f}, p = {p_size_spat:.4f}")
lines.append(f"     Spearman correlation (cluster size vs mean-cell ES₉₅):   "
             f"ρ = {rho_size_mean:.3f}, p = {p_size_mean:.4f}")
if abs(rho_size_spat) > 0.5:
    lines.append(f"     → Spatial-max ES is {'' if rho_size_spat > 0 else 'negatively '}"
                 f"correlated with cluster size — SIZE BIAS CONFIRMED")
if abs(rho_size_mean) < 0.3:
    lines.append(f"     → Mean-cell ES is weakly correlated with size "
                 f"— size bias successfully removed")
lines.append("")

# Worst year distribution
lines.append("  2. WORST YEAR DISTRIBUTION")
lines.append("     Which calendar years produce the most 'worst events' "
             "across clusters?")
from collections import Counter
worst_years = Counter(s["worst_year"] for s in cluster_info)
for year, count in worst_years.most_common():
    lines.append(f"     {year}: {count} cluster(s)")
lines.append("")

# ES plausibility
lines.append("  3. ES PLAUSIBILITY CHECK")
lines.append(f"     Overall block maxima range: "
             f"[{bm.min():.1f}, {bm.max():.1f}] mm/day")
lines.append(f"     Per-cell ES₉₅ range: "
             f"[{cell_es_mm.min():.1f}, {cell_es_mm.max():.1f}] mm/day")
lines.append(f"     Per-cell ES₉₅ median: {np.median(cell_es_mm):.1f} mm/day")
lines.append(f"     Note: With only 20 annual maxima, ES₉₅ is computed "
             f"from the top 1–2 values.")
lines.append(f"     This is inherently noisy but unbiased.")

lines.append("")
lines.append("=" * 85)

table_path = os.path.join(OUT_DIR, "cluster_risk_report.txt")
with open(table_path, "w") as f:
    f.write("\n".join(lines))
print(f"    ✓ {os.path.basename(table_path)}")

elapsed = time.time() - t0
print(f"\n{'=' * 65}")
print(f"  DONE in {elapsed:.0f}s")
print(f"  Maps:   {OUT_DIR}/")
print(f"  Report: {table_path}")
print(f"{'=' * 65}")

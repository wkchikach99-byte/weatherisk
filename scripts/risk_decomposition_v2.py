#!/usr/bin/env python3
"""Risk Model 2 — Empirical spatial risk decomposition.

Instead of using the parametric extremal coefficient θ (which saturates
near 2 at coarse resolution), this script measures spatial coherence
DIRECTLY from the 20-year block maxima data:

  1. Per-cell hazard intensity (ES₉₅ in mm/day)
  2. Empirical co-exceedance rate: when one cell is extreme, how many
     others in the cluster are simultaneously extreme?
  3. Empirical mean pairwise correlation of block maxima within cluster
  4. Event footprint: during the cluster's worst year, what fraction
     of cells are above their own 90th-percentile?
  5. GDP-weighted risk combining hazard, footprint, and exposure

Maps saved to docs/figures/risk_model_2/
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import run_cpc_pipeline, PipelineConfig
from weatherisk.risk import compute_var, compute_es

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

OUT_DIR = "docs/figures/risk_model_2"
os.makedirs(OUT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
#  1. Run pipeline
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("  RISK MODEL 2 — Empirical Spatial Risk Decomposition")
print("=" * 70)

t0 = time.time()

cfg = PipelineConfig(
    gdp_path="data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc",
    gdp_year=2015,
)
print("\n[1/7] Running CPC pipeline...")
result = run_cpc_pipeline(cfg, verbose=False)

bm = result["bm"]              # (20, 384) mm/day
labels = result["labels_lec"]
k = result["k_lec"]
lat = result["lat"]
lon = result["lon"]
gdp = result["gdp_per_cell"]
land_idx = result["land_idx"]
lats_1d = result["lats_1d"]
lons_1d = result["lons_1d"]
n_lat = result["n_lat"]
n_lon = result["n_lon"]
smoothed = result["smoothed"]

n_years, n_cells = bm.shape
clusters = sorted(np.unique(labels))
n_clusters = len(clusters)

print(f"    {n_years} years, {n_cells} cells, {n_clusters} LEC clusters")
print(f"    Block maxima range: [{bm.min():.1f}, {bm.max():.1f}] mm/day")


# ═══════════════════════════════════════════════════════════════
#  2. Per-cell hazard: ES₉₅ in mm/day
# ═══════════════════════════════════════════════════════════════

print("\n[2/7] Per-cell ES₉₅ (mm/day)...")

cell_es = np.array([compute_es(bm[:, c], 0.95) for c in range(n_cells)])
# Per-cell 90th percentile threshold (for co-exceedance)
cell_q90 = np.quantile(bm, 0.90, axis=0)  # (384,) per-cell 90th percentile

print(f"    ES₉₅: [{cell_es.min():.1f}, {cell_es.max():.1f}], "
      f"median={np.median(cell_es):.1f} mm/day")


# ═══════════════════════════════════════════════════════════════
#  3. Empirical spatial coherence metrics per cluster
# ═══════════════════════════════════════════════════════════════

print("\n[3/7] Computing empirical spatial coherence from actual data...")

cluster_stats = []

for cl in clusters:
    mask = labels == cl
    idx = np.where(mask)[0]
    n_cl = len(idx)

    bm_cl = bm[:, mask]  # (20, n_cl) block maxima for this cluster

    # --- Per-cluster spatial max ES in mm/day ---
    cmax = bm_cl.max(axis=1)  # spatial max per year
    es_spatial = compute_es(cmax, 0.95)
    var_spatial = compute_var(cmax, 0.95)

    # --- Mean per-cell ES ---
    mean_cell_es = cell_es[mask].mean()
    max_cell_es = cell_es[mask].max()

    # --- METRIC 1: Mean pairwise Spearman correlation of BM time series ---
    #     (How correlated are annual extremes between cells in this cluster?)
    if n_cl >= 2:
        from scipy.stats import spearmanr
        # Compute all pairwise Spearman correlations
        rho_matrix, _ = spearmanr(bm_cl)  # (n_cl, n_cl) if n_cl > 2
        if n_cl == 2:
            # spearmanr returns scalar for 2 columns
            mean_rho = float(rho_matrix)
        else:
            # Mean of upper triangle
            tri_idx = np.triu_indices(n_cl, k=1)
            mean_rho = rho_matrix[tri_idx].mean()
            min_rho = rho_matrix[tri_idx].min()
            max_rho = rho_matrix[tri_idx].max()
    else:
        mean_rho = 1.0  # single cell
        min_rho = max_rho = 1.0

    # --- METRIC 2: Co-exceedance rate ---
    #     For each year, count fraction of cells exceeding their own q90.
    #     Then average across years.
    #     = "When it's extreme somewhere in this cluster, how widespread is the
    #       extreme across the cluster?"
    exceed = bm_cl >= cell_q90[mask][None, :]  # (20, n_cl) boolean
    co_exceedance_per_year = exceed.mean(axis=1)  # fraction of cells in exceedance per year
    # Focus on the worst years (when at least one cell exceeds q90)
    any_exceed = exceed.any(axis=1)  # (20,) — years with at least one tail event
    if any_exceed.sum() > 0:
        # Mean fraction of cells exceeding during years with at least one exceedance
        co_exceedance = co_exceedance_per_year[any_exceed].mean()
        # Max co-exceedance across all years
        max_co_exceedance = co_exceedance_per_year.max()
    else:
        co_exceedance = 0.0
        max_co_exceedance = 0.0

    # --- METRIC 3: Event footprint ---
    #     In the year of the spatial maximum, what fraction of cells
    #     exceed their own MEDIAN? (= how much of the cluster "activates")
    worst_year = np.argmax(cmax)
    cell_median = np.median(bm_cl, axis=0)
    footprint = (bm_cl[worst_year, :] >= cell_median).mean()

    # Also: in the worst year, what's the mean normalised intensity?
    # (each cell's BM / its own ES₉₅, averaged across cells)
    cell_es_cl = cell_es[mask]
    if cell_es_cl.min() > 0:
        norm_intensity_worst = (bm_cl[worst_year, :] / cell_es_cl).mean()
    else:
        norm_intensity_worst = 0.0

    # --- METRIC 4: Simultaneous tail count ---
    #     Average number of cells simultaneously in their top-2 years
    top2_threshold = np.quantile(bm_cl, 0.90, axis=0)  # 90th percentile per cell
    simul_tail = (bm_cl >= top2_threshold[None, :]).sum(axis=1)  # cells in tail per year
    mean_simul = simul_tail.mean()
    max_simul = simul_tail.max()

    # --- GDP ---
    gdp_total = gdp[mask].sum() if gdp is not None else 0.0

    cluster_stats.append(dict(
        cluster=cl,
        n_cells=n_cl,
        es_spatial_mm=es_spatial,
        var_spatial_mm=var_spatial,
        mean_cell_es=mean_cell_es,
        max_cell_es=max_cell_es,
        mean_rho=mean_rho,
        co_exceedance=co_exceedance,
        max_co_exceedance=max_co_exceedance,
        footprint=footprint,
        norm_intensity_worst=norm_intensity_worst,
        mean_simul=mean_simul,
        max_simul=max_simul,
        gdp_total=gdp_total,
    ))

print("\n    Cluster-level empirical coherence:")
print(f"    {'Cl':>3s} {'N':>4s} {'ES_sp':>6s} {'ρ̄':>6s} "
      f"{'CoExc':>6s} {'Foot':>5s} {'SimTail':>7s} {'GDP':>8s}")
print("    " + "-" * 55)
for s in sorted(cluster_stats, key=lambda x: x["es_spatial_mm"], reverse=True):
    def _fmt(v):
        if v >= 1e12: return f"${v/1e12:.0f}T"
        if v >= 1e9: return f"${v/1e9:.0f}B"
        if v >= 1e6: return f"${v/1e6:.0f}M"
        return f"${v:,.0f}"
    print(f"    {s['cluster']:3d} {s['n_cells']:4d} {s['es_spatial_mm']:6.1f} "
          f"{s['mean_rho']:6.3f} {s['co_exceedance']:6.3f} "
          f"{s['footprint']:5.2f} {s['max_simul']:4.0f}/{s['n_cells']:<3d} "
          f"{_fmt(s['gdp_total']):>8s}")


# ═══════════════════════════════════════════════════════════════
#  4. Key diagnostic: Does ρ̄ vary across clusters?
# ═══════════════════════════════════════════════════════════════

print("\n[4/7] Coherence variation analysis...")

multi = [s for s in cluster_stats if s["n_cells"] >= 3]
rhos = np.array([s["mean_rho"] for s in multi])
coexcs = np.array([s["co_exceedance"] for s in multi])
footprints = np.array([s["footprint"] for s in multi])

print(f"\n    Mean Spearman ρ̄ (clusters ≥ 3 cells):")
print(f"      range: [{rhos.min():.3f}, {rhos.max():.3f}]")
print(f"      variation: {rhos.max()/max(rhos.min(), 1e-6):.1f}x")
print(f"    Co-exceedance rate:")
print(f"      range: [{coexcs.min():.3f}, {coexcs.max():.3f}]")
print(f"    Footprint (worst year):")
print(f"      range: [{footprints.min():.2f}, {footprints.max():.2f}]")

# Does coherence change the ranking vs pure marginal?
from scipy.stats import spearmanr as sp_corr
marg = np.array([s["mean_cell_es"] for s in multi])
agg = np.array([s["es_spatial_mm"] * (1 + s["mean_rho"]) / 2 for s in multi])
rho_rank, p_rank = sp_corr(marg, agg)
print(f"\n    Spearman rank-corr (marginal vs ρ-adjusted aggregate): "
      f"{rho_rank:.3f} (p={p_rank:.4f})")


# ═══════════════════════════════════════════════════════════════
#  5. Composite Risk Score
# ═══════════════════════════════════════════════════════════════

print("\n[5/7] Computing composite risk scores...")

# Risk score per cluster:
#   Hazard = spatial-max ES₉₅ (mm/day)
#   Coherence factor = (1 + ρ̄) / 2 ∈ [0.5, 1.0]
#     ρ̄ = 1 → coherence = 1.0 (worst case)
#     ρ̄ = 0 → coherence = 0.5 (independent)
#     ρ̄ < 0 → coherence < 0.5 (anti-correlated, risk-reducing)
#   Exposure = total GDP

for s in cluster_stats:
    s["coherence_factor"] = (1.0 + s["mean_rho"]) / 2.0
    s["hazard_coherence"] = s["es_spatial_mm"] * s["coherence_factor"]
    s["full_risk"] = s["hazard_coherence"] * s["gdp_total"]


# ═══════════════════════════════════════════════════════════════
#  6. Generate maps
# ═══════════════════════════════════════════════════════════════

print("\n[6/7] Generating maps...")

def to_grid(values):
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
extent = [cfg.lon_range[0] - 2, cfg.lon_range[1] + 2,
          cfg.lat_range[0] - 2, cfg.lat_range[1] + 2]

import cartopy.crs as ccrs
import cartopy.feature as cfeature


def make_ax(fig, pos):
    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax.add_feature(cfeature.LAND, facecolor="wheat", alpha=0.10)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    return ax


# ── Map 1: Per-cell ES₉₅ (mm/day) ──────────────────────────

fig = plt.figure(figsize=(11, 7))
ax = make_ax(fig, 111)
grid = to_grid(cell_es)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(mesh, ax=ax, shrink=0.7, label="ES₉₅ (mm/day)")
ax.set_title("(a) Hazard Intensity — Per-Cell Expected Shortfall\n"
             "95th-Percentile Tail Mean of Annual Block Maxima (2000–2019)",
             fontsize=12, fontweight="bold")
p = os.path.join(OUT_DIR, "map1_hazard_intensity.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── Map 2: Empirical Spearman ρ̄ per cluster ─────────────────

rho_per_cell = np.zeros(n_cells)
for s in cluster_stats:
    rho_per_cell[labels == s["cluster"]] = s["mean_rho"]

fig = plt.figure(figsize=(11, 7))
ax = make_ax(fig, 111)
grid = to_grid(rho_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="RdYlBu_r",
                     vmin=-0.2, vmax=0.8,
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7)
cb.set_label("Mean Spearman ρ̄ (within-cluster)")
cb.ax.annotate("High: extremes co-occur →",
               xy=(0.5, 0.98), fontsize=7, ha="center",
               xycoords="axes fraction", color="darkred")
cb.ax.annotate("← Low: independent extremes",
               xy=(0.5, 0.02), fontsize=7, ha="center",
               xycoords="axes fraction", color="darkblue")
ax.set_title("(b) Spatial Coherence — Empirical Correlation of Extremes\n"
             f"Mean Pairwise Spearman ρ of Annual Block Maxima per LEC Cluster (k={k})",
             fontsize=12, fontweight="bold")
p = os.path.join(OUT_DIR, "map2_empirical_coherence.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── Map 3: Co-exceedance rate per cluster ────────────────────

coexc_per_cell = np.zeros(n_cells)
for s in cluster_stats:
    coexc_per_cell[labels == s["cluster"]] = s["co_exceedance"]

fig = plt.figure(figsize=(11, 7))
ax = make_ax(fig, 111)
grid = to_grid(coexc_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="OrRd",
                     vmin=0.0, vmax=0.6,
                     transform=ccrs.PlateCarree(), shading="flat")
cb = plt.colorbar(mesh, ax=ax, shrink=0.7)
cb.set_label("Co-exceedance rate")
ax.set_title("(c) Simultaneous Extreme Footprint\n"
             "Fraction of Cluster Cells Exceeding 90th Percentile During Tail Events",
             fontsize=12, fontweight="bold")
p = os.path.join(OUT_DIR, "map3_coexceedance.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── Map 4: Composite hazard × coherence ─────────────────────

hc_per_cell = np.zeros(n_cells)
for s in cluster_stats:
    hc_per_cell[labels == s["cluster"]] = s["hazard_coherence"]

fig = plt.figure(figsize=(11, 7))
ax = make_ax(fig, 111)
grid = to_grid(hc_per_cell)
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="hot_r",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(mesh, ax=ax, shrink=0.7,
             label="ES₉₅ × Coherence Factor (mm/day)")
ax.set_title("(d) Coherence-Adjusted Hazard\n"
             "Spatial-Max ES₉₅ × (1+ρ̄)/2  — Higher When Extremes Are Widespread",
             fontsize=12, fontweight="bold")
p = os.path.join(OUT_DIR, "map4_coherence_hazard.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── Map 5: GDP-weighted full risk ────────────────────────────

risk_per_cell = np.zeros(n_cells)
for s in cluster_stats:
    mask = labels == s["cluster"]
    # Each cell: cluster hazard_coherence × cell GDP
    risk_per_cell[mask] = s["hazard_coherence"] * gdp[mask]

fig = plt.figure(figsize=(11, 7))
ax = make_ax(fig, 111)
grid = to_grid(np.log10(np.maximum(risk_per_cell, 1.0)))
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(mesh, ax=ax, shrink=0.7, label="log₁₀(Risk Score)")
ax.set_title("(e) Exposure-Weighted Precipitation Risk\n"
             "Coherence-Adjusted ES₉₅ × GDP PPP per Cell (2000–2019, GDP 2015)",
             fontsize=12, fontweight="bold")
p = os.path.join(OUT_DIR, "map5_gdp_risk.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ── Summary panel: 2×3 ──────────────────────────────────────

fig = plt.figure(figsize=(20, 13))

# (1) LEC clusters
ax1 = fig.add_subplot(2, 3, 1, projection=ccrs.PlateCarree())
ax1.set_extent(extent, crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax1.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
base_colors = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
cmap_cl = ListedColormap([base_colors[i % len(base_colors)] for i in range(k)])
norm_cl = BoundaryNorm(np.arange(0.5, k + 1.5), cmap_cl.N)
grid = to_grid(labels.astype(float))
masked = np.ma.masked_invalid(grid)
ax1.pcolormesh(lon_e, lat_e, masked, cmap=cmap_cl, norm=norm_cl,
               transform=ccrs.PlateCarree(), shading="flat")
ax1.set_title(f"(a) LEC Clusters (k={k})", fontsize=11, fontweight="bold")

# (2) Per-cell ES₉₅
ax2 = fig.add_subplot(2, 3, 2, projection=ccrs.PlateCarree())
ax2.set_extent(extent, crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax2.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
grid = to_grid(cell_es)
masked = np.ma.masked_invalid(grid)
m2 = ax2.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(m2, ax=ax2, shrink=0.55, label="mm/day")
ax2.set_title("(b) Per-Cell ES₉₅ (mm/day)", fontsize=11, fontweight="bold")

# (3) Empirical ρ̄
ax3 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
ax3.set_extent(extent, crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax3.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
grid = to_grid(rho_per_cell)
masked = np.ma.masked_invalid(grid)
m3 = ax3.pcolormesh(lon_e, lat_e, masked, cmap="RdYlBu_r",
                     vmin=-0.2, vmax=0.8,
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(m3, ax=ax3, shrink=0.55, label="Spearman ρ̄")
ax3.set_title("(c) Empirical Coherence (ρ̄)", fontsize=11, fontweight="bold")

# (4) Co-exceedance
ax4 = fig.add_subplot(2, 3, 4, projection=ccrs.PlateCarree())
ax4.set_extent(extent, crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax4.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
grid = to_grid(coexc_per_cell)
masked = np.ma.masked_invalid(grid)
m4 = ax4.pcolormesh(lon_e, lat_e, masked, cmap="OrRd",
                     vmin=0.0, vmax=0.6,
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(m4, ax=ax4, shrink=0.55, label="Co-exceedance rate")
ax4.set_title("(d) Co-Exceedance Rate", fontsize=11, fontweight="bold")

# (5) Hazard × Coherence
ax5 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
ax5.set_extent(extent, crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax5.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax5.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
grid = to_grid(hc_per_cell)
masked = np.ma.masked_invalid(grid)
m5 = ax5.pcolormesh(lon_e, lat_e, masked, cmap="hot_r",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(m5, ax=ax5, shrink=0.55, label="ES₉₅ × (1+ρ̄)/2")
ax5.set_title("(e) Hazard × Coherence", fontsize=11, fontweight="bold")

# (6) GDP risk
ax6 = fig.add_subplot(2, 3, 6, projection=ccrs.PlateCarree())
ax6.set_extent(extent, crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
ax6.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax6.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
grid = to_grid(np.log10(np.maximum(risk_per_cell, 1.0)))
masked = np.ma.masked_invalid(grid)
m6 = ax6.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                     transform=ccrs.PlateCarree(), shading="flat")
plt.colorbar(m6, ax=ax6, shrink=0.55, label="log₁₀(Risk)")
ax6.set_title("(f) GDP-Weighted Risk", fontsize=11, fontweight="bold")

plt.suptitle("Spatial Risk Decomposition — CPC Precipitation Extremes (2000–2019)\n"
             "LEC Clusters + Empirical Coherence + GDP Exposure",
             fontsize=14, fontweight="bold", y=0.99)
plt.tight_layout(rect=[0, 0, 1, 0.95])

p = os.path.join(OUT_DIR, "summary_panel_6maps.png")
fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
print(f"    ✓ {os.path.basename(p)}")


# ═══════════════════════════════════════════════════════════════
#  7. Full summary table + findings
# ═══════════════════════════════════════════════════════════════

print("\n[7/7] Summary table\n")

def _fmt(v):
    if v >= 1e12: return f"${v/1e12:.1f}T"
    if v >= 1e9: return f"${v/1e9:.1f}B"
    if v >= 1e6: return f"${v/1e6:.1f}M"
    if v >= 1e3: return f"${v/1e3:.0f}K"
    return f"${v:,.0f}"

sorted_stats = sorted(cluster_stats,
    key=lambda s: s["hazard_coherence"], reverse=True)

header = (f"{'Cl':>3s} {'N':>4s} {'ES_sp':>6s} {'CellES':>6s} "
          f"{'ρ̄':>6s} {'CoExc':>5s} {'Foot':>5s} {'Coh.F':>5s} "
          f"{'ES×Coh':>7s} {'GDP':>8s} {'FullRisk':>10s}")
print(header)
print("-" * len(header))

for s in sorted_stats:
    print(
        f"{s['cluster']:3d} {s['n_cells']:4d} "
        f"{s['es_spatial_mm']:6.1f} {s['mean_cell_es']:6.1f} "
        f"{s['mean_rho']:6.3f} {s['co_exceedance']:5.3f} "
        f"{s['footprint']:5.2f} {s['coherence_factor']:5.3f} "
        f"{s['hazard_coherence']:7.1f} "
        f"{_fmt(s['gdp_total']):>8s} {_fmt(s['full_risk']):>10s}"
    )

# Key findings
print(f"\n{'='*70}")
print("  KEY FINDINGS — RISK MODEL 2")
print(f"{'='*70}")

multi = [s for s in cluster_stats if s["n_cells"] >= 3]
rhos = [s["mean_rho"] for s in multi]
coexcs = [s["co_exceedance"] for s in multi]

print(f"\n  Empirical Spearman ρ̄ across clusters (≥3 cells):")
print(f"    Range: [{min(rhos):.3f}, {max(rhos):.3f}]")
print(f"    → Factor of variation: {max(rhos)/max(min(rhos), 0.001):.1f}x")

print(f"\n  Co-exceedance rate:")
print(f"    Range: [{min(coexcs):.3f}, {max(coexcs):.3f}]")

# Most coherent
most_coh = max(multi, key=lambda s: s["mean_rho"])
print(f"\n  Most spatially coherent cluster: #{most_coh['cluster']} "
      f"({most_coh['n_cells']} cells, ρ̄={most_coh['mean_rho']:.3f})")

# Least coherent
least_coh = min(multi, key=lambda s: s["mean_rho"])
print(f"  Least coherent cluster: #{least_coh['cluster']} "
      f"({least_coh['n_cells']} cells, ρ̄={least_coh['mean_rho']:.3f})")

# Risk ranking impact
multi_sorted_marg = sorted(multi, key=lambda s: s["mean_cell_es"], reverse=True)
multi_sorted_coh = sorted(multi, key=lambda s: s["hazard_coherence"], reverse=True)
print(f"\n  Risk ranking comparison (marginal vs coherence-adjusted):")
for i, s in enumerate(multi_sorted_coh[:5]):
    r_marg = next(j+1 for j, m in enumerate(multi_sorted_marg)
                  if m["cluster"] == s["cluster"])
    print(f"    #{i+1} by coherence-adjusted: Cl {s['cluster']} "
          f"(marginal rank #{r_marg})")

elapsed = time.time() - t0
print(f"\n  Total time: {elapsed:.0f} s ({elapsed/60:.1f} min)")
print(f"  Maps saved to {OUT_DIR}/")

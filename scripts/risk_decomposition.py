#!/usr/bin/env python3
"""Spatial risk decomposition using LEC clusters + extremal coefficients.

Computes three layers of risk information:
  1. Per-cell hazard intensity (ES₉₅ in mm/day)
  2. Within-cluster spatial coherence (mean extremal coefficient θ̄)
  3. GDP exposure per cluster

Generates publication-quality maps and a summary table.
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import run_cpc_pipeline, PipelineConfig
from weatherisk.risk import compute_var, compute_es
from weatherisk.covariance import (
    cov_fkt_2d_nonstat2,
    build_nonstat_cov_matrix,
    cov_to_ec,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


# ═══════════════════════════════════════════════════════════════
#  1. Run pipeline
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("  SPATIAL RISK DECOMPOSITION — LEC Clusters + Extremal Coefficients")
print("=" * 70)

t0 = time.time()

cfg = PipelineConfig(
    gdp_path="data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc",
    gdp_year=2015,
)
print("\n[1/6] Running CPC pipeline (Steps 1–8)...")
result = run_cpc_pipeline(cfg, verbose=False)

bm = result["bm"]              # (20, 384) mm/day block maxima
frechet = result["frechet"]     # (20, 384) Fréchet
labels = result["labels_lec"]   # (384,) cluster assignments
k = result["k_lec"]
smoothed = result["smoothed"]   # (384, 3) per-cell (a, b, γ)
incl_lec = result["incl_lec"]   # dict: cluster → (a, b, γ)
gev_p = result["gev_params"]    # (384, 3) per-cell GEV (μ, σ, ξ)
lat = result["lat"]
lon = result["lon"]
gdp = result["gdp_per_cell"]    # (384,) GDP PPP per cell
land_idx = result["land_idx"]
lats_1d = result["lats_1d"]
lons_1d = result["lons_1d"]
n_lat = result["n_lat"]
n_lon = result["n_lon"]
coords = np.column_stack([
    -5.0 + 10.0 * (lat - lat.min()) / max(lat.max() - lat.min(), 1e-6),
    -5.0 + 10.0 * (lon - lon.min()) / max(lon.max() - lon.min(), 1e-6),
])

n_years, n_cells = bm.shape
clusters = sorted(np.unique(labels))
n_clusters = len(clusters)

print(f"    {n_years} years, {n_cells} cells, {n_clusters} LEC clusters")
print(f"    Pipeline took {time.time() - t0:.0f} s")


# ═══════════════════════════════════════════════════════════════
#  2. Per-cell hazard: ES₉₅ in mm/day
# ═══════════════════════════════════════════════════════════════

print("\n[2/6] Computing per-cell ES₉₅ (mm/day)...")

cell_es = np.array([compute_es(bm[:, c], 0.95) for c in range(n_cells)])

print(f"    Per-cell ES₉₅: min={cell_es.min():.1f}, "
      f"max={cell_es.max():.1f}, median={np.median(cell_es):.1f} mm/day")


# ═══════════════════════════════════════════════════════════════
#  3. Within-cluster extremal coefficients
# ═══════════════════════════════════════════════════════════════

print("\n[3/6] Computing pairwise extremal coefficients from fitted model...")

# Use smoothed per-cell parameters to build the full covariance matrix
a_flat = smoothed[:, 0]
b_flat = smoothed[:, 1]
g_flat = smoothed[:, 2]
X_flat = coords[:, 1]  # normalised lon
Y_flat = coords[:, 0]  # normalised lat

# Build full non-stationary covariance matrix
cov_matrix = build_nonstat_cov_matrix(X_flat, Y_flat, cfg.alpha, a_flat, b_flat, g_flat)

# Convert to extremal coefficients
ec_matrix = np.zeros_like(cov_matrix)
for i in range(n_cells):
    for j in range(i + 1, n_cells):
        ec_matrix[i, j] = cov_to_ec(cfg.df, cov_matrix[i, j])
        ec_matrix[j, i] = ec_matrix[i, j]
# Diagonal: θ(s,s) = 1 (perfect dependence with self)
np.fill_diagonal(ec_matrix, 1.0)

print(f"    EC matrix: shape {ec_matrix.shape}")
print(f"    θ range (off-diag): [{ec_matrix[np.triu_indices(n_cells, k=1)].min():.3f}, "
      f"{ec_matrix[np.triu_indices(n_cells, k=1)].max():.3f}]")
print(f"    θ median (off-diag): {np.median(ec_matrix[np.triu_indices(n_cells, k=1)]):.3f}")

# Within-cluster mean θ
cluster_stats = []
for cl in clusters:
    mask = labels == cl
    idx = np.where(mask)[0]
    n_cl = len(idx)

    # Per-cluster spatial max ES in mm/day
    cmax_mm = bm[:, mask].max(axis=1)
    es_spatial = compute_es(cmax_mm, 0.95)
    var_spatial = compute_var(cmax_mm, 0.95)

    # Mean per-cell ES in mm/day
    mean_cell_es = cell_es[mask].mean()
    max_cell_es = cell_es[mask].max()

    # Within-cluster mean θ
    if n_cl >= 2:
        # Extract sub-matrix of within-cluster extremal coefficients
        ec_sub = ec_matrix[np.ix_(idx, idx)]
        # Mean of upper triangle (excluding diagonal)
        tri_vals = ec_sub[np.triu_indices(n_cl, k=1)]
        theta_mean = tri_vals.mean()
        theta_min = tri_vals.min()
        theta_max = tri_vals.max()
    else:
        theta_mean = 1.0  # single cell: perfect self-dependence
        theta_min = theta_max = 1.0

    # Effective number of independent sites
    # n_eff ≈ n_cl * (θ̄ - 1) for θ̄ ∈ [1, 2]
    # When θ̄ = 1: all perfectly dependent → n_eff = 0 (acts as 1 site)
    # When θ̄ = 2: all independent → n_eff = n_cl
    # Simpler: n_eff = n_cl * (θ̄ / 2), where θ̄/2 is the "independence fraction"
    n_eff = n_cl * (theta_mean / 2.0)

    # Coherence score: 1/θ̄ ∈ [0.5, 1.0]
    # Higher = more coherent = more dangerous for aggregate risk
    coherence = 1.0 / theta_mean

    # GDP
    gdp_total = gdp[mask].sum() if gdp is not None else 0.0

    cluster_stats.append(dict(
        cluster=cl,
        n_cells=n_cl,
        es_spatial_mm=es_spatial,
        var_spatial_mm=var_spatial,
        mean_cell_es=mean_cell_es,
        max_cell_es=max_cell_es,
        theta_mean=theta_mean,
        theta_min=theta_min,
        theta_max=theta_max,
        n_eff=n_eff,
        coherence=coherence,
        gdp_total=gdp_total,
    ))

print("\n    Within-cluster θ̄ by cluster:")
for s in cluster_stats:
    print(f"    Cl {s['cluster']:2d} ({s['n_cells']:3d} cells): "
          f"θ̄ = {s['theta_mean']:.3f}  "
          f"[{s['theta_min']:.3f}, {s['theta_max']:.3f}]  "
          f"n_eff = {s['n_eff']:.1f}")


# ═══════════════════════════════════════════════════════════════
#  4. Risk ranking comparison
# ═══════════════════════════════════════════════════════════════

print("\n[4/6] Risk ranking comparison...")

# Only consider clusters with >= 2 cells for meaningful comparison
multi = [s for s in cluster_stats if s["n_cells"] >= 2]

# Rank by marginal intensity (mean per-cell ES)
rank_marginal = sorted(multi, key=lambda s: s["mean_cell_es"], reverse=True)
# Rank by spatial aggregate (ES_spatial × coherence)
rank_aggregate = sorted(multi,
    key=lambda s: s["es_spatial_mm"] * s["coherence"], reverse=True)
# Rank by ES_spatial × coherence × GDP
rank_full = sorted(multi,
    key=lambda s: s["es_spatial_mm"] * s["coherence"] * s["gdp_total"],
    reverse=True)

print("\n    ┌─────┬───────┬─────────────┬──────────────────┬──────────────────────┐")
print("    │  Cl │ Cells │ Mean ES     │ Spatial ES × 1/θ̄ │ Full (×GDP)          │")
print("    │     │       │ (mm/day)    │ (mm/day)         │                      │")
print("    ├─────┼───────┼─────────────┼──────────────────┼──────────────────────┤")
for i, s in enumerate(rank_full[:15]):
    marg_rank = next(j+1 for j, m in enumerate(rank_marginal) if m["cluster"] == s["cluster"])
    agg_rank = next(j+1 for j, m in enumerate(rank_aggregate) if m["cluster"] == s["cluster"])
    full_rank = i + 1

    def _fmt_gdp(v):
        if v >= 1e12: return f"${v/1e12:.1f}T"
        if v >= 1e9: return f"${v/1e9:.1f}B"
        if v >= 1e6: return f"${v/1e6:.1f}M"
        return f"${v:,.0f}"

    risk_score = s["es_spatial_mm"] * s["coherence"] * s["gdp_total"]
    print(f"    │ {s['cluster']:3d} │ {s['n_cells']:5d} │ "
          f"{s['mean_cell_es']:6.1f} (#{marg_rank:<2d}) │ "
          f"{s['es_spatial_mm'] * s['coherence']:7.1f} (#{agg_rank:<2d})  │ "
          f"{_fmt_gdp(risk_score):>10s} (#{full_rank})     │")
print("    └─────┴───────┴─────────────┴──────────────────┴──────────────────────┘")


# === Does θ̄ actually vary enough to change rankings? ===
thetas = np.array([s["theta_mean"] for s in multi])
print(f"\n    θ̄ across multi-cell clusters: "
      f"min={thetas.min():.3f}, max={thetas.max():.3f}, "
      f"std={thetas.std():.3f}, CV={thetas.std()/thetas.mean():.1%}")
print(f"    → θ̄ varies by a factor of {thetas.max()/thetas.min():.2f}x")

# Spearman rank correlation between marginal rank and aggregate rank
from scipy.stats import spearmanr
marg_vals = np.array([s["mean_cell_es"] for s in multi])
agg_vals = np.array([s["es_spatial_mm"] * s["coherence"] for s in multi])
rho, pval = spearmanr(marg_vals, agg_vals)
print(f"    Spearman ρ (marginal vs aggregate rank): {rho:.3f}  (p={pval:.4f})")
if abs(rho) < 0.8:
    print("    → GOOD: Spatial dependence is changing the risk ranking substantially.")
elif abs(rho) < 0.95:
    print("    → MODERATE: Some reranking due to spatial dependence.")
else:
    print("    → WEAK: Rankings are very similar — spatial dependence adds little.")


# ═══════════════════════════════════════════════════════════════
#  5. Generate maps
# ═══════════════════════════════════════════════════════════════

print("\n[5/6] Generating maps...")

out_dir = "docs/figures"
os.makedirs(out_dir, exist_ok=True)

# Grid reconstruction helper
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

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("    WARNING: cartopy not available, skipping geographic maps")


if HAS_CARTOPY:

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
    ax.set_title("Hazard Intensity — Per-Cell Expected Shortfall (ES₉₅)\n"
                 "Annual Precipitation Block Maxima, 2000–2019",
                 fontsize=12, fontweight="bold")
    p = os.path.join(out_dir, "risk_hazard_intensity.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ {os.path.basename(p)}")

    # ── Map 2: Within-cluster mean θ̄ ────────────────────────────

    # Assign θ̄ to each cell based on its cluster
    theta_per_cell = np.zeros(n_cells)
    for s in cluster_stats:
        theta_per_cell[labels == s["cluster"]] = s["theta_mean"]

    fig = plt.figure(figsize=(11, 7))
    ax = make_ax(fig, 111)
    grid = to_grid(theta_per_cell)
    masked = np.ma.masked_invalid(grid)
    # Use reversed colormap: low θ = high coherence = warm colors (dangerous)
    mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="RdYlGn",
                         vmin=1.0, vmax=2.0,
                         transform=ccrs.PlateCarree(), shading="flat")
    cb = plt.colorbar(mesh, ax=ax, shrink=0.7)
    cb.set_label("Mean Extremal Coefficient θ̄")
    cb.ax.annotate("← More coherent\n   (simultaneous extremes)",
                   xy=(0.5, 0.02), fontsize=7, ha="center",
                   xycoords="axes fraction", color="darkred")
    cb.ax.annotate("More independent →\n(diversified risk)",
                   xy=(0.5, 0.95), fontsize=7, ha="center",
                   xycoords="axes fraction", color="darkgreen")
    ax.set_title("Spatial Coherence — Within-Cluster Extremal Coefficient (θ̄)\n"
                 f"LEC Clusters (k={k}), max-stable process model, 2000–2019",
                 fontsize=12, fontweight="bold")
    p = os.path.join(out_dir, "risk_spatial_coherence.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ {os.path.basename(p)}")

    # ── Map 3: Per-cluster aggregate risk score ─────────────────

    # Score = ES_spatial × coherence (no GDP, pure hazard+dependence)
    agg_score_per_cell = np.zeros(n_cells)
    for s in cluster_stats:
        agg_score_per_cell[labels == s["cluster"]] = (
            s["es_spatial_mm"] * s["coherence"]
        )

    fig = plt.figure(figsize=(11, 7))
    ax = make_ax(fig, 111)
    grid = to_grid(agg_score_per_cell)
    masked = np.ma.masked_invalid(grid)
    mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="hot_r",
                         transform=ccrs.PlateCarree(), shading="flat")
    plt.colorbar(mesh, ax=ax, shrink=0.7,
                 label="Aggregate Hazard Score (ES₉₅ × 1/θ̄)")
    ax.set_title("Aggregate Spatial Risk Score\n"
                 "Spatial-Max ES₉₅ (mm/day) × Coherence Factor (1/θ̄)",
                 fontsize=12, fontweight="bold")
    p = os.path.join(out_dir, "risk_aggregate_score.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ {os.path.basename(p)}")

    # ── Map 4: GDP-weighted risk ────────────────────────────────

    if gdp is not None:
        gdp_risk_per_cell = np.zeros(n_cells)
        for s in cluster_stats:
            mask = labels == s["cluster"]
            # Each cell gets: (cluster aggregate hazard score) × (cell GDP)
            gdp_risk_per_cell[mask] = (
                s["es_spatial_mm"] * s["coherence"] * gdp[mask]
            )

        fig = plt.figure(figsize=(11, 7))
        ax = make_ax(fig, 111)
        grid = to_grid(np.log10(np.maximum(gdp_risk_per_cell, 1.0)))
        masked = np.ma.masked_invalid(grid)
        mesh = ax.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                             transform=ccrs.PlateCarree(), shading="flat")
        cb = plt.colorbar(mesh, ax=ax, shrink=0.7,
                          label="log₁₀(Risk Score)")
        ax.set_title("Economic Exposure-Weighted Precipitation Risk\n"
                     "ES₉₅ × Coherence (1/θ̄) × GDP PPP per Cell  (2000–2019, GDP 2015)",
                     fontsize=12, fontweight="bold")
        p = os.path.join(out_dir, "risk_gdp_weighted.png")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"    ✓ {os.path.basename(p)}")

    # ── Summary panel: 2×2 ──────────────────────────────────────

    fig = plt.figure(figsize=(18, 14))

    # (1) LEC clusters
    ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
    ax1.set_extent(extent, crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    from matplotlib.colors import ListedColormap, BoundaryNorm
    base_colors = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    cmap_cl = ListedColormap([base_colors[i % len(base_colors)] for i in range(k)])
    norm_cl = BoundaryNorm(np.arange(0.5, k + 1.5), cmap_cl.N)
    grid = to_grid(labels.astype(float))
    masked = np.ma.masked_invalid(grid)
    ax1.pcolormesh(lon_e, lat_e, masked, cmap=cmap_cl, norm=norm_cl,
                   transform=ccrs.PlateCarree(), shading="flat")
    ax1.set_title(f"(a) LEC Spatial Dependence Clusters (k={k})",
                  fontsize=11, fontweight="bold")

    # (2) Per-cell ES₉₅
    ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    grid = to_grid(cell_es)
    masked = np.ma.masked_invalid(grid)
    m2 = ax2.pcolormesh(lon_e, lat_e, masked, cmap="YlOrRd",
                        transform=ccrs.PlateCarree(), shading="flat")
    plt.colorbar(m2, ax=ax2, shrink=0.6, label="mm/day")
    ax2.set_title("(b) Hazard Intensity — Per-Cell ES₉₅ (mm/day)",
                  fontsize=11, fontweight="bold")

    # (3) Within-cluster θ̄
    ax3 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
    ax3.set_extent(extent, crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax3.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    grid = to_grid(theta_per_cell)
    masked = np.ma.masked_invalid(grid)
    m3 = ax3.pcolormesh(lon_e, lat_e, masked, cmap="RdYlGn",
                        vmin=1.0, vmax=2.0,
                        transform=ccrs.PlateCarree(), shading="flat")
    cb3 = plt.colorbar(m3, ax=ax3, shrink=0.6, label="θ̄")
    ax3.set_title("(c) Spatial Coherence — Extremal Coefficient θ̄",
                  fontsize=11, fontweight="bold")

    # (4) Aggregate score
    ax4 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
    ax4.set_extent(extent, crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax4.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax4.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    grid = to_grid(agg_score_per_cell)
    masked = np.ma.masked_invalid(grid)
    m4 = ax4.pcolormesh(lon_e, lat_e, masked, cmap="hot_r",
                        transform=ccrs.PlateCarree(), shading="flat")
    plt.colorbar(m4, ax=ax4, shrink=0.6, label="ES₉₅ × 1/θ̄")
    ax4.set_title("(d) Aggregate Risk — Hazard × Coherence",
                  fontsize=11, fontweight="bold")

    plt.suptitle("Spatial Risk Decomposition — CPC Precipitation Extremes (2000–2019)\n"
                 "LEC Clustering with Extremal Coefficient Dependence Analysis",
                 fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    p = os.path.join(out_dir, "risk_decomposition_panel.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ {os.path.basename(p)}")


# ═══════════════════════════════════════════════════════════════
#  6. Full summary table
# ═══════════════════════════════════════════════════════════════

print("\n[6/6] Full cluster summary table\n")

def _fmt(v):
    if v >= 1e12: return f"${v/1e12:.1f}T"
    if v >= 1e9: return f"${v/1e9:.1f}B"
    if v >= 1e6: return f"${v/1e6:.1f}M"
    if v >= 1e3: return f"${v/1e3:.0f}K"
    return f"${v:,.0f}"

# Sort by aggregate hazard score (ES × coherence), descending
sorted_stats = sorted(cluster_stats,
    key=lambda s: s["es_spatial_mm"] * s["coherence"], reverse=True)

header = (
    f"{'Cl':>3s} {'Cells':>5s} "
    f"{'ES_sp':>7s} {'θ̄':>6s} {'θ_min':>6s} {'θ_max':>6s} "
    f"{'n_eff':>6s} {'1/θ̄':>5s} "
    f"{'ES×1/θ̄':>8s} {'CellES':>7s} "
    f"{'GDP':>10s} {'Full Risk':>12s}"
)
print(header)
print("-" * len(header))

for s in sorted_stats:
    agg_haz = s["es_spatial_mm"] * s["coherence"]
    full_risk = agg_haz * s["gdp_total"]
    print(
        f"{s['cluster']:3d} {s['n_cells']:5d} "
        f"{s['es_spatial_mm']:7.1f} {s['theta_mean']:6.3f} "
        f"{s['theta_min']:6.3f} {s['theta_max']:6.3f} "
        f"{s['n_eff']:6.1f} {s['coherence']:5.3f} "
        f"{agg_haz:8.1f} {s['mean_cell_es']:7.1f} "
        f"{_fmt(s['gdp_total']):>10s} {_fmt(full_risk):>12s}"
    )

# Key finding
print(f"\n{'='*70}")
print("  KEY FINDINGS")
print(f"{'='*70}")

thetas_all = [s["theta_mean"] for s in cluster_stats if s["n_cells"] >= 2]
print(f"\n  θ̄ range across clusters: "
      f"[{min(thetas_all):.3f}, {max(thetas_all):.3f}]")
print(f"  → Factor of variation: {max(thetas_all)/min(thetas_all):.2f}x")

# Most coherent cluster
most_coherent = min(
    [s for s in cluster_stats if s["n_cells"] >= 2],
    key=lambda s: s["theta_mean"])
print(f"\n  Most spatially coherent cluster: #{most_coherent['cluster']} "
      f"({most_coherent['n_cells']} cells, θ̄={most_coherent['theta_mean']:.3f})")
print(f"    → Extremes in this region tend to occur simultaneously")

# Most independent cluster
most_indep = max(
    [s for s in cluster_stats if s["n_cells"] >= 2],
    key=lambda s: s["theta_mean"])
print(f"\n  Most independent cluster: #{most_indep['cluster']} "
      f"({most_indep['n_cells']} cells, θ̄={most_indep['theta_mean']:.3f})")
print(f"    → Extremes in this region are more localised/diversified")

# Biggest risk reranking
if len(multi) > 5:
    marg_order = [s["cluster"] for s in sorted(multi,
        key=lambda s: s["mean_cell_es"], reverse=True)]
    agg_order = [s["cluster"] for s in sorted(multi,
        key=lambda s: s["es_spatial_mm"] * s["coherence"], reverse=True)]
    # Find biggest rank change
    max_change = 0
    change_cl = None
    for cl in marg_order:
        r1 = marg_order.index(cl) + 1
        r2 = agg_order.index(cl) + 1
        if abs(r1 - r2) > max_change:
            max_change = abs(r1 - r2)
            change_cl = cl
    if change_cl:
        r1 = marg_order.index(change_cl) + 1
        r2 = agg_order.index(change_cl) + 1
        print(f"\n  Largest rank change from spatial dependence: "
              f"Cluster #{change_cl}")
        print(f"    Marginal rank: #{r1} → Aggregate rank: #{r2} "
              f"(moved {abs(r1-r2)} places)")

elapsed = time.time() - t0
print(f"\n  Total time: {elapsed:.0f} s ({elapsed/60:.1f} min)")
print(f"\n  Maps saved to {out_dir}/")
print(f"    - risk_hazard_intensity.png")
print(f"    - risk_spatial_coherence.png")
print(f"    - risk_aggregate_score.png")
if gdp is not None:
    print(f"    - risk_gdp_weighted.png")
print(f"    - risk_decomposition_panel.png")

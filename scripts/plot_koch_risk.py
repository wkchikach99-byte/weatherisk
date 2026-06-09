#!/usr/bin/env python3
"""Produce Koch (2016) risk plots — E=1 pure hazard baseline.

Outputs (saved to risk/outputs/):
  1. temporal_loss_LEC.pdf/png  — L_N time series for LEC clusters
  2. temporal_loss_EDC.pdf/png  — L_N time series for EDC clusters
  3. var_map_LEC.pdf/png        — VaR world map (LEC)
  4. var_map_EDC.pdf/png        — VaR world map (EDC)
  5. es_map_LEC.pdf/png         — ES world map (LEC)
  6. es_map_EDC.pdf/png         — ES world map (EDC)
  7. threshold_sensitivity.pdf/png — VaR & ES vs p for top-k clusters

Usage:
    python scripts/plot_koch_risk.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from weatherisk.koch_risk import (
    compute_koch_risk,
    compute_cluster_loss,
    compute_var,
    compute_es,
    compute_trend,
    compute_cross_cluster_dependence,
    _cosine_weights,
    _frechet_threshold,
)

# ── paths ────────────────────────────────────────────────────────────
DATA_PATH = os.path.join(ROOT, "risk", "inputs", "pipeline_results.npz")
OUT_DIR   = os.path.join(ROOT, "risk", "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

DPI = 200

# ── helpers ──────────────────────────────────────────────────────────
def make_map_ax(fig, pos=111, projection=None):
    """GeoAxes with coastlines and subtle gridlines."""
    if projection is None:
        projection = ccrs.Robinson()
    ax = fig.add_subplot(pos, projection=projection)
    ax.set_global()
    ax.coastlines(linewidth=0.4, color="0.3")
    ax.gridlines(draw_labels=False, linewidth=0.15, color="grey", alpha=0.5)
    return ax


def _build_grid_field(
    values_per_cell: np.ndarray,
    valid_idx: np.ndarray,
    n_lat: int,
    n_lon: int,
) -> np.ndarray:
    """Map 1-D cell values back to a (n_lat, n_lon) grid; NaN elsewhere."""
    grid = np.full(n_lat * n_lon, np.nan)
    grid[valid_idx.astype(int)] = values_per_cell
    return grid.reshape(n_lat, n_lon)


# ═════════════════════════════════════════════════════════════════════
#  PLOT 1 — Temporal L_N time series
# ═════════════════════════════════════════════════════════════════════
def plot_temporal_loss(results, cname, top_k=8):
    """Line chart of L_N(t) for the top-k clusters by ES, with VaR line."""
    res = results[cname]
    # Sort by ES descending, take top_k
    sorted_cl = sorted(res.clusters, key=lambda c: c["ES"], reverse=True)[:top_k]

    fig, axes = plt.subplots(
        (top_k + 1) // 2, 2,
        figsize=(14, 2.5 * ((top_k + 1) // 2)),
        sharex=True,
    )
    axes = axes.flatten()

    for i, cl in enumerate(sorted_cl):
        ax = axes[i]
        years = res.years
        L = cl["L_N"]
        var_q = cl["VaR"]

        ax.plot(years, L, lw=0.7, color="steelblue", alpha=0.85)
        ax.axhline(var_q, ls="--", lw=0.8, color="firebrick", alpha=0.7,
                    label=f"VaR$_{{0.95}}$={var_q:.3f}")
        ax.fill_between(years, L, var_q, where=(L >= var_q),
                        color="firebrick", alpha=0.15)
        ax.set_ylabel("$L_N$", fontsize=9)
        ax.set_title(
            f"Cluster {cl['cluster']} ({cl['n_cells']} cells)   "
            f"ES$_{{0.95}}$={cl['ES']:.3f}",
            fontsize=9,
        )
        ax.legend(fontsize=7, loc="upper right")
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"Koch Loss $L_N$ — {cname} clusters (top {top_k} by ES$_{{0.95}}$, "
        f"$p$={res.p}, $q$={res.q})",
        fontsize=12,
        y=1.01,
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"temporal_loss_{cname}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ temporal_loss_{cname}")


# ═════════════════════════════════════════════════════════════════════
#  PLOT 2 & 3 — VaR / ES world maps
# ═════════════════════════════════════════════════════════════════════
def _plot_risk_map(
    metric_name: str,
    cname: str,
    results,
    valid_idx: np.ndarray,
    labels: np.ndarray,
    lats_1d: np.ndarray,
    lons_1d: np.ndarray,
):
    """Colour every cell by its cluster's VaR or ES."""
    res = results[cname]
    n_lat, n_lon = len(lats_1d), len(lons_1d)

    # Map cluster -> metric value
    cl_metric = {}
    for cl in res.clusters:
        cl_metric[cl["cluster"]] = cl[metric_name]

    # Assign metric to each cell
    cell_values = np.array([cl_metric[labels[i]] for i in range(len(labels))])
    grid = _build_grid_field(cell_values, valid_idx, n_lat, n_lon)

    fig = plt.figure(figsize=(14, 7))
    ax = make_map_ax(fig)

    v = cell_values[np.isfinite(cell_values)]
    vmin, vmax = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    if vmin == vmax:
        vmin, vmax = 0, 1

    im = ax.pcolormesh(
        lons_1d, lats_1d, grid,
        transform=ccrs.PlateCarree(),
        cmap="YlOrRd",
        vmin=vmin, vmax=vmax,
        shading="auto",
        rasterized=True,
    )
    label_map = {"VaR": f"VaR$_{{0.95}}$", "ES": f"ES$_{{0.95}}$"}
    cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                      pad=0.04, shrink=0.7, aspect=35)
    cb.set_label(f"{label_map[metric_name]} — Koch indicator fraction", fontsize=12)
    ax.set_title(
        f"{label_map[metric_name]} per {cname} cluster  "
        f"($p$={res.p}, $q$={res.q})",
        fontsize=14, pad=12,
    )

    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{metric_name.lower()}_map_{cname}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {metric_name.lower()}_map_{cname}")


def plot_var_map(cname, results, valid_idx, labels, lats_1d, lons_1d):
    _plot_risk_map("VaR", cname, results, valid_idx, labels, lats_1d, lons_1d)


def plot_es_map(cname, results, valid_idx, labels, lats_1d, lons_1d):
    _plot_risk_map("ES", cname, results, valid_idx, labels, lats_1d, lons_1d)


# ═════════════════════════════════════════════════════════════════════
#  PLOT 4 — Threshold sensitivity: VaR & ES vs p
# ═════════════════════════════════════════════════════════════════════
def plot_threshold_sensitivity(top_k=5):
    """Plot VaR_0.95 and ES_0.95 as a function of p for top clusters."""
    d = np.load(DATA_PATH)
    frechet = d["frechet"]
    labels_lec = d["labels_lec"]
    lats = d["lats"]
    valid_idx = d["valid_idx"]
    n_lon = int(d["lons"].shape[0])
    cos_weights = _cosine_weights(lats, valid_idx, n_lon)

    # Identify top-k LEC clusters by ES at p=0.95
    base_results = compute_koch_risk(DATA_PATH, p=0.95, q=0.95, clusterings=["LEC"])
    top_ids = [
        c["cluster"]
        for c in sorted(base_results["LEC"].clusters, key=lambda x: x["ES"], reverse=True)[:top_k]
    ]

    p_values = np.arange(0.80, 0.991, 0.01)
    q = 0.95

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    colours = matplotlib.colormaps.get_cmap("tab10").resampled(top_k)

    for idx, cl_id in enumerate(top_ids):
        mask = labels_lec == cl_id
        n_cells = int(mask.sum())
        vars_, ess_ = [], []
        for p in p_values:
            u_p = _frechet_threshold(p)
            L_N = (frechet[:, mask] > u_p).astype(np.float64) @ (
                cos_weights[mask] / cos_weights[mask].sum()
            )
            vars_.append(compute_var(L_N, q))
            ess_.append(compute_es(L_N, q))
        c = colours(idx)
        label = f"Cl {cl_id} ({n_cells} cells)"
        ax1.plot(p_values, vars_, "-o", ms=2, lw=1.2, color=c, label=label)
        ax2.plot(p_values, ess_, "-s", ms=2, lw=1.2, color=c, label=label)

    for ax, title in [(ax1, "VaR$_{0.95}$"), (ax2, "ES$_{0.95}$")]:
        ax.set_xlabel("Threshold probability $p$", fontsize=11)
        ax.set_ylabel(title, fontsize=11)
        ax.set_title(f"{title} vs $p$ — top {top_k} LEC clusters", fontsize=12)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"threshold_sensitivity.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  ✔ threshold_sensitivity")


# ═════════════════════════════════════════════════════════════════════
#  PLOT 5a — Tail shape: sorted L_N with VaR/ES annotated
# ═════════════════════════════════════════════════════════════════════
def plot_tail_shape(results, cname, top_k=8):
    """For each cluster, plot the sorted L_N values (empirical CDF / order
    statistics) with VaR and ES clearly marked.  This shows the full
    distribution shape and how the tail beyond VaR generates ES."""
    res = results[cname]
    sorted_cl = sorted(res.clusters, key=lambda c: c["ES"], reverse=True)[:top_k]

    fig, axes = plt.subplots(
        (top_k + 1) // 2, 2,
        figsize=(14, 2.8 * ((top_k + 1) // 2)),
    )
    axes = axes.flatten()

    for i, cl in enumerate(sorted_cl):
        ax = axes[i]
        L = cl["L_N"]
        n = len(L)
        srt = np.sort(L)
        q = res.q
        k = int(np.ceil(n * q))  # 1-indexed VaR position

        # Empirical quantile positions (probability axis)
        probs = np.arange(1, n + 1) / n

        # Plot all order statistics
        ax.step(probs, srt, where="post", lw=0.9, color="steelblue",
                label="Sorted $L_N$")

        # Shade the tail beyond VaR
        tail_mask = np.arange(n) >= (k - 1)
        ax.fill_between(probs, 0, srt, where=tail_mask, step="post",
                        color="firebrick", alpha=0.15,
                        label=f"Tail (worst {n - k + 1} yrs)")

        # VaR line
        ax.axhline(cl["VaR"], ls="--", lw=1.0, color="firebrick", alpha=0.7)
        ax.annotate(
            f"VaR$_{{0.95}}$={cl['VaR']:.3f}",
            xy=(q, cl["VaR"]), xytext=(q - 0.18, cl["VaR"] + 0.03),
            fontsize=7, color="firebrick",
            arrowprops=dict(arrowstyle="->", color="firebrick", lw=0.6),
        )

        # ES line
        ax.axhline(cl["ES"], ls="-.", lw=1.0, color="darkred", alpha=0.7)
        ax.annotate(
            f"ES$_{{0.95}}$={cl['ES']:.3f}",
            xy=(0.98, cl["ES"]), xytext=(0.75, cl["ES"] + 0.04),
            fontsize=7, color="darkred",
            arrowprops=dict(arrowstyle="->", color="darkred", lw=0.6),
        )

        # VaR vertical marker
        ax.axvline(q, ls=":", lw=0.6, color="grey", alpha=0.5)

        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, max(srt[-1] * 1.15, 0.05))
        ax.set_xlabel("Empirical probability", fontsize=8)
        ax.set_ylabel("$L_N$ (area fraction)", fontsize=8)
        ax.set_title(
            f"Cluster {cl['cluster']} ({cl['n_cells']} cells)",
            fontsize=9,
        )
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper left")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"Tail Shape of $L_N$ — {cname} (top {top_k} by ES, "
        f"$p$={res.p}, $q$={res.q})",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    cost_tag = getattr(res, 'cost', 'indicator')
    suffix = f"_{cost_tag}" if cost_tag != "indicator" else ""
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"tail_shape_{cname}{suffix}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ tail_shape_{cname}{suffix}")


# ═════════════════════════════════════════════════════════════════════
#  PLOT 5a′ — Focused tail comparison: 3 clusters side by side
# ═════════════════════════════════════════════════════════════════════
def plot_tail_comparison(results, cname, focus_ids=(5, 6, 10)):
    """Paper-quality 3-panel figure comparing the tail shapes of three
    carefully chosen clusters that illustrate different risk regimes:
      - Heavy tail (small region, strong co-occurrence)
      - Strong dependence at scale (large region, still risky)
      - Spatial diversification (very large, diluted risk)
    """
    res = results[cname]
    cl_map = {c["cluster"]: c for c in res.clusters}
    focus = [cl_map[fid] for fid in focus_ids]

    subtitles = {
        5:  "Heavy tail\n(95 cells, ~14°S 171°E)",
        6:  "Strong dependence at scale\n(2770 cells, ~13°S 197°E)",
        10: "Spatial diversification\n(5028 cells, ~5°N 172°E)",
    }

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharey=False)
    colours = ["#D62728", "#2CA02C", "#1F77B4"]  # red, green, blue

    for ax, cl, col in zip(axes, focus, colours):
        L = cl["L_N"]
        n = len(L)
        srt = np.sort(L)
        q = res.q
        k = int(np.ceil(n * q))
        probs = np.arange(1, n + 1) / n

        # Full distribution
        ax.step(probs, srt, where="post", lw=1.2, color=col, alpha=0.9,
                label="Sorted $L_N$")

        # Shade tail
        tail_mask = np.arange(n) >= (k - 1)
        ax.fill_between(probs, 0, srt, where=tail_mask, step="post",
                        color=col, alpha=0.12)

        # VaR line + label
        ax.axhline(cl["VaR"], ls="--", lw=1.0, color="0.3", alpha=0.7)
        # Position VaR label to the left of the threshold line
        txt_y_var = cl["VaR"] + max(srt[-1] * 0.04, 0.005)
        ax.text(0.50, txt_y_var, f"VaR$_{{0.95}}$ = {cl['VaR']:.3f}",
                fontsize=9, color="0.2", ha="center")

        # ES line + label
        ax.axhline(cl["ES"], ls="-.", lw=1.0, color="0.3", alpha=0.7)
        txt_y_es = cl["ES"] + max(srt[-1] * 0.04, 0.005)
        ax.text(0.50, txt_y_es, f"ES$_{{0.95}}$ = {cl['ES']:.3f}",
                fontsize=9, color="0.2", ha="center",
                fontweight="bold")

        # Shade gap between VaR and ES to emphasise tail heaviness
        ax.axhspan(cl["VaR"], cl["ES"], color=col, alpha=0.06)

        # q vertical line
        ax.axvline(q, ls=":", lw=0.6, color="grey", alpha=0.5)
        ax.text(q + 0.005, srt[-1] * 0.02, "$q$=0.95",
                fontsize=7, color="grey", rotation=90, va="bottom")

        # Ratio annotation
        if cl["VaR"] > 0:
            ratio = cl["ES"] / cl["VaR"]
            ax.text(0.03, srt[-1] * 0.92,
                    f"ES/VaR = {ratio:.1f}×",
                    fontsize=10, color=col, fontweight="bold",
                    transform=ax.get_yaxis_transform())

        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, max(srt[-1] * 1.18, 0.05))
        ax.set_xlabel("Empirical probability", fontsize=10)
        ax.set_ylabel("$L_N$ (fraction of area in exceedance)", fontsize=10)
        ax.set_title(
            f"Cluster {cl['cluster']}  —  "
            + subtitles.get(cl["cluster"], ""),
            fontsize=10, pad=8,
        )
        ax.tick_params(labelsize=8)

    fig.suptitle(
        f"Three Risk Regimes — Tail Shape of $L_N$  "
        f"($p$={res.p}, $q$={res.q})",
        fontsize=13, y=1.03,
    )
    fig.tight_layout()
    cost_tag = getattr(res, 'cost', 'indicator')
    suffix = f"_{cost_tag}" if cost_tag != "indicator" else ""
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"tail_comparison_{cname}{suffix}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ tail_comparison_{cname}{suffix}")


# ═════════════════════════════════════════════════════════════════════
#  PLOT 5b — Bar chart: VaR & ES per cluster
# ═════════════════════════════════════════════════════════════════════
def plot_cluster_bar_chart(results, cname):
    """Horizontal bar chart showing VaR and ES for every cluster,
    sorted by ES descending.  Each bar is labelled with cluster ID
    and cell count so the reader can cross-reference Figure 9."""
    res = results[cname]
    clusters = sorted(res.clusters, key=lambda c: c["ES"], reverse=True)

    ids    = [c["cluster"] for c in clusters]
    n_c    = [c["n_cells"] for c in clusters]
    vars_  = [c["VaR"] for c in clusters]
    ess_   = [c["ES"]  for c in clusters]
    labels = [f"Cl {i}  ({nc} cells)" for i, nc in zip(ids, n_c)]

    y = np.arange(len(clusters))
    h = 0.35

    fig, ax = plt.subplots(figsize=(9, max(5, 0.45 * len(clusters))))
    bars_es  = ax.barh(y - h/2, ess_,  h, label=f"ES$_{{0.95}}$",
                        color="firebrick", alpha=0.8)
    bars_var = ax.barh(y + h/2, vars_, h, label=f"VaR$_{{0.95}}$",
                        color="steelblue", alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()  # highest ES on top
    ax.set_xlabel("Risk measure value (fraction of cluster area in exceedance)",
                  fontsize=10)
    ax.set_title(
        f"Koch Risk per {cname} Cluster  ($p$={res.p}, $q$={res.q})",
        fontsize=12,
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(0, min(1.05, max(ess_) * 1.15))
    ax.grid(axis="x", alpha=0.3)

    fig.tight_layout()
    cost_tag = getattr(res, 'cost', 'indicator')
    suffix = f"_{cost_tag}" if cost_tag != "indicator" else ""
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"cluster_bar_{cname}{suffix}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ cluster_bar_{cname}{suffix}")


# ═════════════════════════════════════════════════════════════════════
#  PLOT 6 — 3-panel summary: Cluster map | VaR map | ES map
# ═════════════════════════════════════════════════════════════════════
def plot_summary_panel(results, cname, valid_idx, labels, lats_1d, lons_1d):
    """Side-by-side panel linking Figure 9 clusters to VaR/ES maps.
    Panel (a): cluster identities (categorical colours, matching Fig 9).
    Panel (b): same cells coloured by cluster VaR.
    Panel (c): same cells coloured by cluster ES.
    The reader sees: cluster → VaR number → ES number in one figure."""
    from matplotlib.colors import ListedColormap, BoundaryNorm

    res = results[cname]
    n_lat, n_lon = len(lats_1d), len(lons_1d)
    K = len(res.clusters)

    # ── Build metric arrays for each cell ──
    cl_var = {}
    cl_es  = {}
    for c in res.clusters:
        cl_var[c["cluster"]] = c["VaR"]
        cl_es[c["cluster"]]  = c["ES"]

    cell_ids  = np.array([labels[i] for i in range(len(labels))], dtype=float)
    cell_var  = np.array([cl_var[labels[i]] for i in range(len(labels))])
    cell_es   = np.array([cl_es[labels[i]]  for i in range(len(labels))])

    grid_id  = _build_grid_field(cell_ids,  valid_idx, n_lat, n_lon)
    grid_var = _build_grid_field(cell_var,   valid_idx, n_lat, n_lon)
    grid_es  = _build_grid_field(cell_es,    valid_idx, n_lat, n_lon)

    # ── Categorical colourmap for cluster IDs ──
    unique_ids = sorted(np.unique(labels).tolist())
    base_cols = (list(plt.get_cmap("tab20").colors)
                 + list(plt.get_cmap("Set3").colors))
    id_colours = [base_cols[i % len(base_cols)] for i in range(len(unique_ids))]
    id_cmap = ListedColormap(id_colours)
    bounds = [uid - 0.5 for uid in unique_ids] + [unique_ids[-1] + 0.5]
    id_norm = BoundaryNorm(bounds, id_cmap.N)

    # ── Common VaR/ES colour range ──
    v = np.concatenate([cell_var, cell_es])
    vmin, vmax = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    if vmin == vmax:
        vmin, vmax = 0, 1

    fig = plt.figure(figsize=(22, 6))

    # Panel (a): cluster map
    ax1 = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson())
    ax1.set_global()
    ax1.coastlines(linewidth=0.4, color="0.3")
    ax1.gridlines(draw_labels=False, linewidth=0.15, color="grey", alpha=0.5)
    im1 = ax1.pcolormesh(
        lons_1d, lats_1d, grid_id,
        transform=ccrs.PlateCarree(),
        cmap=id_cmap, norm=id_norm,
        shading="auto", rasterized=True,
    )
    cb1 = fig.colorbar(im1, ax=ax1, orientation="horizontal",
                       pad=0.04, shrink=0.7, aspect=25,
                       ticks=unique_ids[::max(1, K // 10)])
    cb1.set_label("Cluster ID", fontsize=10)
    ax1.set_title(f"(a) {cname} Clusters (Figure 9)", fontsize=11, pad=10)

    # Panel (b): VaR map
    ax2 = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson())
    ax2.set_global()
    ax2.coastlines(linewidth=0.4, color="0.3")
    ax2.gridlines(draw_labels=False, linewidth=0.15, color="grey", alpha=0.5)
    im2 = ax2.pcolormesh(
        lons_1d, lats_1d, grid_var,
        transform=ccrs.PlateCarree(),
        cmap="YlOrRd", vmin=vmin, vmax=vmax,
        shading="auto", rasterized=True,
    )
    cb2 = fig.colorbar(im2, ax=ax2, orientation="horizontal",
                       pad=0.04, shrink=0.7, aspect=25)
    cb2.set_label(f"VaR$_{{0.95}}$", fontsize=10)
    ax2.set_title("(b) VaR$_{0.95}$ per cluster", fontsize=11, pad=10)

    # Panel (c): ES map
    ax3 = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson())
    ax3.set_global()
    ax3.coastlines(linewidth=0.4, color="0.3")
    ax3.gridlines(draw_labels=False, linewidth=0.15, color="grey", alpha=0.5)
    im3 = ax3.pcolormesh(
        lons_1d, lats_1d, grid_es,
        transform=ccrs.PlateCarree(),
        cmap="YlOrRd", vmin=vmin, vmax=vmax,
        shading="auto", rasterized=True,
    )
    cb3 = fig.colorbar(im3, ax=ax3, orientation="horizontal",
                       pad=0.04, shrink=0.7, aspect=25)
    cb3.set_label(f"ES$_{{0.95}}$", fontsize=10)
    ax3.set_title("(c) ES$_{0.95}$ per cluster", fontsize=11, pad=10)

    fig.suptitle(
        f"From {cname} Clusters to Spatial Risk  "
        f"($p$={res.p}, $q$={res.q})",
        fontsize=13, y=1.02,
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"summary_panel_{cname}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ summary_panel_{cname}")


def plot_cluster_lookup_panel(results, cname, valid_idx, labels, lats_1d, lons_1d):
    """Figure-9-like cluster map + explicit VaR/ES lookup table.

    This avoids the common confusion that VaR/ES are per grid cell:
    the left panel keeps the categorical cluster colours, while the
    right panel lists exactly one VaR and one ES per cluster ID.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm

    res = results[cname]
    n_lat, n_lon = len(lats_1d), len(lons_1d)

    # Build cluster ID grid
    cell_ids = np.array([labels[i] for i in range(len(labels))], dtype=float)
    grid_id = _build_grid_field(cell_ids, valid_idx, n_lat, n_lon)

    # Categorical colour map matching cluster IDs
    unique_ids = sorted(np.unique(labels).tolist())
    base_cols = (list(plt.get_cmap("tab20").colors)
                 + list(plt.get_cmap("Set3").colors))
    id_colours = [base_cols[i % len(base_cols)] for i in range(len(unique_ids))]
    id_cmap = ListedColormap(id_colours)
    bounds = [uid - 0.5 for uid in unique_ids] + [unique_ids[-1] + 0.5]
    id_norm = BoundaryNorm(bounds, id_cmap.N)

    fig = plt.figure(figsize=(18, 7))

    # Left: Figure-9-like cluster map
    ax_map = fig.add_subplot(1, 2, 1, projection=ccrs.Robinson())
    ax_map.set_global()
    ax_map.coastlines(linewidth=0.4, color="0.3")
    ax_map.gridlines(draw_labels=False, linewidth=0.15, color="grey", alpha=0.5)
    im = ax_map.pcolormesh(
        lons_1d, lats_1d, grid_id,
        transform=ccrs.PlateCarree(),
        cmap=id_cmap, norm=id_norm,
        shading="auto", rasterized=True,
    )
    cb = fig.colorbar(
        im, ax=ax_map, orientation="horizontal",
        pad=0.04, shrink=0.8, aspect=28,
        ticks=unique_ids,
    )
    cb.set_label("Cluster ID", fontsize=10)
    ax_map.set_title(f"(a) {cname} clusters (Figure 9 style)", fontsize=11, pad=10)

    # Overlay cluster IDs at approximate centroids
    for cid in unique_ids:
        mask = labels == cid
        idx = valid_idx[mask]
        lat_idx = (idx // n_lon).astype(int)
        lon_idx = (idx % n_lon).astype(int)
        lat_c = float(np.mean(lats_1d[lat_idx]))
        lon_c = float(np.mean(lons_1d[lon_idx]))
        ax_map.text(
            lon_c, lat_c, str(int(cid)), transform=ccrs.PlateCarree(),
            ha="center", va="center", fontsize=6.5, color="black",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="0.3", alpha=0.75),
        )

    # Right: explicit lookup table (one VaR/ES per cluster)
    ax_tbl = fig.add_subplot(1, 2, 2)
    ax_tbl.axis("off")
    ax_tbl.set_title(
        "(b) Cluster-level lookup (one VaR and one ES per cluster)",
        fontsize=11, pad=10,
    )

    rows = []
    for c in sorted(res.clusters, key=lambda x: x["cluster"]):
        rows.append([
            int(c["cluster"]),
            int(c["n_cells"]),
            f"{c['VaR']:.4f}",
            f"{c['ES']:.4f}",
        ])

    table = ax_tbl.table(
        cellText=rows,
        colLabels=["Cluster", "Cells", "VaR$_{0.95}$", "ES$_{0.95}$"],
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.2)

    # Header style
    for j in range(4):
        table[(0, j)].set_facecolor("#efefef")
        table[(0, j)].set_text_props(weight="bold")

    fig.suptitle(
        f"{cname} clusters with explicit per-cluster risk values  "
        f"($p$={res.p}, $q$={res.q})",
        fontsize=13, y=0.98,
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"cluster_lookup_{cname}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ cluster_lookup_{cname}")


# ═════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════
def main():
    print("Computing Koch risk measures (E=1, p=0.95, q=0.95) ...")
    results = compute_koch_risk(DATA_PATH, p=0.95, q=0.95)

    d = np.load(DATA_PATH)
    valid_idx = d["valid_idx"]
    lats_1d = d["lats"]
    lons_1d = d["lons"]

    label_map = {"LEC": d["labels_lec"], "EDC": d["labels_edc"]}

    for cname in ["LEC", "EDC"]:
        k = len(results[cname].clusters)
        top = sorted(results[cname].clusters, key=lambda c: c["ES"], reverse=True)[0]
        print(f"\n  {cname}: {k} clusters, "
              f"highest ES={top['ES']:.4f} (cluster {top['cluster']}, {top['n_cells']} cells)")

    print("\n── Temporal L_N plots ──")
    for cname in ["LEC", "EDC"]:
        plot_temporal_loss(results, cname)

    print("\n── VaR world maps ──")
    for cname in ["LEC", "EDC"]:
        plot_var_map(cname, results, valid_idx, label_map[cname], lats_1d, lons_1d)

    print("\n── ES world maps ──")
    for cname in ["LEC", "EDC"]:
        plot_es_map(cname, results, valid_idx, label_map[cname], lats_1d, lons_1d)

    print("\n── Threshold sensitivity ──")
    plot_threshold_sensitivity()

    print("\n── Tail shape plots ──")
    for cname in ["LEC", "EDC"]:
        plot_tail_shape(results, cname)

    print("\n── Focused tail comparison (LEC clusters 5, 6, 10) ──")
    plot_tail_comparison(results, "LEC", focus_ids=(5, 6, 10))

    print("\n── Cluster bar charts ──")
    for cname in ["LEC", "EDC"]:
        plot_cluster_bar_chart(results, cname)

    print("\n── Summary panels (clusters → VaR → ES) ──")
    for cname in ["LEC", "EDC"]:
        plot_summary_panel(results, cname, valid_idx, label_map[cname], lats_1d, lons_1d)

    print("\n── Figure-9-style cluster lookup panels ──")
    for cname in ["LEC", "EDC"]:
        plot_cluster_lookup_panel(results, cname, valid_idx, label_map[cname], lats_1d, lons_1d)

    # ── Extension 1: Excess-based loss ──────────────────────────────
    print("\n══ Excess-based cost E(z) = (z − u_p)⁺ ══")
    print("Computing Koch risk measures (cost=excess) ...")
    results_ex = compute_koch_risk(DATA_PATH, p=0.95, q=0.95, cost="excess")
    for cname in ["LEC", "EDC"]:
        k = len(results_ex[cname].clusters)
        top = sorted(results_ex[cname].clusters, key=lambda c: c["ES"], reverse=True)[0]
        print(f"  {cname}: highest ES={top['ES']:.4f} "
              f"(cluster {top['cluster']}, {top['n_cells']} cells)")

    print("\n── Excess: tail comparison (LEC clusters 5, 6, 10) ──")
    plot_tail_comparison(results_ex, "LEC", focus_ids=(5, 6, 10))

    print("\n── Excess: cluster bar charts ──")
    for cname in ["LEC", "EDC"]:
        plot_cluster_bar_chart(results_ex, cname)

    print("\n── Excess: tail shape (top 8) ──")
    for cname in ["LEC", "EDC"]:
        plot_tail_shape(results_ex, cname)

    print("\n── Indicator vs Excess comparison panel ──")
    plot_indicator_vs_excess(results, results_ex)

    # ── Extension 2: Temporal trend ─────────────────────────────────
    print("\n══ Temporal trend analysis ══")
    print("\n── Trend: indicator-based ──")
    plot_temporal_trend(results, cost_label="indicator")
    print("\n── Trend: excess-based ──")
    plot_temporal_trend(results_ex, cost_label="excess")

    # ── Extension 3: Cross-cluster co-exceedance ───────────────────
    print("\n══ Cross-cluster co-exceedance analysis ══")
    for cost_label, res_dict in [("indicator", results), ("excess", results_ex)]:
        cross = compute_cross_cluster_dependence(res_dict["LEC"])
        print(f"\n── Co-exceedance heatmap ({cost_label}) ──")
        plot_coexceedance_heatmap(cross, cost_label=cost_label)
        print(f"\n── Top co-exceedance pairs ({cost_label}) ──")
        plot_top_coexceedance_pairs(cross, res_dict, cost_label=cost_label)

        # Print summary of top-5 pairs
        K = len(cross.cluster_ids)
        pairs = []
        for i in range(K):
            for j in range(i + 1, K):
                pairs.append((
                    cross.cluster_ids[i],
                    cross.cluster_ids[j],
                    cross.coexceedance_ratio[i, j],
                    cross.n_joint_years[i, j],
                ))
        pairs.sort(key=lambda p: p[2], reverse=True)
        print(f"\n  Top-5 systemic pairs ({cost_label}):")
        for a, b, r, n in pairs[:5]:
            print(f"    Cl {a:2d}–{b:2d}: ratio={r:.1f}× ({n} joint years)")

    print(f"\nAll outputs saved to {OUT_DIR}/")


# ═════════════════════════════════════════════════════════════════════
#  EXTENSION 1 — Indicator vs Excess side-by-side comparison
# ═════════════════════════════════════════════════════════════════════
def plot_indicator_vs_excess(results_ind, results_exc, cname="LEC", focus_ids=(5, 6, 10)):
    """2×3 panel: top row = indicator cost, bottom row = excess cost,
    columns = three focus clusters. Shows how tail shape changes when
    moving from binary to severity-weighted loss.
    """
    res_ind = results_ind[cname]
    res_exc = results_exc[cname]
    cl_ind = {c["cluster"]: c for c in res_ind.clusters}
    cl_exc = {c["cluster"]: c for c in res_exc.clusters}

    fig, axes = plt.subplots(2, 3, figsize=(16, 7))

    titles = {
        5:  "Cl 5 — Heavy tail",
        6:  "Cl 6 — Dependence at scale",
        10: "Cl 10 — Diversification",
    }
    row_labels = [
        r"Indicator: $E(z) = \mathbf{1}\{z > u_p\}$",
        r"Excess: $E(z) = (z - u_p)^+$",
    ]
    colours = ["#D62728", "#2CA02C", "#1F77B4"]

    for row, (res, cl_map, label) in enumerate(
        [(res_ind, cl_ind, row_labels[0]), (res_exc, cl_exc, row_labels[1])]
    ):
        q = res.q
        for col, (fid, colour) in enumerate(zip(focus_ids, colours)):
            ax = axes[row, col]
            cl = cl_map[fid]
            L = cl["L_N"]
            n = len(L)
            srt = np.sort(L)
            k = int(np.ceil(n * q))
            probs = np.arange(1, n + 1) / n

            ax.step(probs, srt, where="post", lw=1.2, color=colour, alpha=0.9)
            tail_mask = np.arange(n) >= (k - 1)
            ax.fill_between(probs, 0, srt, where=tail_mask, step="post",
                            color=colour, alpha=0.12)

            ax.axhline(cl["VaR"], ls="--", lw=0.8, color="0.3", alpha=0.6)
            ax.axhline(cl["ES"], ls="-.", lw=0.8, color="0.3", alpha=0.6)

            # Labels
            y_top = max(srt[-1] * 1.15, 0.01)
            ax.text(0.50, cl["VaR"] + y_top * 0.03,
                    f"VaR={cl['VaR']:.3f}", fontsize=8, ha="center", color="0.3")
            ax.text(0.50, cl["ES"] + y_top * 0.03,
                    f"ES={cl['ES']:.3f}", fontsize=8, ha="center", color="0.3",
                    fontweight="bold")

            if cl["VaR"] > 0:
                ratio = cl["ES"] / cl["VaR"]
                ax.text(0.03, srt[-1] * 0.88, f"ES/VaR={ratio:.1f}×",
                        fontsize=9, color=colour, fontweight="bold",
                        transform=ax.get_yaxis_transform())

            ax.set_xlim(0, 1.02)
            ax.set_ylim(0, y_top)
            if col == 0:
                ax.set_ylabel(label, fontsize=9.5)
            if row == 0:
                ax.set_title(titles[fid], fontsize=10, pad=6)
            if row == 1:
                ax.set_xlabel("Empirical probability", fontsize=9)
            ax.tick_params(labelsize=7)

    fig.suptitle(
        f"Indicator vs Excess Cost — Tail Shape of $L_N$  "
        f"({cname}, $p$={res_ind.p}, $q$={res_ind.q})",
        fontsize=13, y=1.02,
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"indicator_vs_excess_{cname}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ indicator_vs_excess_{cname}")


# ═════════════════════════════════════════════════════════════════════
#  EXTENSION 2 — Temporal trend in L_N
# ═════════════════════════════════════════════════════════════════════
def plot_temporal_trend(results, cname="LEC", top_k=8, cost_label="indicator"):
    """Time series of L_N with OLS trend line + Mann-Kendall significance.

    For each of the top-k clusters (by ES), shows:
      - annual L_N in grey
      - OLS trend line in colour
      - slope per decade annotated
      - Mann-Kendall p-value: ** (p<0.01), * (p<0.05), n.s.
    """
    res = results[cname]
    sorted_cl = sorted(res.clusters, key=lambda c: c["ES"], reverse=True)[:top_k]

    n_cols = 2
    n_rows = (top_k + 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 2.8 * n_rows), sharex=True)
    axes = axes.flatten()
    cmap = matplotlib.colormaps.get_cmap("tab10").resampled(top_k)

    for i, cl in enumerate(sorted_cl):
        ax = axes[i]
        years = res.years
        L = cl["L_N"]
        colour = cmap(i)

        trend = compute_trend(years, L)
        t = years - years[0]
        fit_line = trend["slope"] * t + trend["intercept"]

        # Data
        ax.plot(years, L, lw=0.5, color="0.6", alpha=0.7)
        # Trend
        ax.plot(years, fit_line, lw=1.8, color=colour, alpha=0.9)

        # Significance marker
        mk_p = trend["mk_pvalue"]
        if mk_p < 0.01:
            sig = "**"
        elif mk_p < 0.05:
            sig = "*"
        else:
            sig = "n.s."

        slope_dec = trend["slope_per_decade"]
        sign = "+" if slope_dec >= 0 else ""

        # Annotation
        ax.text(
            0.02, 0.95,
            f"Cl {cl['cluster']}  ({cl['n_cells']} cells)\n"
            f"slope = {sign}{slope_dec:.4f}/decade  {sig}\n"
            f"(MK p = {mk_p:.3f})",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.85),
        )

        ax.set_ylabel("$L_N$", fontsize=9)
        ax.tick_params(labelsize=7)

    # Hide unused panels
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    axes[-2].set_xlabel("Year", fontsize=10)
    axes[-1].set_xlabel("Year", fontsize=10)

    cost_str = "indicator" if cost_label == "indicator" else "excess"
    fig.suptitle(
        f"Temporal Trend in $L_N$ — {cname} top-{top_k} clusters\n"
        f"(cost = {cost_str}, $p$={res.p}, $q$={res.q})",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    tag = f"trend_{cost_label}_{cname}"
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{tag}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {tag}")


# ═════════════════════════════════════════════════════════════════════
#  EXTENSION 3 — Cross-cluster co-exceedance heatmap
# ═════════════════════════════════════════════════════════════════════
def plot_coexceedance_heatmap(cross, cost_label="indicator"):
    """Heatmap of co-exceedance ratio for all cluster pairs.

    Values > 1 (red) mean joint extreme years occur more often than
    independence predicts; ≈ 1 (white) means independent; < 1 (blue)
    would mean anti-correlated extremes.
    """
    ratio = cross.coexceedance_ratio.copy()
    K = ratio.shape[0]
    ids = cross.cluster_ids

    # Mask diagonal (self-comparison is always 1/P_marginal, not useful)
    np.fill_diagonal(ratio, np.nan)

    fig, ax = plt.subplots(figsize=(10, 8.5))
    cmap = matplotlib.colormaps.get_cmap("RdBu_r").copy()
    cmap.set_bad("0.85")  # grey for diagonal

    # Clip for colour scale — ratios can be very large for rare pairs
    vmax = min(np.nanmax(ratio), 20.0)
    im = ax.imshow(ratio, cmap=cmap, vmin=0, vmax=vmax,
                   interpolation="nearest", aspect="equal")

    ax.set_xticks(range(K))
    ax.set_yticks(range(K))
    ax.set_xticklabels(ids, fontsize=7, rotation=90)
    ax.set_yticklabels(ids, fontsize=7)
    ax.set_xlabel("Cluster ID", fontsize=10)
    ax.set_ylabel("Cluster ID", fontsize=10)

    cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label("Co-exceedance ratio  "
                 r"$\hat{P}_{\mathrm{joint}} \,/\, P_{\mathrm{indep}}$",
                 fontsize=10)

    ax.set_title(
        f"Cross-Cluster Co-exceedance — {cross.clustering_name}\n"
        f"(cost={cost_label}, $q$={cross.q})    "
        r"ratio $> 1$ $\Rightarrow$ systemic risk",
        fontsize=12, pad=10,
    )

    # Annotate cells with ratio values for small K
    if K <= 25:
        for i in range(K):
            for j in range(K):
                if i == j:
                    continue
                val = cross.coexceedance_ratio[i, j]
                n = cross.n_joint_years[i, j]
                colour = "white" if val > vmax * 0.65 else "0.2"
                ax.text(j, i, f"{val:.1f}\n({n})",
                        ha="center", va="center", fontsize=5.5,
                        color=colour)

    fig.tight_layout()
    tag = f"coexceedance_{cross.clustering_name}_{cost_label}"
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{tag}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {tag}")


def plot_top_coexceedance_pairs(cross, results, cost_label="indicator", top_n=15):
    """Bar chart of the top-N cluster pairs by co-exceedance ratio,
    with annotations showing joint years count and cluster sizes.
    """
    K = len(cross.cluster_ids)
    cl_map = {c["cluster"]: c for c in results[cross.clustering_name].clusters}

    # Collect off-diagonal pairs
    pairs = []
    for i in range(K):
        for j in range(i + 1, K):
            cid_a = cross.cluster_ids[i]
            cid_b = cross.cluster_ids[j]
            pairs.append({
                "pair": f"Cl {cid_a}–{cid_b}",
                "ratio": cross.coexceedance_ratio[i, j],
                "n_joint": cross.n_joint_years[i, j],
                "p_joint": cross.joint_prob[i, j],
                "p_indep": cross.indep_prob[i, j],
                "cells_a": cl_map[cid_a]["n_cells"],
                "cells_b": cl_map[cid_b]["n_cells"],
            })

    pairs.sort(key=lambda p: p["ratio"], reverse=True)
    top = pairs[:top_n]

    fig, ax = plt.subplots(figsize=(10, max(5, 0.45 * top_n)))
    y = np.arange(len(top))
    ratios = [p["ratio"] for p in top]
    labels = [p["pair"] for p in top]

    colours = ["#D62728" if r > 3 else "#FF7F0E" if r > 1.5 else "#1F77B4"
               for r in ratios]
    bars = ax.barh(y, ratios, color=colours, alpha=0.85, edgecolor="0.3", lw=0.3)

    # Independence line
    ax.axvline(1.0, ls="--", lw=1.0, color="0.4", alpha=0.7)
    ax.text(1.05, -0.8, "independence", fontsize=8, color="0.4", va="top")

    # Annotations
    for i, p in enumerate(top):
        ax.text(
            p["ratio"] + 0.15, i,
            f"{p['n_joint']} joint yrs  "
            f"({p['cells_a']}+{p['cells_b']} cells)  "
            f"$\\hat{{P}}$={p['p_joint']:.3f} vs {p['p_indep']:.4f}",
            va="center", fontsize=7, color="0.3",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Co-exceedance ratio  "
                  r"($\hat{P}_{\mathrm{joint}} / P_{\mathrm{indep}}$)",
                  fontsize=10)
    ax.set_title(
        f"Top {top_n} Cluster Pairs by Systemic Co-exceedance — "
        f"{cross.clustering_name}\n"
        f"(cost={cost_label}, $q$={cross.q})",
        fontsize=12, pad=8,
    )
    ax.set_xlim(0, max(ratios) * 1.35)
    ax.grid(axis="x", alpha=0.2)

    fig.tight_layout()
    tag = f"top_coexceedance_{cross.clustering_name}_{cost_label}"
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{tag}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ {tag}")


if __name__ == "__main__":
    main()

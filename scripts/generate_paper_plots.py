#!/usr/bin/env python3
"""Generate all figures for the SoftwareX paper.

Produces plots similar to the Justus (2024) reference paper:
- Fig 1: Extremal coefficient maps (4 reference points, stripes preset)
- Fig 3: Stripes preset — EDC clusters (b), LEC clusters (b), true b
- Fig 6: Rotate preset — EDC clusters (gamma), LEC clusters (gamma), true gamma
- Fig 7: Bigsmall preset — EDC clusters (a), LEC clusters (a), true a
- Fig 4/5: EICI better/worse comparison heatmap

Usage:
    python scripts/generate_paper_plots.py
"""

from __future__ import annotations

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import fcluster

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.grid import Grid
from weatherisk.covariance import cov_fkt_2d, cov_fkt_2d_nonstat2, cov_to_ec
from weatherisk.simulation import sim_expt_2d_nonstat
from weatherisk.density import pairwise_density_optim_local
from weatherisk.estimation import smooth_local_estimates, llh_in_cluster
from weatherisk.clustering import (
    calc_distance_ellipses,
    c_extrcoeff_matrix,
    clustering,
)
from weatherisk.parameters import get_preset
from weatherisk.pipeline import run_nonstationary_pipeline
from weatherisk.plotting import plot_cluster_comparison

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Colour helpers ──────────────────────────────────────────────────

def blues_cmap(n=101):
    """Reversed 'Blues' palette similar to R's brewer.pal(9,'Blues')."""
    return plt.cm.Blues_r


def rdbu_cmap():
    """RdBu palette for parameter maps."""
    return plt.cm.RdBu_r


def _make_cluster_cmap(k):
    """Qualitative colormap with k distinct colours."""
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    return ListedColormap([base[i % len(base)] for i in range(max(k, 1))])


# ── Plot functions matching R paper style ───────────────────────────

def plot_parameter_map(grid, values_2d, vmin, vmax, title="", cmap=None,
                       cluster_contours=None, filename=None, label=""):
    """Plot a heatmap of parameter values on a grid with optional cluster contours.

    Similar to R's plot_map() with contours overlay.
    """
    if cmap is None:
        cmap = blues_cmap()

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.pcolormesh(
        grid.x_ax, grid.y_ax, values_2d,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    plt.colorbar(im, ax=ax, label=label, shrink=0.8)

    # Draw cluster contours if provided
    if cluster_contours is not None:
        cl_2d = cluster_contours.reshape(grid.nrow, grid.ncol)
        ax.contour(
            grid.x_ax, grid.y_ax, cl_2d,
            colors="black", linewidths=0.8, alpha=0.6,
        )

    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()

    if filename:
        fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300, bbox_inches="tight")
        print(f"  Saved {filename}")
    plt.close(fig)
    return fig


def plot_ec_map(grid, ec_values, title="", filename=None):
    """Plot extremal coefficient map for a single reference point."""
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.pcolormesh(
        grid.x_ax, grid.y_ax, ec_values,
        shading="auto", cmap="RdBu_r", vmin=1.0, vmax=2.0,
    )
    plt.colorbar(im, ax=ax, label=r"$\theta$", shrink=0.8)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()
    if filename:
        fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300, bbox_inches="tight")
        print(f"  Saved {filename}")
    plt.close(fig)
    return fig


def plot_eici_map(grid, eici_2d, title="", filename=None):
    """Plot EICI better/worse heatmap.

    0 = LEC better, 1 = EDC better. Mean across runs gives [0,1].
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    cmap = plt.cm.RdBu  # Blue = LEC better, Red = EDC better
    im = ax.pcolormesh(
        grid.x_ax, grid.y_ax, eici_2d * 100,
        shading="auto", cmap=cmap, vmin=0, vmax=100,
    )
    cb = plt.colorbar(im, ax=ax, label="% runs where EDC better", shrink=0.8)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()
    if filename:
        fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300, bbox_inches="tight")
        print(f"  Saved {filename}")
    plt.close(fig)
    return fig


# ── Full pipeline for one preset ────────────────────────────────────

def run_full_preset(preset_name, resolution=None, n_sim=None, seed=42):
    """Run the full pipeline via weatherisk.pipeline.run_nonstationary_pipeline.

    This is a thin wrapper that delegates to the library function.
    """
    return run_nonstationary_pipeline(
        preset=preset_name,
        resolution=resolution,
        n_sim=n_sim,
        seed=seed,
        threshold_quantile=0.30,
        verbose=True,
    )


# ── Figure generators ──────────────────────────────────────────────

def gen_fig_ec_maps(result, preset_name):
    """Generate extremal coefficient maps for 4 reference points (cf. Fig. 1)."""
    grid = result["grid"]
    p = result["preset"]
    a_matrix = result["a_matrix"]
    b_matrix = result["b_matrix"]
    g_matrix = result["g_matrix"]

    ref_points = [(-3, 2), (-3, -2), (3, 2), (3, -2)]

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    a_flat = a_matrix.ravel()
    b_flat = b_matrix.ravel()
    g_flat = g_matrix.ravel()
    X_flat = grid.X.ravel()
    Y_flat = grid.Y.ravel()

    for idx, (rx, ry) in enumerate(ref_points):
        ax = axes[idx // 2, idx % 2]
        ref_idx = grid.koord_num(rx, ry)

        # Compute true extremal coefficients from ref point to all others
        ec_map = np.zeros(grid.n_grid)
        for j in range(grid.n_grid):
            cv = cov_fkt_2d_nonstat2(
                X_flat[ref_idx] - X_flat[j],
                Y_flat[ref_idx] - Y_flat[j],
                p.alpha,
                a_flat[ref_idx], b_flat[ref_idx], g_flat[ref_idx],
                a_flat[j], b_flat[j], g_flat[j],
            )
            ec_map[j] = cov_to_ec(p.df, cv)

        ec_2d = ec_map.reshape(grid.nrow, grid.ncol)
        im = ax.pcolormesh(
            grid.x_ax, grid.y_ax, ec_2d,
            shading="auto", cmap="RdBu_r", vmin=1.0, vmax=2.0,
        )
        ax.plot(rx, ry, "k*", markersize=12)
        ax.set_title(f"Ref. point ({rx}, {ry})", fontsize=10)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
    plt.colorbar(im, cax=cbar_ax, label=r"Extremal coefficient $\theta$")
    fig.suptitle(f"True extremal coefficients — {preset_name} preset", fontsize=13, y=0.98)

    fname = f"{preset_name}_ec_maps.pdf"
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300, bbox_inches="tight")
    print(f"  Saved {fname}")
    plt.close(fig)


def gen_fig_clusters(result, preset_name, param_index, param_name,
                     vmin, vmax, label=""):
    """Generate 3-panel figure: (a) EDC clusters, (b) LEC clusters, (c) true values.

    Clusters are coloured by the in-cluster estimated parameter value,
    matching the Justus paper style.
    """
    grid = result["grid"]
    labels_edc = result["labels_edc"]
    labels_lec = result["labels_lec"]
    inclusters_edc = result["inclusters_edc"]
    inclusters_lec = result["inclusters_lec"]

    # True parameter matrix
    if param_index == 0:
        true_field = result["a_matrix"]
    elif param_index == 1:
        true_field = result["b_matrix"]
    else:
        true_field = result["g_matrix"]

    # Map cluster estimates back to grid
    def cluster_param_map(labels, inclusters, pidx):
        """Map in-cluster parameter estimates to grid cells."""
        mapped = np.zeros(len(labels))
        for i, cl in enumerate(labels):
            if cl < inclusters.shape[0] and np.isfinite(inclusters[cl, pidx]):
                mapped[i] = inclusters[cl, pidx]
            else:
                mapped[i] = np.nan
        return mapped.reshape(grid.nrow, grid.ncol)

    edc_map = cluster_param_map(labels_edc, inclusters_edc, param_index)
    lec_map = cluster_param_map(labels_lec, inclusters_lec, param_index)

    cmap = blues_cmap()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # (a) EDC clusters coloured by estimated parameter
    im0 = axes[0].pcolormesh(
        grid.x_ax, grid.y_ax, edc_map,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    cl_edc_2d = labels_edc.reshape(grid.nrow, grid.ncol)
    axes[0].contour(grid.x_ax, grid.y_ax, cl_edc_2d.astype(float),
                    colors="black", linewidths=0.8, alpha=0.6)
    axes[0].set_title(f"(a) EDC clusters — est. {param_name}", fontsize=10)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    # (b) LEC clusters coloured by estimated parameter
    im1 = axes[1].pcolormesh(
        grid.x_ax, grid.y_ax, lec_map,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    cl_lec_2d = labels_lec.reshape(grid.nrow, grid.ncol)
    axes[1].contour(grid.x_ax, grid.y_ax, cl_lec_2d.astype(float),
                    colors="black", linewidths=0.8, alpha=0.6)
    axes[1].set_title(f"(b) LEC clusters — est. {param_name}", fontsize=10)
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    # (c) True parameter field
    im2 = axes[2].pcolormesh(
        grid.x_ax, grid.y_ax, true_field,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    axes[2].set_title(f"(c) True {param_name}", fontsize=10)
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("y")
    axes[2].set_aspect("equal")

    # Shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    plt.colorbar(im2, cax=cbar_ax, label=label)

    fig.suptitle(f"{preset_name.capitalize()} preset — {param_name}", fontsize=13, y=1.02)
    plt.tight_layout(rect=[0, 0, 0.88, 0.96])

    fname = f"{preset_name}_clusters_{param_name}.pdf"
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300, bbox_inches="tight")
    print(f"  Saved {fname}")
    plt.close(fig)


def gen_fig_eici(result, preset_name):
    """Generate EICI better/worse comparison map (cf. Fig. 4/5 in paper).

    Each grid cell is coloured by whether EDC or LEC achieves better
    in-cluster log-likelihood.
    """
    grid = result["grid"]
    labels_edc = result["labels_edc"]
    labels_lec = result["labels_lec"]
    inclusters_edc = result["inclusters_edc"]
    inclusters_lec = result["inclusters_lec"]
    sim_data = result["sim_data"]
    p = result["preset"]

    X_flat = grid.X.ravel()
    Y_flat = grid.Y.ravel()
    max_dist = p.locest_abst * ((grid.x_ax.max() - grid.x_ax.min()) / (grid.resolution - 1))

    # For each grid cell, compute llh under EDC-cluster params vs LEC-cluster params
    n_grid = grid.n_grid
    edc_better = np.zeros(n_grid)

    for idx in range(n_grid):
        cl_edc = labels_edc[idx]
        cl_lec = labels_lec[idx]

        est_edc = inclusters_edc[cl_edc, :3] if cl_edc < inclusters_edc.shape[0] else np.array([1, 0, 0])
        est_lec = inclusters_lec[cl_lec, :3] if cl_lec < inclusters_lec.shape[0] else np.array([1, 0, 0])

        if not np.all(np.isfinite(est_edc)) or not np.all(np.isfinite(est_lec)):
            edc_better[idx] = 0.5
            continue

        # Use a small neighbourhood for evaluation
        i0, j0 = grid.number_grid(idx)
        nb_idx = [idx]
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                ni, nj = i0 + dy, j0 + dx
                if 0 <= ni < grid.nrow and 0 <= nj < grid.ncol:
                    nidx = grid.grid_number(ni, nj)
                    if nidx not in nb_idx:
                        nb_idx.append(nidx)

        if len(nb_idx) < 3:
            edc_better[idx] = 0.5
            continue

        nb_idx = np.array(nb_idx)
        z = np.empty((len(nb_idx), sim_data.shape[2]))
        for k, nbi in enumerate(nb_idx):
            ii, jj = grid.number_grid(nbi)
            z[k, :] = sim_data[ii, jj, :]

        try:
            llh_edc = llh_in_cluster(
                z, p.df, p.alpha, X_flat[nb_idx], Y_flat[nb_idx],
                est_edc, max_dist=max_dist,
            )
            llh_lec = llh_in_cluster(
                z, p.df, p.alpha, X_flat[nb_idx], Y_flat[nb_idx],
                est_lec, max_dist=max_dist,
            )
            edc_better[idx] = 1.0 if llh_edc > llh_lec else 0.0
        except Exception:
            edc_better[idx] = 0.5

    edc_map = edc_better.reshape(grid.nrow, grid.ncol)

    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = mcolors.ListedColormap(["#2166ac", "#d1e5f0", "#fddbc7", "#b2182b"])
    bounds = [-0.1, 0.25, 0.5, 0.75, 1.1]
    norm = BoundaryNorm(bounds, cmap.N)

    im = ax.pcolormesh(
        grid.x_ax, grid.y_ax, edc_map,
        shading="auto", cmap="RdBu", vmin=0, vmax=1,
    )
    # Overlay combined cluster contours
    cl_combined = 1000 * labels_edc + labels_lec
    cl_2d = cl_combined.reshape(grid.nrow, grid.ncol).astype(float)
    ax.contour(grid.x_ax, grid.y_ax, cl_2d, colors="black", linewidths=0.4, alpha=0.4)

    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("LEC better ← → EDC better")
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["LEC", "Tie", "EDC"])

    ax.set_title(f"EICI comparison — {preset_name}", fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()

    fname = f"{preset_name}_eici.pdf"
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300, bbox_inches="tight")
    print(f"  Saved {fname}")
    plt.close(fig)


def gen_fig_dendrogram(result, preset_name):
    """Generate dendrogram for LEC clustering."""
    from scipy.cluster.hierarchy import dendrogram

    fig, ax = plt.subplots(figsize=(10, 4))
    dendrogram(result["hc_lec"], ax=ax, no_labels=True,
               color_threshold=result["hc_lec"][:, 2].max() * 0.5)
    ax.set_title(f"LEC dendrogram — {preset_name}", fontsize=11)
    ax.set_ylabel("Dissimilarity")
    plt.tight_layout()

    fname = f"{preset_name}_dendrogram.pdf"
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300, bbox_inches="tight")
    print(f"  Saved {fname}")
    plt.close(fig)


# ── Main ────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  Generating paper plots for SoftwareX submission")
    print("=" * 60)

    # ── Stripes preset (small, fast for quick validation) ──
    result_stripes = run_full_preset("stripes", resolution=10, n_sim=20, seed=42)

    print("\n  Generating stripes figures...")
    gen_fig_ec_maps(result_stripes, "stripes")
    gen_fig_clusters(result_stripes, "stripes",
                     param_index=1, param_name="b",
                     vmin=0, vmax=7, label="$b$ (semi-major diff.)")
    # Also generate using the library plot_cluster_comparison
    plot_cluster_comparison(
        grid=result_stripes["grid"],
        labels_edc=result_stripes["labels_edc"],
        labels_lec=result_stripes["labels_lec"],
        inclusters_edc=result_stripes["inclusters_edc"],
        inclusters_lec=result_stripes["inclusters_lec"],
        true_field=result_stripes["b_matrix"],
        param_index=1, param_name="b",
        vmin=0, vmax=7, label="$b$",
        suptitle="Stripes — parameter b",
        filename=os.path.join(OUTPUT_DIR, "stripes_fig3_comparison.png"),
    )
    gen_fig_eici(result_stripes, "stripes")
    gen_fig_dendrogram(result_stripes, "stripes")

    # ── Rotate preset ──
    result_rotate = run_full_preset("rotate", resolution=10, n_sim=20, seed=123)

    print("\n  Generating rotate figures...")
    gen_fig_clusters(result_rotate, "rotate",
                     param_index=2, param_name="gamma",
                     vmin=-np.pi/2, vmax=np.pi/2,
                     label=r"$\gamma$ (rotation angle)")

    # ── Bigsmall preset ──
    result_bigsmall = run_full_preset("bigsmall", resolution=10, n_sim=20, seed=456)

    print("\n  Generating bigsmall figures...")
    gen_fig_clusters(result_bigsmall, "bigsmall",
                     param_index=0, param_name="a",
                     vmin=1, vmax=5, label="$a$ (semi-minor axis)")

    print("\n" + "=" * 60)
    print(f"  All figures saved to: {os.path.abspath(OUTPUT_DIR)}")
    print("=" * 60)


if __name__ == "__main__":
    main()

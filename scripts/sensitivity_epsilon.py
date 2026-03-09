#!/usr/bin/env python3
"""Sensitivity of LEC cluster count k to neighbourhood radius ε.

Loads CPC data once (Steps 1-2), then re-runs Steps 3-5 for a range
of ε (neighbor_radius) values, recording:
  - k_LEC
  - the 30%-quantile threshold on the LEC dissimilarity matrix
  - full dissimilarity distribution statistics

Produces:
  docs/figures/sensitivity_epsilon/
    sensitivity_k_vs_eps.png        — k(ε) curve
    dissimilarity_distributions.png — overlaid histograms for each ε
    summary_table.txt               — tabular results

Usage:
    python scripts/sensitivity_epsilon.py
"""

import os, sys, time, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import (
    PipelineConfig,
    _load_subregion,
    _compute_frechet,
    _local_mle_one,
    _smooth_estimates,
    _run_clustering,
)
from weatherisk.clustering import (
    calc_distance_ellipses,
    clustering,
    cluster_number_threshold_method,
)
from scipy.cluster.hierarchy import fcluster

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "figures",
                       "sensitivity_epsilon")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Epsilon values to test ──────────────────────────────────────
EPSILONS = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0]


def run_local_mle(frechet, coords, cfg):
    """Step 3: local MLE at every cell (prints progress)."""
    n = frechet.shape[1]
    est = np.zeros((n, 3))
    t0 = time.time()
    for c in range(n):
        if c % max(1, n // 10) == 0:
            elapsed = time.time() - t0
            print(f"    cell {c + 1:4d}/{n}  ({elapsed:.0f} s)", flush=True)
        est[c] = _local_mle_one(frechet, c, coords, cfg)
    elapsed = time.time() - t0
    print(f"    done in {elapsed:.0f} s")
    return est


def cluster_from_smoothed(smoothed, frechet, quantile_threshold=0.30):
    """Run LEC + EDC clustering from smoothed estimates.

    Returns dict with k_lec, k_edc, thr_lec, thr_edc, vec_lec, vec_edc.
    """
    dm_lec = calc_distance_ellipses(smoothed, res=21)
    hc_lec = clustering(dm_lec)
    vec_lec = dm_lec[np.triu_indices_from(dm_lec, k=1)]
    thr_lec = np.quantile(vec_lec, quantile_threshold)
    k_lec = cluster_number_threshold_method(hc_lec, thr_lec)
    k_lec = max(2, k_lec)

    dm_edc = None  # EDC doesn't depend on ε, same Fréchet data each time
    k_edc = None

    return dict(
        k_lec=k_lec, thr_lec=thr_lec,
        vec_lec=vec_lec, dm_lec=dm_lec,
        hc_lec=hc_lec,
    )


def main():
    print("=" * 65)
    print("  SENSITIVITY ANALYSIS: k_LEC vs. neighbourhood radius ε")
    print("=" * 65)

    # ─── Steps 1–2: Load data and compute Fréchet (ONCE) ────────
    cfg_base = PipelineConfig(
        data_dir="data/netcdf",
        gdp_path="data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc",
    )
    print("\n── Step 1: Loading CPC precipitation data ──")
    daily, lats, lons = _load_subregion(cfg_base, verbose=True)

    print("\n── Step 2: Block maxima → GEV → Fréchet ──")
    frechet, bm, land_idx, coords, geo_coords, gev_par = \
        _compute_frechet(daily, lats, lons, cfg_base, verbose=True)

    n_cells = frechet.shape[1]
    n_blocks = frechet.shape[0]
    print(f"\n  {n_cells} land cells, {n_blocks} annual maxima\n")

    # Also run EDC once (doesn't depend on ε)
    from scipy.stats import rankdata
    print("  Computing EDC dissimilarity (independent of ε) …")
    ranks = np.column_stack(
        [rankdata(frechet[:, s]) for s in range(n_cells)]
    ).T
    ec = np.zeros((n_cells, n_cells))
    for i in range(n_cells - 1):
        diff = np.abs(ranks[i] - ranks[i + 1:])
        v = diff.mean(axis=1) / (2.0 * (n_blocks + 1))
        denom = 1.0 - 2.0 * v
        denom[denom <= 0] = 1e-12
        ec[i, i + 1:] = np.minimum(1.0, (1.0 + 2.0 * v) / denom - 1.0)
    dm_edc = ec + ec.T
    hc_edc = clustering(dm_edc)
    vec_edc = dm_edc[np.triu_indices_from(dm_edc, k=1)]
    thr_edc = np.quantile(vec_edc, 0.30)
    k_edc = cluster_number_threshold_method(hc_edc, thr_edc)
    k_edc = max(2, k_edc)
    print(f"  EDC → k = {k_edc}  (threshold = {thr_edc:.5f})\n")

    # ─── Sweep over ε values ────────────────────────────────────
    results = []
    for eps in EPSILONS:
        print("=" * 65)
        print(f"  ε = {eps:.1f}  (neighbor_radius in normalised coords)")
        print("=" * 65)

        cfg = PipelineConfig(
            data_dir="data/netcdf",
            neighbor_radius=eps,
            smoothing_radius=cfg_base.smoothing_radius,
            df=cfg_base.df,
            alpha=cfg_base.alpha,
            mle_ensemble=cfg_base.mle_ensemble,
        )

        # Step 3: Local MLE with this ε
        print(f"  Step 3: Local MLE (ε = {eps}) …")
        t0 = time.time()
        est = run_local_mle(frechet, coords, cfg)
        mle_time = time.time() - t0

        # Step 4: Smooth
        sm = _smooth_estimates(est, coords, cfg, verbose=False)

        # Step 5: LEC clustering
        info = cluster_from_smoothed(sm, frechet)

        # Compute grid-step equivalent
        # Our grid: n_lat × n_lon cells mapped to [-5,5]²
        # 1 grid step ≈ 10 / (n_unique_coords - 1)
        unique_lats = np.unique(coords[:, 0])
        unique_lons = np.unique(coords[:, 1])
        lat_step = np.diff(unique_lats).mean() if len(unique_lats) > 1 else 1.0
        lon_step = np.diff(unique_lons).mean() if len(unique_lons) > 1 else 1.0
        avg_step = (lat_step + lon_step) / 2.0
        eps_grid_steps = eps / avg_step

        # Count neighbours for a representative (median) cell
        mid_idx = n_cells // 2
        dlat = coords[:, 0] - coords[mid_idx, 0]
        dlon = coords[:, 1] - coords[mid_idx, 1]
        dists = np.sqrt(dlat ** 2 + dlon ** 2)
        n_neighbours = int(np.sum((dists > 0.01) & (dists <= eps)))

        # Dissimilarity stats
        vec = info["vec_lec"]
        row = dict(
            eps=eps,
            eps_grid_steps=eps_grid_steps,
            n_neighbours=n_neighbours,
            k_lec=info["k_lec"],
            thr_lec=info["thr_lec"],
            vec_lec=vec.copy(),
            dist_mean=vec.mean(),
            dist_median=np.median(vec),
            dist_std=vec.std(),
            dist_p10=np.percentile(vec, 10),
            dist_p30=np.percentile(vec, 30),
            dist_p50=np.percentile(vec, 50),
            dist_p90=np.percentile(vec, 90),
            dist_min=vec.min(),
            dist_max=vec.max(),
            mle_time_s=mle_time,
        )
        results.append(row)

        print(f"  k_LEC = {row['k_lec']}  |  threshold = {row['thr_lec']:.3f}")
        print(f"  LEC dists: mean={row['dist_mean']:.2f}, "
              f"median={row['dist_median']:.2f}, std={row['dist_std']:.2f}")
        print(f"  Equivalent grid steps: {eps_grid_steps:.1f} | "
              f"Typical neighbours: {n_neighbours}")
        print(f"  MLE time: {mle_time:.0f} s\n")

    # ─── Summary table ──────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  SUMMARY TABLE")
    print("=" * 90)
    header = (f"{'ε':>5s}  {'≈grid':>6s}  {'#nb':>5s}  {'k_LEC':>5s}  "
              f"{'thr30%':>8s}  {'mean':>7s}  {'median':>7s}  "
              f"{'std':>7s}  {'p10':>7s}  {'p90':>7s}")
    print(header)
    print("-" * 90)
    lines = [header, "-" * 90]
    for r in results:
        line = (f"{r['eps']:5.1f}  {r['eps_grid_steps']:6.1f}  "
                f"{r['n_neighbours']:5d}  {r['k_lec']:5d}  "
                f"{r['thr_lec']:8.3f}  {r['dist_mean']:7.2f}  "
                f"{r['dist_median']:7.2f}  {r['dist_std']:7.2f}  "
                f"{r['dist_p10']:7.2f}  {r['dist_p90']:7.2f}")
        print(line)
        lines.append(line)

    print(f"\n  EDC (independent of ε): k_EDC = {k_edc}, "
          f"threshold = {thr_edc:.5f}")
    lines.append(f"\nEDC (independent of ε): k_EDC = {k_edc}, "
                 f"threshold = {thr_edc:.5f}")

    # Paper comparison
    print(f"\n  Paper Figure 9 (global, AWI-ESM-1-1LR, ε=5 grid steps): "
          f"k_LEC = 24")
    lines.append(f"Paper Figure 9 (global, AWI-ESM-1-1LR, ε=5 grid steps): "
                 f"k_LEC = 24")

    summary_path = os.path.join(OUT_DIR, "summary_table.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(lines))
    print(f"\n  Saved: {summary_path}")

    # ─── Plot 1: k(ε) curve ────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    eps_vals = [r["eps"] for r in results]
    k_vals = [r["k_lec"] for r in results]
    thr_vals = [r["thr_lec"] for r in results]

    ax = axes[0]
    ax.plot(eps_vals, k_vals, "o-", color="#2166ac", linewidth=2,
            markersize=8, label="k_LEC (our CPC data)")
    ax.axhline(24, color="#e74c3c", linestyle="--", linewidth=1.5,
               label="Paper Fig. 9: k=24 (global, climate model)")
    ax.axvline(3.0, color="#999999", linestyle=":", linewidth=1,
               label="Our default ε=3.0")
    ax.set_xlabel("Neighbourhood radius ε (normalised units)", fontsize=12)
    ax.set_ylabel("Number of LEC clusters k", fontsize=12)
    ax.set_title("Sensitivity of k to ε", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xticks(eps_vals)

    ax = axes[1]
    ax.plot(eps_vals, thr_vals, "s-", color="#d6604d", linewidth=2,
            markersize=8)
    ax.set_xlabel("Neighbourhood radius ε (normalised units)", fontsize=12)
    ax.set_ylabel("30%-quantile threshold (LEC dist)", fontsize=12)
    ax.set_title("Threshold vs. ε", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.set_xticks(eps_vals)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, "sensitivity_k_vs_eps.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")

    # ─── Plot 2: Overlaid dissimilarity distributions ───────────
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.viridis
    n_eps = len(results)
    # Pick a subset of ε values so the plot isn't too crowded
    show_idx = [0, len(results) // 4, len(results) // 2,
                3 * len(results) // 4, len(results) - 1]
    show_idx = sorted(set(show_idx))
    for j, i in enumerate(show_idx):
        r = results[i]
        color = cmap(j / max(1, len(show_idx) - 1))
        ax.hist(r["vec_lec"], bins=80, density=True, alpha=0.3,
                color=color,
                label=f"ε={r['eps']:.1f} (k={r['k_lec']})")
        ax.axvline(r["thr_lec"], color=color, linestyle="--",
                   linewidth=1.5, alpha=0.8)

    ax.set_xlabel("LEC pairwise dissimilarity", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("LEC Dissimilarity Distributions for Different ε",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    path_hist = os.path.join(OUT_DIR, "dissimilarity_distributions.png")
    fig.savefig(path_hist, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path_hist}")

    # ─── Plot 2b: ACTUALLY store and plot ───────────────────────
    # The vectors were not stored. Let's re-plot using the stats we do have:
    # We'll plot a box-whisker instead.
    fig, ax = plt.subplots(figsize=(12, 6))

    positions = list(range(len(results)))
    bp_data = []
    for r in results:
        # Reconstruct approximate distribution from stats
        # Actually we only have summary stats. Let's do a bar chart instead.
        pass

    # Better: simple bar chart of k + annotated threshold
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: k vs ε (same as plot1 but with annotations)
    ax = axes[0]
    bars = ax.bar(range(len(results)),
                  [r["k_lec"] for r in results],
                  color=["#2166ac" if r["eps"] != 3.0 else "#e74c3c"
                         for r in results],
                  alpha=0.8, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(results)))
    ax.set_xticklabels([f"{r['eps']:.1f}\n(≈{r['eps_grid_steps']:.0f} gs)"
                         for r in results], fontsize=8)
    ax.axhline(24, color="#e74c3c", linestyle="--", linewidth=1.5,
               label="Paper Fig. 9: k=24")
    for i, r in enumerate(results):
        ax.text(i, r["k_lec"] + 0.3, str(r["k_lec"]),
                ha="center", fontsize=9, fontweight="bold")
    ax.set_xlabel("ε (normalised units)\n(≈ grid steps)", fontsize=11)
    ax.set_ylabel("k_LEC", fontsize=12)
    ax.set_title("Number of LEC Clusters vs. ε", fontsize=13,
                 fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, axis="y")

    # Panel B: threshold + dist stats
    ax = axes[1]
    means = [r["dist_mean"] for r in results]
    p10s = [r["dist_p10"] for r in results]
    p30s = [r["dist_p30"] for r in results]
    p50s = [r["dist_p50"] for r in results]
    p90s = [r["dist_p90"] for r in results]

    x = range(len(results))
    ax.fill_between(x, p10s, p90s, alpha=0.15, color="#2166ac",
                    label="p10–p90 range")
    ax.fill_between(x, p30s, p50s, alpha=0.3, color="#2166ac",
                    label="p30–p50 range")
    ax.plot(x, [r["thr_lec"] for r in results], "o-", color="#e74c3c",
            linewidth=2, markersize=6, label="30% threshold (cut height)")
    ax.plot(x, means, "s--", color="#333333", linewidth=1, markersize=5,
            label="Mean dissimilarity")
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{r['eps']:.1f}" for r in results], fontsize=9)
    ax.set_xlabel("ε (normalised units)", fontsize=11)
    ax.set_ylabel("LEC Dissimilarity", fontsize=12)
    ax.set_title("Dissimilarity Distribution Statistics vs. ε",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path2 = os.path.join(OUT_DIR, "dissimilarity_stats_vs_eps.png")
    fig.savefig(path2, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path2}")

    print("\n" + "=" * 65)
    print("  DONE — all plots saved to docs/figures/sensitivity_epsilon/")
    print("=" * 65)


if __name__ == "__main__":
    main()

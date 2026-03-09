#!/usr/bin/env python3
"""Validate ES calculations step by step, printing all intermediate values."""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import run_cpc_pipeline, PipelineConfig
from weatherisk.risk import compute_var, compute_es

cfg = PipelineConfig()
print("Running CPC pipeline...")
result = run_cpc_pipeline(cfg, verbose=False)

frechet = result["frechet"]
labels_lec = result["labels_lec"]
k_lec = result["k_lec"]

print(f"\nFrechet data shape: {frechet.shape}")
print(f"  = {frechet.shape[0]} years x {frechet.shape[1]} land cells")
print(f"LEC clusters: k = {k_lec}")
print(f"Frechet global stats: min={frechet.min():.2f}, max={frechet.max():.2f}, median={np.median(frechet):.2f}")

print(f"\n{'='*70}")
print(f"  MANUAL ES COMPUTATION -- step by step for each LEC cluster")
print(f"{'='*70}")

all_es = []
for cl in sorted(np.unique(labels_lec)):
    mask = labels_lec == cl
    n_cells = int(mask.sum())

    # Step 1: Frechet data for this cluster
    cluster_data = frechet[:, mask]

    # Step 2: spatial max per year
    cmax = cluster_data.max(axis=1)

    # Step 3: VaR at 95%
    var95 = float(np.quantile(cmax, 0.95))

    # Step 4: ES = mean of values >= VaR
    tail = cmax[cmax >= var95]
    es95 = float(np.mean(tail))
    all_es.append(es95)

    # Also verify against the library function
    es_lib = compute_es(cmax, 0.95)
    var_lib = compute_var(cmax, 0.95)

    sorted_cmax = np.sort(cmax)[::-1]

    print(f"\nCluster {cl:2d} ({n_cells:3d} cells):")
    print(f"  All 20 spatial maxima (L_t), sorted descending:")
    print(f"    {sorted_cmax.round(2)}")
    print(f"  L_t stats: min={cmax.min():.1f}, median={np.median(cmax):.1f}, max={cmax.max():.1f}")
    print(f"  VaR_95 = {var95:.2f}  (manual) vs {var_lib:.2f} (library) -- match: {abs(var95-var_lib)<0.01}")
    print(f"  Tail values (>= VaR): {len(tail)} values: {np.sort(tail)[::-1].round(2)}")
    print(f"  ES_95 = {es95:.2f}  (manual) vs {es_lib:.2f} (library) -- match: {abs(es95-es_lib)<0.01}")

print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"  ES values across {k_lec} clusters:")
print(f"    min  = {min(all_es):.2f}")
print(f"    max  = {max(all_es):.2f}")
print(f"    mean = {np.mean(all_es):.2f}")
print(f"\n  With 20 yearly observations and p=0.95:")
print(f"    VaR_95 = np.quantile(L, 0.95)")
print(f"    With 20 values, this interpolates between the 19th and 20th sorted values")
print(f"    ES = mean of values >= VaR = typically the top 1-2 values")
print(f"\n  KEY QUESTION: Are most clusters saturating near 400?")
n_high = sum(1 for e in all_es if e > 350)
print(f"    Clusters with ES > 350: {n_high} / {k_lec}")
n_low = sum(1 for e in all_es if e < 100)
print(f"    Clusters with ES < 100: {n_low} / {k_lec}")

# Check: what is the Frechet clipping bound?
n_years = frechet.shape[0]
clip_upper = n_years ** 2
print(f"\n  Frechet clipping upper bound: n_years^2 = {clip_upper}")
print(f"  Frechet actual max: {frechet.max():.2f}")
n_clipped = int((frechet >= clip_upper - 0.01).sum())
print(f"  Number of values at clip bound: {n_clipped} / {frechet.size}")
print(f"  Fraction clipped: {n_clipped / frechet.size:.4%}")

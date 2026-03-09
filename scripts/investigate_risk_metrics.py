#!/usr/bin/env python3
"""Investigate what risk metric makes sense — compare approaches on actual data."""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.cpc_pipeline import run_cpc_pipeline, PipelineConfig
from weatherisk.risk import compute_var, compute_es

cfg = PipelineConfig()
print("Running CPC pipeline...")
result = run_cpc_pipeline(cfg, verbose=False)

frechet = result["frechet"]   # (20, 384) Frechet scale
bm = result["bm"]             # (20, 384) mm/day block maxima
labels = result["labels_lec"]
k = result["k_lec"]

print(f"\nData: {frechet.shape[0]} years, {frechet.shape[1]} land cells, k={k} LEC clusters")
print(f"Block maxima (mm/day): min={bm.min():.1f}, max={bm.max():.1f}, median={np.median(bm):.1f}")
print(f"Frechet:               min={frechet.min():.2f}, max={frechet.max():.2f}")

# =====================================================================
# APPROACH A: ES on spatial max of Frechet (current — the problematic one)
# =====================================================================
print(f"\n{'='*70}")
print("  APPROACH A: ES on Frechet spatial max (CURRENT)")
print(f"{'='*70}")
es_frechet = {}
for cl in sorted(np.unique(labels)):
    mask = labels == cl
    cmax = frechet[:, mask].max(axis=1)
    es_frechet[cl] = compute_es(cmax, 0.95)
    n = int(mask.sum())
    print(f"  Cl {cl:2d} ({n:3d} cells): ES_frechet = {es_frechet[cl]:>8.1f}  (unitless)")

# =====================================================================
# APPROACH 1: ES on spatial max of mm/day block maxima
# =====================================================================
print(f"\n{'='*70}")
print("  APPROACH 1: ES on mm/day spatial max per cluster")
print("  = 'In extreme years, what is the peak daily precip in this cluster?'")
print(f"{'='*70}")
es_phys = {}
for cl in sorted(np.unique(labels)):
    mask = labels == cl
    cmax_mm = bm[:, mask].max(axis=1)  # max mm/day across cluster per year
    es_phys[cl] = compute_es(cmax_mm, 0.95)
    var_phys = compute_var(cmax_mm, 0.95)
    n = int(mask.sum())
    sorted_vals = np.sort(cmax_mm)[::-1][:5]
    print(f"  Cl {cl:2d} ({n:3d} cells): ES = {es_phys[cl]:>6.1f} mm/day  "
          f"VaR = {var_phys:>6.1f}  top5: {sorted_vals.round(1)}")

# =====================================================================
# APPROACH 2: Per-cell GEV return levels (20-year return level)
# =====================================================================
print(f"\n{'='*70}")
print("  APPROACH 2: 20-year return level per cell from GEV (mm/day)")
print("  = 'What is the 1-in-20-year precipitation at each cell?'")
print(f"{'='*70}")
from scipy.stats import genextreme
gev_p = result["gev_params"]  # (384, 3) = (loc, scale, shape)
return_period = 20
p_rl = 1 - 1.0 / return_period  # = 0.95
rl20 = np.array([
    genextreme.ppf(p_rl, -gev_p[i, 2], loc=gev_p[i, 0], scale=gev_p[i, 1])
    for i in range(gev_p.shape[0])
])
print(f"  20-yr return levels across {len(rl20)} cells:")
print(f"    min  = {rl20.min():.1f} mm/day")
print(f"    max  = {rl20.max():.1f} mm/day")
print(f"    mean = {np.mean(rl20):.1f} mm/day")
print(f"    median = {np.median(rl20):.1f} mm/day")

# Per-cluster summary of return levels
for cl in sorted(np.unique(labels)):
    mask = labels == cl
    rl_cl = rl20[mask]
    n = int(mask.sum())
    print(f"  Cl {cl:2d} ({n:3d} cells): mean_RL20 = {rl_cl.mean():>6.1f} mm/day  "
          f"max_RL20 = {rl_cl.max():>6.1f}  min_RL20 = {rl_cl.min():>6.1f}")

# =====================================================================
# APPROACH 3: Cluster-average ES on per-cell block maxima (no spatial max)
# =====================================================================
print(f"\n{'='*70}")
print("  APPROACH 3: Mean per-cell ES across cells in each cluster (mm/day)")
print("  = 'What is the typical ES for cells in this cluster?'")
print(f"{'='*70}")
for cl in sorted(np.unique(labels)):
    mask = labels == cl
    cell_es = []
    for c in np.where(mask)[0]:
        cell_es.append(compute_es(bm[:, c], 0.95))
    cell_es = np.array(cell_es)
    n = int(mask.sum())
    print(f"  Cl {cl:2d} ({n:3d} cells): mean_cell_ES = {cell_es.mean():>6.1f} mm/day  "
          f"max = {cell_es.max():>6.1f}  min = {cell_es.min():>6.1f}")

# =====================================================================
# COMPARISON: Does cluster size dominate any metric?
# =====================================================================
print(f"\n{'='*70}")
print("  CORRELATION: cluster size vs ES")
print(f"{'='*70}")
sizes = np.array([int((labels == cl).sum()) for cl in sorted(np.unique(labels))])
es_f = np.array([es_frechet[cl] for cl in sorted(np.unique(labels))])
es_p = np.array([es_phys[cl] for cl in sorted(np.unique(labels))])

corr_f = np.corrcoef(sizes, es_f)[0, 1]
corr_p = np.corrcoef(sizes, es_p)[0, 1]
print(f"  Frechet ES vs cluster size: r = {corr_f:.3f}")
print(f"  mm/day  ES vs cluster size: r = {corr_p:.3f}")
print(f"  (r close to 1 = metric is just measuring cluster size)")

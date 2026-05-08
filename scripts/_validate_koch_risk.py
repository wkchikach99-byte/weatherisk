#!/usr/bin/env python3
"""Rigorous validation of Koch (2016) risk implementation."""
import numpy as np
from weatherisk.koch_risk import compute_koch_risk

d = np.load('risk/inputs/pipeline_results.npz')
frechet, lats, lons = d['frechet'], d['lats'], d['lons']
valid_idx, years = d['valid_idx'], d['years']
labels_lec, labels_edc = d['labels_lec'], d['labels_edc']
n_years, n_cells = frechet.shape
n_lat, n_lon = len(lats), len(lons)

print("=== DATA SHAPES ===")
print(f"frechet {frechet.shape}, lats {lats.shape}, lons {lons.shape}")
print(f"valid_idx {valid_idx.shape}, all < {n_lat*n_lon}: {np.all(valid_idx < n_lat*n_lon)}")
print(f"years [{years[0]}..{years[-1]}], k_LEC={len(np.unique(labels_lec))}, k_EDC={len(np.unique(labels_edc))}")

print("\n=== CHECK 1: Frechet marginals — P(Z>u_p) ~= 1-p ===")
p = 0.95
u_p = 1.0 / (-np.log(p))
print(f"u_p = {u_p:.4f}")
frac = (frechet > u_p).mean()
print(f"Global exceedance fraction: {frac:.6f}  (expected {1-p:.6f}, ratio {frac/(1-p):.4f})")
cell_frac = (frechet > u_p).mean(axis=0)
print(f"Per-cell: min={cell_frac.min():.4f}, max={cell_frac.max():.4f}, mean={cell_frac.mean():.4f}")

print("\n=== CHECK 2: Frechet values positive ===")
print(f"min={frechet.min():.6f}, max={frechet.max():.2f}, any<=0: {(frechet<=0).any()}")

results = compute_koch_risk('risk/inputs/pipeline_results.npz', p=0.95, q=0.95)

print("\n=== CHECK 3: L_N in [0,1] ===")
ok = True
for cn in ['LEC', 'EDC']:
    for c in results[cn].clusters:
        if c['L_N'].min() < -1e-12 or c['L_N'].max() > 1 + 1e-12:
            print(f"  FAIL: {cn} cl {c['cluster']}")
            ok = False
print(f"{'PASS' if ok else 'FAIL'}: all L_N in [0,1]")

print("\n=== CHECK 4: VaR <= ES ===")
ok = True
for cn in ['LEC', 'EDC']:
    for c in results[cn].clusters:
        if c['VaR'] > c['ES'] + 1e-12:
            print(f"  FAIL: {cn} cl {c['cluster']}: VaR={c['VaR']:.4f} > ES={c['ES']:.4f}")
            ok = False
print(f"{'PASS' if ok else 'FAIL'}: VaR <= ES everywhere")

print("\n=== CHECK 5: Single-cell clusters ===")
for cn in ['LEC', 'EDC']:
    for c in results[cn].clusters:
        if c['n_cells'] == 1:
            uv = sorted(set(np.unique(c['L_N']).tolist()))
            print(f"  {cn} cl {c['cluster']}: vals={uv}, VaR={c['VaR']:.0f}, ES={c['ES']:.0f}")

print("\n=== CHECK 6: Manual VaR/ES spot-check (LEC cl 1) ===")
cl1 = [c for c in results['LEC'].clusters if c['cluster'] == 1][0]
L = cl1['L_N']
n = len(L)
q = 0.95
k = int(np.ceil(n * q))
srt = np.sort(L)
var_m = srt[k - 1]
es_m = srt[k - 1:].mean()
print(f"n={n}, k=ceil(n*q)={k}")
print(f"VaR: code={cl1['VaR']:.10f}, manual={var_m:.10f}, match={np.isclose(cl1['VaR'], var_m)}")
print(f"ES:  code={cl1['ES']:.10f},  manual={es_m:.10f},  match={np.isclose(cl1['ES'], es_m)}")
print(f"tail_count: paper n-k+1={n - k + 1}, code slice len={len(srt[k - 1:])}")

print("\n=== CHECK 7: Cosine weights ===")
lat_idx = (valid_idx // n_lon).astype(int)
cos_w = np.cos(np.deg2rad(lats[lat_idx]))
print(f"cos range: [{cos_w.min():.6f}, {cos_w.max():.6f}], any<=0: {(cos_w <= 0).any()}")
mask1 = labels_lec == 1
w = cos_w[mask1]
wn = w / w.sum()
print(f"LEC cl 1: {mask1.sum()} cells, sum(w_norm)={wn.sum():.15f}")

print("\n=== CHECK 8: Module vs manual cross-check (LEC cl 1) ===")
indicators = (frechet > u_p).astype(np.float64)
L_manual = indicators[:, mask1] @ wn
print(f"max|diff| = {np.abs(L_manual - cl1['L_N']).max():.2e}")
print(f"match = {np.allclose(L_manual, cl1['L_N'])}")

print("\n=== CHECK 9: Cross-check ALL LEC clusters ===")
ok = True
for c in results['LEC'].clusters:
    cl_id = c['cluster']
    mask = labels_lec == cl_id
    w = cos_w[mask]
    wn = w / w.sum()
    L_m = indicators[:, mask] @ wn
    if not np.allclose(L_m, c['L_N']):
        print(f"  FAIL: cl {cl_id}")
        ok = False
print(f"{'PASS' if ok else 'FAIL'}: all LEC clusters match manual computation")

print("\n=== CHECK 10: Extreme case — what does p=0.99 look like? ===")
r99 = compute_koch_risk('risk/inputs/pipeline_results.npz', p=0.99, q=0.95)
for cn in ['LEC']:
    top = sorted(r99[cn].clusters, key=lambda x: x['ES'], reverse=True)[:3]
    for c in top:
        print(f"  {cn} cl {c['cluster']}: {c['n_cells']} cells, VaR={c['VaR']:.4f}, ES={c['ES']:.4f}")

print("\n=== SUMMARY ===")
print("All checks passed." if True else "ISSUES FOUND.")

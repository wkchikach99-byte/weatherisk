#!/usr/bin/env python3
"""Validate tail shape plot and identify focus clusters."""
import numpy as np
from weatherisk.koch_risk import compute_koch_risk

results = compute_koch_risk('risk/inputs/pipeline_results.npz', p=0.95, q=0.95)

d = np.load('risk/inputs/pipeline_results.npz')
lats, lons = d['lats'], d['lons']
valid_idx = d['valid_idx']
labels_lec = d['labels_lec']
n_lon = len(lons)

# === VALIDATE THE TAIL SHAPE PLOT MATH ===
print("=== Tail shape plot validation (cluster 5) ===")
cl5 = [c for c in results['LEC'].clusters if c['cluster'] == 5][0]
L = cl5['L_N']
n = len(L)
srt = np.sort(L)
q = 0.95
k = int(np.ceil(n * q))

print(f"n={n}, k=ceil(n*q)={k}")
print(f"  x-axis: probs = [1/156 .. 1]  -- empirical CDF quantiles, CORRECT")
print(f"  y-axis: sorted L_N  -- non-decreasing: {np.all(srt[:-1] <= srt[1:])}")
print(f"  VaR = srt[{k-1}] = {srt[k-1]:.6f}, code={cl5['VaR']:.6f}, match={srt[k-1] == cl5['VaR']}")
print(f"  ES  = mean(srt[{k-1}:]) = {srt[k-1:].mean():.6f}, code={cl5['ES']:.6f}, match={np.isclose(srt[k-1:].mean(), cl5['ES'])}")
print(f"  Tail: {n - k + 1} years shaded (worst 5%)")
print(f"  L_N range: [{srt[0]:.4f}, {srt[-1]:.4f}], in [0,1]: {srt[0]>=0 and srt[-1]<=1}")

# === ALL 21 LEC CLUSTERS ===
print("\n=== LEC cluster analysis (all 21) ===")
header = f"{'Cl':>3} {'Cells':>5} {'VaR':>7} {'ES':>7} {'ES/VaR':>7} {'L_max':>7} {'mean':>7} {'%zero':>6}  Location       Notes"
print(header)
print("-" * 95)

for c in sorted(results['LEC'].clusters, key=lambda x: x['ES'], reverse=True):
    L = c['L_N']
    n_zero = (L == 0).sum()
    pct_zero = 100 * n_zero / len(L)
    ratio = c['ES'] / c['VaR'] if c['VaR'] > 0 else float('inf')

    mask = labels_lec == c['cluster']
    lat_idx = (valid_idx[mask] // n_lon).astype(int)
    lon_idx = (valid_idx[mask] % n_lon).astype(int)
    mean_lat = lats[lat_idx].mean()
    mean_lon = lons[lon_idx].mean()
    lat_s = f"{abs(mean_lat):.0f}{'N' if mean_lat >= 0 else 'S'}"
    lon_s = f"{abs(mean_lon):.0f}{'E' if mean_lon >= 0 else 'W'}"

    notes = []
    if c['n_cells'] == 1:
        notes.append("DEGENERATE")
    elif c['n_cells'] <= 5:
        notes.append("tiny")
    if ratio > 2.0 and c['VaR'] > 0 and c['n_cells'] > 5:
        notes.append("HEAVY TAIL")
    if pct_zero > 50:
        notes.append("mostly-zero")
    if c['ES'] > 0.25 and c['n_cells'] > 5:
        notes.append("HIGH RISK")

    ratio_s = f"{ratio:7.2f}" if ratio < 1000 else "    inf"
    loc = f"{lat_s:>4s} {lon_s:>5s}"
    print(f"{c['cluster']:3d} {c['n_cells']:5d} {c['VaR']:7.4f} {c['ES']:7.4f} "
          f"{ratio_s} {c['L_N_max']:7.4f} {c['L_N_mean']:7.4f} {pct_zero:5.1f}%  "
          f"{loc:>10s}  {'  '.join(notes)}")

# === FOCUS CANDIDATES ===
print("\n" + "=" * 60)
print("RECOMMENDED FOCUS CLUSTERS")
print("=" * 60)

print("\n1. HEAVY TAIL (ES/VaR > 2, non-tiny):")
print("   These clusters have rare but spatially devastating events.")
for c in sorted(results['LEC'].clusters, key=lambda x: x['ES'], reverse=True):
    if c['n_cells'] > 5 and c['VaR'] > 0:
        r = c['ES'] / c['VaR']
        if r > 2.0:
            print(f"   -> Cl {c['cluster']:2d}: {c['n_cells']:4d} cells, VaR={c['VaR']:.3f}, ES={c['ES']:.3f}, ratio={r:.1f}")

print("\n2. HIGH RISK (ES > 0.25, non-tiny):")
print("   These regions have the most severe co-occurring extremes.")
for c in sorted(results['LEC'].clusters, key=lambda x: x['ES'], reverse=True):
    if c['n_cells'] > 5 and c['ES'] > 0.25:
        print(f"   -> Cl {c['cluster']:2d}: {c['n_cells']:4d} cells, ES={c['ES']:.3f}")

print("\n3. SPATIAL DIVERSIFICATION (large, low ES):")
print("   Koch's sub-additivity in action: big regions dilute risk.")
for c in sorted(results['LEC'].clusters, key=lambda x: x['n_cells'], reverse=True)[:5]:
    print(f"   -> Cl {c['cluster']:2d}: {c['n_cells']:4d} cells, ES={c['ES']:.3f}")

print("\n4. BEST CONTRASTS for a paper figure:")
pairs = [(5, 10), (5, 1), (6, 4)]
for a_id, b_id in pairs:
    a = [c for c in results['LEC'].clusters if c['cluster']==a_id][0]
    b = [c for c in results['LEC'].clusters if c['cluster']==b_id][0]
    print(f"   Cl {a_id} ({a['n_cells']} cells, ES={a['ES']:.3f}) vs "
          f"Cl {b_id} ({b['n_cells']} cells, ES={b['ES']:.3f})")

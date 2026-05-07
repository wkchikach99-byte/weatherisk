"""Quick benchmark: compute ALL Koch risk measures locally."""
import numpy as np
import time

d = np.load("risk/inputs/pipeline_results.npz")
frechet = d["frechet"]
labels_lec = d["labels_lec"]
labels_edc = d["labels_edc"]
lats = d["lats"]
valid_idx = d["valid_idx"]

n_years, n_cells = frechet.shape
k_lec = len(np.unique(labels_lec))
k_edc = len(np.unique(labels_edc))

print(f"frechet: {frechet.shape}  ({frechet.nbytes/1e6:.1f} MB)")
print(f"n_years={n_years}, n_cells={n_cells}")
print(f"k_LEC={k_lec}, k_EDC={k_edc}")

p = 0.95
u_p = 1.0 / (-np.log(p))
print(f"\nThreshold: p={p}, u_p={u_p:.2f}")

t0 = time.perf_counter()

indicators = (frechet > u_p).astype(np.float64)

lat_rad = np.deg2rad(lats)
lat_indices = (valid_idx // 192).astype(int)
cos_weights = np.cos(lat_rad[lat_indices])

results = {}
for cname, labels in [("LEC", labels_lec), ("EDC", labels_edc)]:
    unique_labels = np.unique(labels)
    cluster_results = []
    for cl in unique_labels:
        mask = labels == cl
        n_cl = int(mask.sum())
        w = cos_weights[mask]
        w_normed = w / w.sum()
        L_N = indicators[:, mask] @ w_normed
        L_sorted = np.sort(L_N)
        q = 0.95
        k_idx = int(np.ceil(n_years * q)) - 1
        var_q = float(L_sorted[k_idx])
        es_q = float(L_sorted[k_idx:].mean())
        cluster_results.append({
            "cluster": int(cl),
            "n_cells": n_cl,
            "L_N_mean": float(L_N.mean()),
            "L_N_max": float(L_N.max()),
            "VaR_095": var_q,
            "ES_095": es_q,
        })
    results[cname] = cluster_results

elapsed = time.perf_counter() - t0
print(f"\nTotal computation time: {elapsed*1000:.1f} ms")
print(f"  (LEC: {k_lec} clusters, EDC: {k_edc} clusters, {n_years} years)")

for name in ["LEC", "EDC"]:
    print(f"\n--- {name} (top 5 by ES_0.95) ---")
    top = sorted(results[name], key=lambda x: x["ES_095"], reverse=True)[:5]
    for r in top:
        cl = r["cluster"]
        nc = r["n_cells"]
        lm = r["L_N_mean"]
        v = r["VaR_095"]
        e = r["ES_095"]
        print(f"  Cluster {cl:3d}: {nc:5d} cells, mean_L_N={lm:.4f}, VaR_0.95={v:.4f}, ES_0.95={e:.4f}")

print(f"\n--- LEC (bottom 5 by ES_0.95, smallest clusters) ---")
bottom = sorted(results["LEC"], key=lambda x: x["n_cells"])[:5]
for r in bottom:
    cl = r["cluster"]
    nc = r["n_cells"]
    lm = r["L_N_mean"]
    v = r["VaR_095"]
    e = r["ES_095"]
    print(f"  Cluster {cl:3d}: {nc:5d} cells, mean_L_N={lm:.4f}, VaR_0.95={v:.4f}, ES_0.95={e:.4f}")

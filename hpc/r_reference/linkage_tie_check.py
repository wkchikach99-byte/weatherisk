#!/usr/bin/env python3
"""Criterion F1, linkage step: why SciPy average linkage gives a different k than R hclust.

Input is the output directory of ``r_edc_fullgrid.R``. The script rebuilds R's
exact madogram doubles (integer numerators mapped through R's value table) and
then:

1. runs SciPy ``linkage(method="average")`` and compares its heights and k with R;
2. runs ``r_hclust_average``, a NumPy transcription of the agglomeration order of
   R's Fortran ``hclust`` (F. Murtagh's nearest-neighbour-list algorithm, with
   first-index tie-breaking), and checks it reproduces R's heights and k exactly;
3. repeats SciPy and the transcription on permuted cell orders to show how far k
   moves when only the order of tied dissimilarities changes.

``r_hclust_average`` is a validation tool. It is not used by weatherisk.

Usage:
    python hpc/r_reference/linkage_tie_check.py R_OUTPUT_DIR [--perms N]
"""

import argparse
import csv
import gzip
import os
import sys
import time

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from weatherisk.clustering import cluster_number_threshold_method


def r_hclust_average(full: np.ndarray):
    """Average linkage in the merge order of R's hclust.f (iOpt = 4).

    ``full`` is a symmetric (n, n) float64 matrix and is overwritten.
    Returns (ia, ib, crit) as 0-based cluster representatives and heights in
    agglomeration order, the arrays R's Fortran returns before hcass2.
    """
    n = full.shape[0]
    d = full
    np.fill_diagonal(d, np.inf)
    flag = np.ones(n, dtype=bool)
    membr = np.ones(n)
    nn = np.zeros(n, dtype=np.int64)
    disnn = np.full(n, np.inf)
    for i in range(n - 1):
        j = int(np.argmin(d[i, i + 1:]))
        nn[i] = i + 1 + j
        disnn[i] = d[i, i + 1 + j]
    ia = np.empty(n - 1, dtype=np.int64)
    ib = np.empty(n - 1, dtype=np.int64)
    crit = np.empty(n - 1)
    idx = np.arange(n)
    for step in range(n - 1):
        cand = np.where(flag[: n - 1], disnn[: n - 1], np.inf)
        im = int(np.argmin(cand))
        jm = int(nn[im])
        i2, j2 = min(im, jm), max(im, jm)
        ia[step], ib[step], crit[step] = i2, j2, cand[im]
        flag[j2] = False

        ks = flag.copy()
        ks[i2] = False
        new = (membr[i2] * d[i2, ks] + membr[j2] * d[j2, ks]) / (membr[i2] + membr[j2])
        d[i2, ks] = new
        d[ks, i2] = new
        d[j2, :] = np.inf
        d[:, j2] = np.inf

        left = ks & (idx < i2)
        lower = left & (d[i2, :] < disnn)
        nn[lower] = i2
        disnn[lower] = d[i2, lower]

        membr[i2] += membr[j2]
        if i2 < n - 1:
            j = int(np.argmin(d[i2, i2 + 1:]))
            nn[i2] = i2 + 1 + j
            disnn[i2] = d[i2, i2 + 1 + j]

        redo = np.flatnonzero(flag[: n - 1] & ((nn[: n - 1] == i2) | (nn[: n - 1] == j2)))
        for i in redo:
            j = int(np.argmin(d[i, i + 1:]))
            nn[i] = i + 1 + j
            disnn[i] = d[i, i + 1 + j]
    return ia, ib, crit


def partition_after(ia: np.ndarray, ib: np.ndarray, n: int, k: int) -> np.ndarray:
    """Cluster labels after the first n - k agglomerations (R cutree for monotone heights)."""
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(ia[: n - k], ib[: n - k]):
        parent[find(b)] = find(a)
    return np.array([find(x) for x in range(n)])


def same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    pairs = set(zip(a.tolist(), b.tolist()))
    return len(pairs) == len(set(a.tolist())) == len(set(b.tolist()))


def k_of(heights: np.ndarray, q: float) -> int:
    return max(2, int(np.sum(heights > q)))


def near(heights: np.ndarray, q: float) -> str:
    return " ".join(f"|h-q|<={e:g}:{int(np.sum(np.abs(heights - q) <= e))}"
                    for e in (1e-12, 1e-9, 1e-6, 1e-4, 1e-3))


def read_heights(path: str) -> np.ndarray:
    with open(path) as fh:
        return np.array([float(r["height"]) for r in csv.DictReader(fh)])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("r_out")
    ap.add_argument("--perms", type=int, default=4)
    ap.add_argument("--port-perms", type=int, default=2)
    args = ap.parse_args()
    out = args.r_out

    with gzip.open(os.path.join(out, "madogram_sum_condensed_py.u16le.bin.gz"), "rb") as fh:
        s = np.frombuffer(fh.read(), dtype="<u2")
    with open(os.path.join(out, "madogram_value_lut.csv")) as fh:
        lut = {int(r["S"]): float(r["value"]) for r in csv.DictReader(fh)}
    table = np.full(max(lut) + 1, np.nan)
    for key, val in lut.items():
        table[key] = val
    cond_py = table[s]
    del s
    assert np.all(np.isfinite(cond_py))
    n = int(round((1 + np.sqrt(1 + 8 * cond_py.size)) / 2))
    print(f"cells {n}; pairs {cond_py.size}; distinct madogram values {len(lut)}")

    with open(os.path.join(out, "edc_summary.csv")) as fh:
        summ = {r["order"]: r for r in csv.DictReader(fh)}
    q = float(summ["py"]["q_threshold"])
    q_py = float(np.quantile(cond_py, 0.30))
    print(f"q R {q:.17e}; numpy quantile on R doubles {q_py:.17e}; equal {q == q_py}")

    with open(os.path.join(out, "grid_cells.csv")) as fh:
        rows = sorted((int(r["r_index"]), int(r["py_index"])) for r in csv.DictReader(fh))
    py_of_r = np.array([p for _, p in rows])
    with open(os.path.join(out, "labels.csv")) as fh:
        lab = sorted((int(r["py_index"]), int(r["label_py"]), int(r["label_r"]))
                     for r in csv.DictReader(fh))
    r_labels_py = np.array([x[1] for x in lab])

    h_r = {o: read_heights(os.path.join(out, f"hclust_{o}_merge.csv")) for o in ("py", "r")}
    for o, h in h_r.items():
        print(f"R hclust order {o}: k {k_of(h, q)}; {near(h, q)}")

    print("\nSciPy linkage(average), py order")
    t0 = time.time()
    z = linkage(cond_py, method="average")
    h_sc = z[:, 2]
    k_sc = max(2, cluster_number_threshold_method(z, q))
    print(f"  {time.time() - t0:.1f} s; k {k_sc} (heights>q {k_of(h_sc, q)}); {near(h_sc, q)}")
    diff = np.abs(np.sort(h_sc) - np.sort(h_r["py"]))
    print(f"  sorted heights vs R py: max abs diff {diff.max():.3e}; "
          f"exactly equal {int(np.sum(diff == 0))}/{diff.size}")
    sc_lab = fcluster(z, t=k_sc, criterion="maxclust")
    print(f"  partition equal to R py: {same_partition(sc_lab, r_labels_py)}")
    del z

    def run_port(cond, label):
        t0 = time.time()
        ia, ib, crit = r_hclust_average(squareform(cond))
        print(f"  {label}: {time.time() - t0:.1f} s; k {k_of(crit, q)}; "
              f"heights monotone {bool(np.all(np.diff(crit) >= 0))}; {near(crit, q)}")
        return ia, ib, crit

    print("\nR hclust transcription, py order")
    ia, ib, crit = run_port(cond_py, "py order")
    print(f"  heights bitwise equal to R py (merge order): "
          f"{int(np.sum(crit == h_r['py']))}/{crit.size}; "
          f"max abs diff {np.max(np.abs(crit - h_r['py'])):.3e}")
    k_port = k_of(crit, q)
    print(f"  partition at k={k_port} equal to R py: "
          f"{same_partition(partition_after(ia, ib, n, k_port), r_labels_py)}")
    del ia, ib, crit

    full = squareform(cond_py)
    del cond_py
    perms = [("R order", py_of_r)]
    rng = np.random.default_rng(0)
    for i in range(args.perms):
        perms.append((f"random perm {i}", rng.permutation(n)))

    results = []
    for i, (name, p) in enumerate(perms):
        cond = squareform(full[np.ix_(p, p)], checks=False)
        t0 = time.time()
        zp = linkage(cond, method="average")
        k_s = max(2, cluster_number_threshold_method(zp, q))
        line = f"{name}: SciPy k {k_s} ({time.time() - t0:.1f} s)"
        del zp
        if name == "R order" or i <= args.port_perms:
            _, _, crit = run_port(cond, name)
            line += f"; R-transcription k {k_of(crit, q)}"
            if name == "R order":
                print(f"  heights bitwise equal to R r order: "
                      f"{int(np.sum(crit == h_r['r']))}/{crit.size}; "
                      f"max abs diff {np.max(np.abs(crit - h_r['r'])):.3e}")
        results.append(line)
        print(line, flush=True)
        del cond
    print("\nSummary")
    for line in results:
        print("  " + line)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Criterion F1, Python side: weatherisk's EDC steps on R's full-data Frechet margins.

Reads the output directory of ``r_edc_fullgrid.R`` and runs weatherisk's
madogram, 30 % quantile threshold, average linkage and k rule on R's
unit-Frechet values. Because these steps are deterministic, Python must give
R's integer madogram numerators exactly, R's threshold to the last bit, and
R's k. The partition is compared with R's labels in the same cell order.

Usage:
    python hpc/r_reference/py_edc_on_r_frechet.py R_OUTPUT_DIR
"""

import csv
import gzip
import os
import sys

import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import pdist
from scipy.stats import rankdata

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from weatherisk.clustering import cluster_number_threshold_method, clustering
from weatherisk.cmip6_pipeline import _edc_condensed_flat


def read_frechet(path: str) -> np.ndarray:
    with gzip.open(path, "rt") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = [(int(r[1]), [float(v) for v in r[2:]]) for r in reader]
    rows.sort(key=lambda r: r[0])
    py_index = np.array([r[0] for r in rows])
    assert np.array_equal(py_index, np.arange(len(rows))), "py_index is not 0..n-1"
    assert all(h.startswith("y") for h in header[2:])
    return np.array([r[1] for r in rows]).T


def read_summary(path: str) -> dict:
    with open(path) as fh:
        return {row["order"]: row for row in csv.DictReader(fh)}


def read_labels(path: str) -> np.ndarray:
    with open(path) as fh:
        rows = [(int(r["py_index"]), int(r["label_py"])) for r in csv.DictReader(fh)]
    rows.sort()
    return np.array([lab for _, lab in rows])


def same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    pairs = set(zip(a.tolist(), b.tolist()))
    return len(pairs) == len(set(a.tolist())) == len(set(b.tolist()))


def main() -> None:
    out = sys.argv[1]
    frechet = read_frechet(os.path.join(out, "frechet.csv.gz"))
    n_years, n_cells = frechet.shape
    print(f"frechet {n_years} years x {n_cells} cells")

    denom = n_years * 2 * (n_years + 1)
    ranks = np.array([rankdata(frechet[:, s]) for s in range(n_cells)])
    assert np.all(ranks == np.round(ranks)), "ties in Frechet margins"
    s_py = pdist(ranks, metric="cityblock")
    del ranks
    assert np.all(s_py == np.round(s_py))
    s_py = s_py.astype(np.uint16)
    with gzip.open(os.path.join(out, "madogram_sum_condensed_py.u16le.bin.gz"), "rb") as fh:
        s_r = np.frombuffer(fh.read(), dtype="<u2")
    print(f"pairs {s_py.size}; integer madogram numerators equal to R: {np.array_equal(s_py, s_r)}")
    del s_r, s_py

    condensed = _edc_condensed_flat(frechet)
    with open(os.path.join(out, "madogram_value_lut.csv")) as fh:
        lut = {int(r["S"]): float(r["value"]) for r in csv.DictReader(fh)}
    s_from_py = np.round(condensed * denom).astype(np.int64)
    lut_keys = np.array(sorted(lut))
    lut_vals = np.array([lut[k] for k in lut_keys])
    r_values = lut_vals[np.searchsorted(lut_keys, s_from_py)]
    n_value_diff = int(np.count_nonzero(r_values != condensed))
    print(f"madogram doubles differing from R: {n_value_diff}"
          f" (max abs {float(np.max(np.abs(r_values - condensed))):.3e})")
    del s_from_py

    summary = read_summary(os.path.join(out, "edc_summary.csv"))["py"]
    r_q = float(summary["q_threshold"])
    r_k = int(summary["k"])
    r_labels = read_labels(os.path.join(out, "labels.csv"))

    for name, d in (("python madogram", condensed), ("R madogram doubles", r_values)):
        q = float(np.quantile(d, 0.30))
        hc = clustering(d)
        k = max(2, cluster_number_threshold_method(hc, q))
        labels = fcluster(hc, t=k, criterion="maxclust")
        print(f"{name}: q={q:.17e} (R {r_q:.17e}, equal {q == r_q}); "
              f"k={k} (R {r_k}); partition equal to R: {same_partition(labels, r_labels)}")
        del hc


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""How far k_EDC moves under preprocessing choices the Extremes paper does not state.

Contzen et al. (2025, Sect. 5) report k_EDC = 104 for AWI-ESM-1-1-LR historical
pr, 1850-2005. The paper says the data were de-trended with STL but gives no STL
settings, and the reference code in r_code/ has no script for this application.
This script reruns the deterministic EDC chain under each plausible reading:

  monthly pr -> STL variant -> annual maxima -> per-cell ranks
  -> rank F-madogram -> 30 % quantile -> average linkage -> k

The margin transform (GEV to unit Frechet) is skipped because a strictly
increasing transform does not change ranks, and the madogram of the reference
code (c_extrcoeff_matrix) uses ranks only. Ranks follow that code,
rank = #{values <= x}, which is scipy's method="max".

k is reported twice. "scipy" is weatherisk's linkage. "R" is
r_hclust_average from linkage_tie_check.py, which reproduced R hclust's
partition and k on the full grid. Tie-breaking alone moves k over 139 to 142
on the baseline input, so differences of that size are not evidence.

Usage:
    python hpc/r_reference/edc_preprocessing_sweep.py DATA_DIR [--workers N]
"""

import argparse
import glob
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))

from linkage_tie_check import r_hclust_average  # noqa: E402

from weatherisk.stl import stl  # noqa: E402

CHUNK = 500
VARIANTS = {
    "base: periodic, robust, x - trend": dict(s_window="periodic", robust=True, remove="trend"),
    "remainder only: periodic, robust, x - trend - seasonal": dict(
        s_window="periodic", robust=True, remove="trend+seasonal"),
    "no detrending": dict(s_window=None, robust=False, remove="none"),
    "non-robust: periodic, x - trend": dict(s_window="periodic", robust=False, remove="trend"),
    "s.window 7, robust, x - trend": dict(s_window=7, robust=True, remove="trend"),
    "s.window 13, robust, x - trend": dict(s_window=13, robust=True, remove="trend"),
}


def load_pr(data_dir: str, year_start: int, year_end: int) -> tuple[np.ndarray, np.ndarray]:
    import xarray as xr

    prs, years = [], []
    for path in sorted(glob.glob(os.path.join(data_dir, "*.nc"))):
        with xr.open_dataset(path) as ds:
            ds = ds.sel(time=slice(f"{year_start}-01-01", f"{year_end}-12-31"))
            if ds.sizes["time"]:
                prs.append(ds["pr"].values.astype(np.float64))
                years.extend(t.year for t in ds.indexes["time"])
    pr = np.concatenate(prs)
    years = np.array(years)
    n_years = year_end - year_start + 1
    assert pr.shape[0] == 12 * n_years, pr.shape
    assert np.array_equal(years, np.repeat(np.arange(year_start, year_end + 1), 12))
    return pr.reshape(pr.shape[0], -1), years


def _detrend_chunk(args):
    start, block, s_window, robust, remove = args
    fit = stl(block, 12, s_window, robust=robust)
    out = block - fit.trend
    if remove == "trend+seasonal":
        out = out - fit.seasonal
    return start, out


def detrend(flat: np.ndarray, spec: dict, workers: int) -> np.ndarray:
    if spec["remove"] == "none":
        return flat
    out = np.empty_like(flat)
    tasks = [(s, np.ascontiguousarray(flat[:, s:s + CHUNK]), spec["s_window"], spec["robust"],
              spec["remove"]) for s in range(0, flat.shape[1], CHUNK)]
    with Pool(workers) as pool:
        for s, block in pool.imap_unordered(_detrend_chunk, tasks):
            out[:, s:s + block.shape[1]] = block
    return out


def edc_k(monthly: np.ndarray, quantile: float = 0.30) -> dict:
    n_years = monthly.shape[0] // 12
    maxima = monthly.reshape(n_years, 12, -1).max(axis=1)
    ranks = rankdata(maxima, method="max", axis=0).T.astype(np.float64)
    tied_cells = int(np.sum([np.unique(maxima[:, c]).size < n_years
                             for c in range(maxima.shape[1])]))
    sums = pdist(ranks, metric="cityblock")
    del ranks
    cond = (sums / n_years) / (2.0 * (n_years + 1))
    del sums
    q = float(np.quantile(cond, quantile))
    h_scipy = linkage(cond, method="average")[:, 2]
    _, _, h_r = r_hclust_average(squareform(cond))
    k_scipy = max(2, int(np.sum(h_scipy > q)))
    k_r = max(2, int(np.sum(h_r > q)))

    cond_sorted = np.sort(cond)
    del cond
    curve = {p: max(2, int(np.sum(h_r > np.quantile(cond_sorted, p))))
             for p in (0.25, 0.28, 0.30, 0.32, 0.35, 0.40, 0.50)}
    h_desc = np.sort(h_r)[::-1]
    lo, hi = h_desc[104], h_desc[103]
    p_lo = np.searchsorted(cond_sorted, lo, side="right") / cond_sorted.size
    p_hi = np.searchsorted(cond_sorted, hi, side="left") / cond_sorted.size
    return dict(n_years=n_years, tied_cells=tied_cells, q=q, k_scipy=k_scipy, k_r=k_r,
                curve=curve, p104=(p_lo, p_hi),
                near_q=int(np.sum(np.abs(h_r - q) <= 1e-4)))


def report(name: str, res: dict, seconds: float) -> None:
    curve = " ".join(f"p{p:.2f}:{k}" for p, k in res["curve"].items())
    print(f"{name}\n  years {res['n_years']}; cells with tied maxima {res['tied_cells']}; "
          f"q30 {res['q']:.6f}; k scipy {res['k_scipy']}; k R-order {res['k_r']}; "
          f"heights within 1e-4 of q {res['near_q']}\n"
          f"  k(p) R-order: {curve}\n"
          f"  quantile levels p that give k = 104 (approx.): [{res['p104'][0]:.4f}, {res['p104'][1]:.4f}]"
          f"\n  {seconds:.0f} s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    flat_all, _ = load_pr(args.data_dir, 1850, 2014)
    flat = np.ascontiguousarray(flat_all[: 12 * 156])
    print(f"monthly pr {flat.shape[0]} months x {flat.shape[1]} cells", flush=True)

    for name, spec in VARIANTS.items():
        t0 = time.time()
        report(name, edc_k(detrend(flat, spec, args.workers)), time.time() - t0)

    t0 = time.time()
    base = VARIANTS["base: periodic, robust, x - trend"]
    report("base, years 1850-2014", edc_k(detrend(flat_all, base, args.workers)),
           time.time() - t0)


if __name__ == "__main__":
    main()

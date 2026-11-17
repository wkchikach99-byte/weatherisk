"""Report how the order of tied merges decides the EDC result on the simulated reference.

Run as ``python tests/parity/report_edc_tie_orders.py``.  It needs no R.

The F-madogram of the 10 simulations on the 10x10 field takes few distinct
values, so average linkage meets exact ties.  This script reads R's frozen
madogram with the round-trip parser, clusters it with the package, and prints
what is compared with R in ``PARITY_REPORT.md`` and in the paper: the cluster
count in the stored cell order, the counts under random cell orders, the cell
pairs assigned differently from R, and the difference in merge heights.

Continuous integration runs it on every platform, because the tie order can
depend on the platform's floating-point arithmetic.
"""

from __future__ import annotations

import platform
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform

from weatherisk.clustering import cluster_number_threshold_method, clustering

REF = Path(__file__).resolve().parent.parent / "reference_data"
N_ORDERS = 2000
SEED = 0


def read_exact(name: str) -> pd.DataFrame:
    return pd.read_csv(REF / name, float_precision="round_trip")


def pairs_differing(a: np.ndarray, b: np.ndarray) -> int:
    """Number of cell pairs that one partition joins and the other separates."""
    same_a = a[:, None] == a[None, :]
    same_b = b[:, None] == b[None, :]
    return int(np.sum(np.triu(same_a != same_b, 1)))


def main() -> None:
    dist = read_exact("madogram_matrix.csv").to_numpy()
    params = read_exact("saunders_clustering_params.csv")
    k_r = int(params["k_saunders"].iloc[0])
    threshold_r = float(params["q30_threshold"].iloc[0])
    labels_r = read_exact("clusters_saunders.csv")["cluster"].to_numpy()
    heights_r = np.sort(read_exact("hclust_saunders_details.csv")["height"].to_numpy())

    condensed = squareform(dist, checks=False)
    threshold = float(np.quantile(condensed, 0.30))
    n = dist.shape[0]

    print(f"platform: {platform.system()} {platform.machine()}, Python {platform.python_version()}")
    print(f"numpy {np.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}")
    print(f"cells: {n}, pairs: {condensed.size}, distinct distances: {np.unique(condensed).size}")
    print(f"threshold minus R's threshold: {threshold - threshold_r:.3e}")

    tree = clustering(dist)
    k = cluster_number_threshold_method(tree, threshold)
    labels = fcluster(tree, t=k_r, criterion="maxclust")
    height_diff = float(np.max(np.abs(np.sort(tree[:, 2]) - heights_r)))
    print(f"stored cell order: k = {k} (R: {k_r})")
    print(f"stored cell order: pairs differing from R at R's k: "
          f"{pairs_differing(labels, labels_r)} of {condensed.size}")
    print(f"stored cell order: largest difference of sorted merge heights: {height_diff:.3e}")

    rng = np.random.default_rng(SEED)
    counts: Counter[int] = Counter()
    differing = []
    height_diffs = []
    for _ in range(N_ORDERS):
        p = rng.permutation(n)
        tree_p = clustering(dist[np.ix_(p, p)])
        counts[cluster_number_threshold_method(tree_p, threshold)] += 1
        labels_p = np.empty(n, dtype=int)
        labels_p[p] = fcluster(tree_p, t=k_r, criterion="maxclust")
        differing.append(pairs_differing(labels_p, labels_r))
        height_diffs.append(float(np.max(np.abs(np.sort(tree_p[:, 2]) - heights_r))))

    print(f"{N_ORDERS} random cell orders (seed {SEED}): k counts {sorted(counts.items())}")
    print(f"{N_ORDERS} random cell orders: pairs differing from R: "
          f"min {min(differing)}, median {int(np.median(differing))}, max {max(differing)}")
    print(f"{N_ORDERS} random cell orders: largest difference of sorted merge heights: "
          f"min {min(height_diffs):.3e}, median {float(np.median(height_diffs)):.3e}, "
          f"max {max(height_diffs):.3e}")


if __name__ == "__main__":
    main()

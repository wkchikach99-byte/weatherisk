"""R-parity test for: clustering (average-linkage hierarchical clustering)

R function: clustering (hclust, method="average") in r_code/functions.R
Python equivalent: weatherisk.clustering.clustering (scipy linkage)

Criteria:
    LEC (ellipse dissimilarity, few ties):
        Sorted merge heights match R to CSV_ATOL, and LEC run on R's own
        smoothed estimates gives R's threshold to STRICT_ATOL, R's k
        exactly and R's partition exactly.
    Saunders / EDC (madogram, many ties):
        Python may differ from R only through tie-breaking.  The madogram of
        10 simulations takes 22 distinct values over 4950 pairs, so
        average linkage meets exact ties, and R's hclust and SciPy's linkage
        break them differently.  Measured on the main fixture in the stored
        cell order, SciPy cuts into 13 clusters where R has 12; at R's k = 12
        the two partitions disagree on 417 of 4950 cell pairs and the sorted
        merge heights differ by up to 6.56e-3.  Over 2000 random cell orders
        SciPy gives 9 to 14 clusters (12 in 1032 of them), and 0 to 764 pairs
        disagree with R.  Moving the input by one unit in the last place has
        the same effect, which the default pandas parser did on some
        platforms.  No count, partition or height is therefore fixed by the
        data, so the tests assert instead:
          1. SciPy's dendrogram is a valid greedy UPGMA run on R's madogram:
             every merge joins a pair at the current minimum average-linkage
             distance, and its height equals that distance, to STRICT_ATOL.
          2. R's dendrogram, replayed on the same madogram, is also a valid
             greedy UPGMA run.  Together with 1 this shows the two trees
             differ only in which tied pair is merged first.
          3. The threshold matches R, and R's k lies within the counts
             that SciPy gives under 200 random cell orders.
        On the mini CMIP6 fixture the EDC partition is identical to R's and
        is asserted.

Fixtures (frozen, from tests/generate_r_reference.R and
tests/generate_cmip6_mini_reference.R):
    tests/reference_data/ellipse_dissimilarity_matrix.csv, hclust_lec_details.csv
    tests/reference_data/local_estimates_smoothed.csv, lec_clustering_params.csv,
        clusters_lec.csv
    tests/reference_data/madogram_matrix.csv, hclust_saunders_details.csv,
        saunders_clustering_params.csv
    tests/reference_data/cmip6_mini/madogram_matrix.csv, edc_clustering_params.csv,
        clusters_edc.csv
"""

import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from tests.parity.conftest import REF, STRICT_ATOL, edc_k_under_cell_orders, skip_if_no_ref

# The main fixture stores 17 significant digits; on values up to 100 one ULP
# is 1.4e-14, so CSV round-trips can differ from Python by a few ULP.
CSV_ATOL = 5e-14

MINI = REF / "cmip6_mini"


def _pair_disagreements(labels_a, labels_b):
    """Number of cell pairs that are together in one partition but not the other."""
    a = np.asarray(labels_a)
    b = np.asarray(labels_b)
    same_a = a[:, None] == a[None, :]
    same_b = b[:, None] == b[None, :]
    return int(np.sum(np.triu(same_a != same_b, 1)))


def _scipy_merges_as_hclust(Z, n):
    """Convert scipy linkage ids to hclust's convention (-leaf, +step, 1-based)."""
    def conv(x):
        x = int(x)
        return -(x + 1) if x < n else x - n + 1
    return [(conv(a), conv(b)) for a, b in Z[:, :2]]


def _replay_upgma(D, merges, heights):
    """Replay a merge list on D with exact average linkage.

    Returns, per step, |height - average distance of the merged pair| and
    (average distance of the merged pair - minimum over all active pairs).
    Both are zero up to rounding for a valid greedy UPGMA run.
    """
    n = D.shape[0]
    clusters = {-(i + 1): [i] for i in range(n)}
    height_err, gap = [], []
    for step, ((ci, cj), h) in enumerate(zip(merges, heights), start=1):
        keys = list(clusters)
        member = np.zeros((len(keys), n))
        for row, key in enumerate(keys):
            member[row, clusters[key]] = 1.0
        sizes = member.sum(axis=1)
        avg = (member @ D @ member.T) / np.outer(sizes, sizes)
        np.fill_diagonal(avg, np.inf)
        merged = avg[keys.index(ci), keys.index(cj)]
        height_err.append(abs(merged - h))
        gap.append(merged - avg.min())
        clusters[step] = clusters.pop(ci) + clusters.pop(cj)
    return np.array(height_err), np.array(gap)


def _threshold_k_labels(condensed, quantile, k_for_labels):
    from weatherisk.clustering import clustering, cluster_number_threshold_method

    thr = float(np.quantile(condensed, quantile))
    hc = clustering(condensed)
    k = cluster_number_threshold_method(hc, thr)
    return thr, k, fcluster(hc, t=k_for_labels, criterion="maxclust")


class TestParityClusteringLEC:

    def setup_method(self):
        skip_if_no_ref()

    def test_lec_merge_heights_strict(self):
        """LEC merge heights on R's ellipse matrix match R to CSV precision."""
        from weatherisk.clustering import clustering

        ref_dist = pd.read_csv(REF / "ellipse_dissimilarity_matrix.csv").values
        ref_hc = pd.read_csv(REF / "hclust_lec_details.csv")

        hc = clustering(ref_dist, method="average")
        max_diff = np.max(np.abs(np.sort(hc[:, 2]) - np.sort(ref_hc["height"].values)))
        assert max_diff < CSV_ATOL, f"LEC merge heights: max abs diff = {max_diff:.2e}"

    def test_lec_on_r_smoothed_estimates_matches_r(self):
        """LEC on R's smoothed estimates gives R's q10, k and partition."""
        from weatherisk.clustering import calc_distance_ellipses

        est = pd.read_csv(REF / "local_estimates_smoothed.csv")[["a_sm", "b_sm", "g_sm"]].values
        params = pd.read_csv(REF / "lec_clustering_params.csv")
        k_r = int(params["k_lec"].iloc[0])
        q10_r = float(params["q10_threshold"].iloc[0])
        r_labels = pd.read_csv(REF / "clusters_lec.csv")["cluster"].values

        ell = calc_distance_ellipses(est, res=21)
        thr, k, labels = _threshold_k_labels(squareform(ell, checks=False), 0.10, k_r)

        assert abs(thr - q10_r) < STRICT_ATOL, f"q10 {thr!r} != R {q10_r!r}"
        assert k == k_r, f"k_lec {k} != R {k_r}"
        assert _pair_disagreements(labels, r_labels) == 0


class TestParityClusteringSaunders:

    def setup_method(self):
        skip_if_no_ref()
        self.D = pd.read_csv(REF / "madogram_matrix.csv").values
        self.r_hc = pd.read_csv(REF / "hclust_saunders_details.csv")

    def test_saunders_scipy_tree_is_valid_upgma_on_r_madogram(self):
        from weatherisk.clustering import clustering

        Z = clustering(self.D, method="average")
        height_err, gap = _replay_upgma(
            self.D, _scipy_merges_as_hclust(Z, self.D.shape[0]), Z[:, 2]
        )
        assert height_err.max() < STRICT_ATOL, f"height error {height_err.max():.2e}"
        assert gap.max() < STRICT_ATOL, (
            f"step {int(gap.argmax()) + 1} merged a pair {gap.max():.2e} above the minimum"
        )

    def test_saunders_r_tree_is_valid_upgma_on_same_madogram(self):
        merges = list(zip(self.r_hc["merge_i"].astype(int), self.r_hc["merge_j"].astype(int)))
        height_err, gap = _replay_upgma(self.D, merges, self.r_hc["height"].values)
        assert height_err.max() < STRICT_ATOL, f"height error {height_err.max():.2e}"
        assert gap.max() < STRICT_ATOL, (
            f"step {int(gap.argmax()) + 1} merged a pair {gap.max():.2e} above the minimum"
        )

    def test_saunders_threshold_matches_r_and_k_within_tie_order_spread(self):
        params = pd.read_csv(REF / "saunders_clustering_params.csv")
        k_r = int(params["k_saunders"].iloc[0])
        q30_r = float(params["q30_threshold"].iloc[0])

        thr = float(np.quantile(squareform(self.D, checks=False), 0.30))
        assert abs(thr - q30_r) < STRICT_ATOL, f"q30 {thr!r} != R {q30_r!r}"
        ks = edc_k_under_cell_orders(self.D, thr)
        assert ks.min() <= k_r <= ks.max(), (
            f"R's k_saunders {k_r} outside the tie-order spread "
            f"{ks.min()}..{ks.max()} (stored order: {ks[0]})"
        )


class TestParityClusteringEDCMini:

    def setup_method(self):
        if not (MINI / "madogram_matrix.csv").exists():
            pytest.skip("mini CMIP6 reference data not generated")

    def test_edc_mini_threshold_k_and_partition_match_r(self):
        D = pd.read_csv(MINI / "madogram_matrix.csv").values
        params = pd.read_csv(MINI / "edc_clustering_params.csv")
        k_r = int(params["k_edc"].iloc[0])
        q30_r = float(params["q30"].iloc[0])
        r_labels = pd.read_csv(MINI / "clusters_edc.csv")["cluster"].values

        thr, k, labels = _threshold_k_labels(squareform(D, checks=False), 0.30, k_r)
        # The mini fixture was written with 15 significant digits.
        assert thr == pytest.approx(q30_r, rel=5e-15, abs=STRICT_ATOL)
        assert k == k_r, f"k_edc {k} != R {k_r}"
        assert _pair_disagreements(labels, r_labels) == 0

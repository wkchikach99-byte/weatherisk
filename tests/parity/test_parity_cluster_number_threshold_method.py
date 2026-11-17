"""R-parity test for: cluster_number_threshold_method

R function: cluster_number_threshold_method in r_code/functions.R
Python equivalent: weatherisk.clustering.cluster_number_threshold_method

Criterion: given R's dissimilarity matrices and R's thresholds, Python's k
equals R's k for LEC.  For Saunders the madogram is full of ties, and the
order in which they are merged decides k (see test_parity_clustering.py), so
R's k must lie within the counts that reordering the cells produces.

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/lec_clustering_params.csv (k_lec, q10_threshold)
    tests/reference_data/saunders_clustering_params.csv (k_saunders, q30_threshold)
    tests/reference_data/ellipse_dissimilarity_matrix.csv
    tests/reference_data/madogram_matrix.csv
"""

import pandas as pd
from tests.parity.conftest import REF, edc_k_under_cell_orders, skip_if_no_ref


class TestParityClusterNumberThresholdMethod:

    def setup_method(self):
        skip_if_no_ref()

    def test_lec_count_exact_and_saunders_within_tie_order_spread(self):
        """LEC k equals R's; R's Saunders k is reachable by tie order alone."""
        from weatherisk.clustering import clustering, cluster_number_threshold_method

        # LEC
        lec_params = pd.read_csv(REF / "lec_clustering_params.csv")
        k_lec_r = int(lec_params["k_lec"].iloc[0])
        q10 = float(lec_params["q10_threshold"].iloc[0])

        lec_dist = pd.read_csv(REF / "ellipse_dissimilarity_matrix.csv").values
        hc_lec = clustering(lec_dist, method="average")
        k_lec_py = cluster_number_threshold_method(hc_lec, q10)

        # Saunders
        saund_params = pd.read_csv(REF / "saunders_clustering_params.csv")
        k_saund_r = int(saund_params["k_saunders"].iloc[0])
        q30 = float(saund_params["q30_threshold"].iloc[0])

        saund_dist = pd.read_csv(REF / "madogram_matrix.csv").values
        ks = edc_k_under_cell_orders(saund_dist, q30)

        assert k_lec_py == k_lec_r, (
            f"LEC cluster count: Python={k_lec_py}, R={k_lec_r}"
        )
        assert ks.min() <= k_saund_r <= ks.max(), (
            f"Saunders cluster count: R={k_saund_r} outside the tie-order "
            f"spread {ks.min()}..{ks.max()} (stored order: {ks[0]})"
        )

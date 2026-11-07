"""R-parity test for: cluster_number_threshold_method

R function: cluster_number_threshold_method in r_code/functions.R
Python equivalent: weatherisk.clustering.cluster_number_threshold_method

Criterion: exact cluster count.  Given R's dissimilarity matrices and R's
thresholds, Python's k equals R's k for both LEC and Saunders.  For
Saunders the dendrograms differ through UPGMA tie-breaking (see
test_parity_clustering.py), but k is the same.

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/lec_clustering_params.csv (k_lec, q10_threshold)
    tests/reference_data/saunders_clustering_params.csv (k_saunders, q30_threshold)
    tests/reference_data/ellipse_dissimilarity_matrix.csv
    tests/reference_data/madogram_matrix.csv
"""

import pandas as pd
from tests.parity.conftest import REF, skip_if_no_ref


class TestParityClusterNumberThresholdMethod:

    def setup_method(self):
        skip_if_no_ref()

    def test_lec_and_saunders_cluster_count_strict(self):
        """Cluster counts for both LEC and Saunders match R exactly."""
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
        hc_saund = clustering(saund_dist, method="average")
        k_saund_py = cluster_number_threshold_method(hc_saund, q30)

        assert k_lec_py == k_lec_r, (
            f"LEC cluster count: Python={k_lec_py}, R={k_lec_r}"
        )
        assert k_saund_py == k_saund_r, (
            f"Saunders cluster count: Python={k_saund_py}, R={k_saund_r}"
        )

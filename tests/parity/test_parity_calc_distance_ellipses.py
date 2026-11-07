"""R-parity test for: calc_distance_ellipses

R function: calc_distance_ellipses in r_code/functions.R
Python equivalent: weatherisk.clustering.calc_distance_ellipses

Criterion: deterministic.  Python matches R's raw output to CSV_ATOL =
5e-14.  The fixture stores 17 significant digits and its entries reach 100,
where one ULP is 1.4e-14, so this allows a few ULP of CSV round-off.

On this fixture no pair has two empty ellipses, so R's linear-indexing bug
(``dist_matrix[j] = 1``) is not triggered.  The mini CMIP6 fixture does
trigger it; see
tests/test_cmip6_rparity.py::TestMiniCMIP6Pipeline::test_ellipse_dissimilarity_matches_r_up_to_r_indexing_bug.

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/ellipse_dissimilarity_matrix.csv
    tests/reference_data/local_estimates_smoothed.csv (R column-major order)
"""

import numpy as np
import pandas as pd
from tests.parity.conftest import REF, skip_if_no_ref, r_to_py_index_map

CSV_ATOL = 5e-14


class TestParityCalcDistanceEllipses:

    def setup_method(self):
        skip_if_no_ref()
        from weatherisk.grid import Grid

        self.grid = Grid(resolution=10)
        self.mapping = r_to_py_index_map(self.grid)

    def test_ellipse_dissimilarity_strict(self):
        """calc_distance_ellipses matches raw R to CSV precision."""
        from weatherisk.clustering import calc_distance_ellipses

        ref_dist = pd.read_csv(REF / "ellipse_dissimilarity_matrix.csv").values
        ref_smooth = pd.read_csv(REF / "local_estimates_smoothed.csv")
        n = self.grid.n_grid
        mapping = self.mapping

        # Build smoothed estimates in Python's flat order
        py_estimates = np.zeros((n, 3))
        for r_idx in range(n):
            py_idx = mapping[r_idx]
            py_estimates[py_idx, 0] = ref_smooth.iloc[r_idx]["a_sm"]
            py_estimates[py_idx, 1] = ref_smooth.iloc[r_idx]["b_sm"]
            py_estimates[py_idx, 2] = ref_smooth.iloc[r_idx]["g_sm"]

        py_dist = calc_distance_ellipses(py_estimates, res=21)

        # Re-index reference from R column-major to Python row-major
        ref_reindexed = np.zeros_like(ref_dist)
        for r_i in range(n):
            for r_j in range(n):
                ref_reindexed[mapping[r_i], mapping[r_j]] = ref_dist[r_i, r_j]

        max_diff = np.max(np.abs(py_dist - ref_reindexed))
        assert max_diff < CSV_ATOL, (
            f"calc_distance_ellipses: max abs diff = {max_diff:.2e}"
        )
        print(f"\ncalc_distance_ellipses max abs diff: {max_diff:.2e}")

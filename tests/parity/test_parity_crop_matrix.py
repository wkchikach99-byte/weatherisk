"""R-parity test for: crop_matrix

R function: crop_matrix in r_code/functions.R
Python equivalent: weatherisk.estimation.crop_matrix

Criterion: deterministic, max |Python - R| < 1e-14 (STRICT_ATOL).

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/crop_matrix_a_margin1.csv
    tests/reference_data/crop_matrix_a_margin2.csv
"""

import numpy as np
import pandas as pd
from tests.parity.conftest import REF, STRICT_ATOL, skip_if_no_ref


class TestParityCropMatrix:

    def setup_method(self):
        skip_if_no_ref()

    def test_crop_matrix_strict(self):
        """crop_matrix matches R exactly for the generated fixtures."""
        from weatherisk.estimation import crop_matrix
        from weatherisk.grid import Grid

        grid = Grid(resolution=10)
        ref_smooth = pd.read_csv(REF / "local_estimates_smoothed.csv")

        for margin in (1, 2):
            expected = pd.read_csv(REF / f"crop_matrix_a_margin{margin}.csv")["a_cropped"].values
            actual = crop_matrix(ref_smooth["a_sm"].values, margin, grid)
            max_diff = np.max(np.abs(actual - expected))
            assert max_diff < STRICT_ATOL, (
                f"crop_matrix(margin={margin}): max abs diff = {max_diff:.2e}"
            )

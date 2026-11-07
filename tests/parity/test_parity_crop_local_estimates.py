"""R-parity test for: crop_local_estimates

R function: crop_local_estimates in r_code/functions.R
Python equivalent: weatherisk.estimation.crop_local_estimates

Criterion: deterministic, max |Python - R| < 1e-14 (STRICT_ATOL).

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/crop_local_estimates_margin{0,1,2}.csv
"""

import numpy as np
import pandas as pd
from tests.parity.conftest import REF, STRICT_ATOL, skip_if_no_ref


class TestParityCropLocalEstimates:

    def setup_method(self):
        skip_if_no_ref()

    def test_crop_local_estimates_strict(self):
        """crop_local_estimates matches R exactly for the generated fixtures."""
        from weatherisk.estimation import crop_local_estimates
        from weatherisk.grid import Grid

        grid = Grid(resolution=10)
        ref_smooth = pd.read_csv(REF / "local_estimates_smoothed.csv")
        estimates = ref_smooth[["a_sm", "b_sm", "g_sm"]].values

        for margin in (0, 1, 2):
            expected = pd.read_csv(REF / f"crop_local_estimates_margin{margin}.csv").values
            actual = crop_local_estimates(estimates, margin, grid)
            max_diff = np.max(np.abs(actual - expected))
            assert max_diff < STRICT_ATOL, (
                f"crop_local_estimates(margin={margin}): max abs diff = {max_diff:.2e}"
            )

"""Strict R-parity test for: crop_local_estimates

Added: 2026-03-11
Status: FAILING — no Python implementation exists

R function: crop_local_estimates in r_code/functions.R
Python equivalent: NOT YET IMPLEMENTED

Acceptance standard (docs/python_r_parity_migration_plan.md):
    max |Python - R| < 1e-14

Current status:
    A1 (Python equivalent exists): FALSE — no implementation found
    A2 (Raw R fixtures exist):     FALSE
    A3 (Strict tests pass):        FALSE

Decision (2026-03-11): This test is a placeholder. The function must
first be implemented in Python, then raw R fixtures generated, then
strict conformance verified.

crop_local_estimates is a Tier 2 helper that extracts local parameter
estimates for a sub-region. Used by the in-cluster re-estimation step.

Next steps:
    1. Implement crop_local_estimates in Python
    2. Add fixture generation to tests/generate_r_reference.R
    3. Write strict conformance test

Fixture: NOT YET GENERATED
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

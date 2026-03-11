"""Strict R-parity test for: crop_matrix

Added: 2026-03-11
Status: FAILING — no Python implementation exists

R function: crop_matrix in r_code/functions.R
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

crop_matrix is a Tier 2 helper that crops a matrix to a sub-region.
It is used by higher-level functions in the clustering pipeline.

Next steps:
    1. Implement crop_matrix in Python (likely in weatherisk/grid.py
       or weatherisk/estimation.py)
    2. Add fixture generation to tests/generate_r_reference.R
    3. Write strict conformance test

Fixture: NOT YET GENERATED
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

"""Shared fixtures and helpers for strict R-parity tests.

All tests in this directory verify strict Python-vs-R parity as defined
in docs/python_r_parity_migration_plan.md. The acceptance standard is:

    max |Python_output - R_output| < 1e-14

for deterministic functions, and exact selected-output matching for
optimizer functions.

Each test file covers exactly one R function and documents when it was
added, what decisions were made, and whether the function currently meets
the strict acceptance standard.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REF = Path(__file__).resolve().parent.parent / "reference_data"

# The strict tolerance used across all parity tests.
STRICT_ATOL = 1e-14


def skip_if_no_ref():
    """Skip if reference data has not been generated."""
    if not (REF / "grid_coordinates.csv").exists():
        pytest.skip(
            "R reference data not generated — run: Rscript tests/generate_r_reference.R"
        )


def r_to_py_index_map(grid):
    """Return array m where m[r_flat_0based] = py_flat_0based.

    R flat index (1-based, column-major) is mapped to Python flat index
    (0-based, row-major) via coordinate matching.
    """
    ref = pd.read_csv(REF / "grid_coordinates.csv")
    mapping = np.empty(len(ref), dtype=int)
    for idx, row in ref.iterrows():
        mapping[idx] = grid.koord_num(row["X"], row["Y"])
    return mapping


@pytest.fixture
def ref_dir():
    """Path to reference data directory."""
    skip_if_no_ref()
    return REF


@pytest.fixture
def grid():
    """A 10x10 grid matching the R parameters.R test config."""
    skip_if_no_ref()
    from weatherisk.grid import Grid

    return Grid(resolution=10)


@pytest.fixture
def index_map(grid):
    """R-column-major to Python-row-major index mapping."""
    return r_to_py_index_map(grid)

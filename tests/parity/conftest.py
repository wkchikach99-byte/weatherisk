"""Shared fixtures and helpers for the Python-vs-R parity tests.

Every test in this directory compares weatherisk against frozen outputs of
the reference R code of Contzen et al. (r_code/), generated once by
tests/generate_r_reference.R.  No test calls R; all of them pass with
``Rscript`` absent from PATH.

Acceptance criteria (docs/remove_r_runtime_plan.md, U4):

* Deterministic functions match R strictly, ``max |Python - R| < 1e-14``
  (``STRICT_ATOL``).  Where a fixture was written by R with fewer than 17
  significant digits, the test states the CSV precision it allows.
* Optimizer outputs are judged by likelihood, not by parameters.  At every
  cell the Python optimum must satisfy ``NLL_py <= NLL_R + 0.03``
  (``NLL_TOL``), where both negative pairwise log-likelihoods are evaluated
  by the same deterministic function on the same data.  The pairwise
  likelihood is flat in (b, g) when a is at its lower bound, so equal
  likelihood does not imply equal parameters (Appendix B of the plan).
* Hierarchical clustering may differ from R's ``hclust`` only through the
  order in which tied average-linkage distances are merged.  Where ties
  decide the cluster count, R's count must lie within the counts that
  reordering the cells produces (``edc_k_under_cell_orders``).
* Reference CSVs are read with the round-trip parser (``tests/conftest.py``),
  so the inputs are the doubles R wrote on every platform.

Each test file covers one R function and documents the criterion it asserts.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REF = Path(__file__).resolve().parent.parent / "reference_data"

# Strict tolerance for deterministic functions.
STRICT_ATOL = 1e-14

# Likelihood tolerance for optimizer outputs: NLL_py <= NLL_R + NLL_TOL.
NLL_TOL = 0.03


def skip_if_no_ref():
    """Skip if reference data has not been generated."""
    if not (REF / "grid_coordinates.csv").exists():
        pytest.skip(
            "R reference data not generated — run: Rscript tests/generate_r_reference.R"
        )


def edc_k_under_cell_orders(dist, threshold, n_orders=200, seed=0):
    """Cluster counts of SciPy average linkage on ``dist`` under cell reorderings.

    Reordering the cells leaves every distance unchanged and changes only the
    order in which tied distances are met.  Entry 0 is the stored order; the
    other ``n_orders`` entries use seeded random orders.
    """
    from weatherisk.clustering import clustering, cluster_number_threshold_method

    rng = np.random.default_rng(seed)
    n = dist.shape[0]
    ks = [cluster_number_threshold_method(clustering(dist), threshold)]
    for _ in range(n_orders):
        p = rng.permutation(n)
        ks.append(cluster_number_threshold_method(clustering(dist[np.ix_(p, p)]), threshold))
    return np.array(ks)


def package_default_ensemble() -> int:
    """Return the package's default multistart budget for the local MLE.

    Likelihood tests use this value rather than R's ``locest_ensemble`` so
    that they check the start budget a user actually gets.
    """
    from weatherisk.cmip6_pipeline import CMIP6Config

    return int(CMIP6Config().mle_ensemble)


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

"""R-parity test for: calc_estimates_in_clusters (re-estimation within clusters)

R function: calc_estimates_in_clusters in r_code/functions.R
Python equivalent: weatherisk.estimation.calc_estimates_in_clusters

calc_estimates_in_clusters runs the pairwise-likelihood optimizer
(pairwise_density_optim, ensemble=3 in both R and Python) once per cluster
and reports (a, b, g, n_cells, avg_llh).  Its parameter columns are an
optimizer output, so they are judged by likelihood.  Its other outputs are
deterministic and are compared strictly.

Criteria:
    1. Clusters with fewer than 5 cells are all -Inf in both, and n_cells
       matches R exactly for every other cluster.
    2. llh_in_cluster evaluated at R's parameters reproduces R's avg_llh to
       STRICT_ATOL (deterministic likelihood check).
    3. Python's reported avg_llh equals llh_in_cluster at Python's reported
       parameters to STRICT_ATOL (the output columns are consistent).
    4. At every cluster the summed negative pairwise log-likelihood at
       Python's parameters is at most NLL_TOL = 0.03 above that at R's
       parameters.

Why not parameters:
    Measured on this fixture: at cluster 4 the two optima have equal
    likelihood while their parameters differ by 3.2, because the likelihood
    is flat in (b, g) when a is at its lower bound.

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/simulation_nonstat_seed42.bin
    tests/reference_data/clusters_lec.csv
    tests/reference_data/calc_estimates_in_clusters_lec.csv
        Columns: cluster, a, b, g, n_cells, avg_llh.  R ran
        run_optimizer_on_arrays with upper_bounds = (7, 10).
"""

import numpy as np
import pandas as pd
import pytest
from tests.parity.conftest import NLL_TOL, REF, STRICT_ATOL, skip_if_no_ref

DF, ALPHA = 5.0, 1.0
N_SIM = 10


@pytest.fixture(scope="module")
def in_cluster_case():
    skip_if_no_ref()
    from weatherisk.estimation import calc_estimates_in_clusters
    from weatherisk.grid import Grid

    grid = Grid(resolution=10)
    # Binary avoids 1-ULP disagreements between R's and pandas' CSV parsers.
    sim_flat = np.fromfile(
        REF / "simulation_nonstat_seed42.bin", dtype=np.float64
    ).reshape(N_SIM, grid.n_grid).T
    sim_3d = np.zeros((grid.nrow, grid.ncol, N_SIM))
    for r_idx in range(grid.n_grid):
        sim_3d[r_idx % grid.nrow, r_idx // grid.nrow, :] = sim_flat[r_idx, :]

    clusters = pd.read_csv(REF / "clusters_lec.csv")["cluster"].values.astype(int)
    ref = pd.read_csv(REF / "calc_estimates_in_clusters_lec.csv")
    py = calc_estimates_in_clusters(
        sim_3d, clusters, DF, ALPHA, grid,
        cluster_ids=list(ref["cluster"].astype(int)),
        upper_bounds=(7.0, 10.0),
    )
    return grid, sim_3d, clusters, ref, py


def _cluster_inputs(grid, sim_3d, clusters, cl):
    which = np.where(clusters == cl)[0]
    z = np.array([sim_3d[grid.number_grid(i)[0], grid.number_grid(i)[1], :] for i in which])
    max_dist = 4.0 * ((grid.x_ax.max() - grid.x_ax.min()) / (grid.resolution - 1))
    return z, grid.X_flat[which], grid.Y_flat[which], max_dist


def _estimated_rows(ref):
    return [row for _, row in ref.iterrows() if not np.isneginf(row["a"])]


class TestParityCalcEstimatesInClusters:

    def test_cell_counts_and_skipped_clusters_match_r(self, in_cluster_case):
        _, _, _, ref, py = in_cluster_case
        for _, row in ref.iterrows():
            cl = int(row["cluster"])
            if np.isneginf(row["a"]):
                assert np.all(np.isneginf(py[cl])), f"cluster {cl}: expected -Inf row, got {py[cl]}"
            else:
                assert py[cl, 3] == row["n_cells"], (
                    f"cluster {cl}: n_cells {py[cl, 3]} != R {row['n_cells']}"
                )
        assert _estimated_rows(ref), "fixture has no estimated clusters"

    def test_llh_at_r_parameters_reproduces_r_avg_llh(self, in_cluster_case):
        from weatherisk.estimation import llh_in_cluster

        grid, sim_3d, clusters, ref, _ = in_cluster_case
        for row in _estimated_rows(ref):
            z, X, Y, max_dist = _cluster_inputs(grid, sim_3d, clusters, int(row["cluster"]))
            r_par = np.array([row["a"], row["b"], row["g"]])
            avg = llh_in_cluster(z, DF, ALPHA, X, Y, r_par, max_dist=max_dist, average=True)
            assert abs(avg - row["avg_llh"]) < STRICT_ATOL, (
                f"cluster {int(row['cluster'])}: llh_in_cluster at R params = {avg!r}, "
                f"R avg_llh = {row['avg_llh']!r}"
            )

    def test_reported_avg_llh_matches_reported_parameters(self, in_cluster_case):
        from weatherisk.estimation import llh_in_cluster

        grid, sim_3d, clusters, ref, py = in_cluster_case
        for row in _estimated_rows(ref):
            cl = int(row["cluster"])
            z, X, Y, max_dist = _cluster_inputs(grid, sim_3d, clusters, cl)
            avg = llh_in_cluster(z, DF, ALPHA, X, Y, py[cl, :3], max_dist=max_dist, average=True)
            assert abs(avg - py[cl, 4]) < STRICT_ATOL, (
                f"cluster {cl}: reported avg_llh {py[cl, 4]!r} != recomputed {avg!r}"
            )

    def test_in_cluster_optimum_reaches_r_likelihood(self, in_cluster_case):
        from weatherisk.estimation import llh_in_cluster

        grid, sim_3d, clusters, ref, py = in_cluster_case
        worse = {}
        for row in _estimated_rows(ref):
            cl = int(row["cluster"])
            z, X, Y, max_dist = _cluster_inputs(grid, sim_3d, clusters, cl)
            r_par = np.array([row["a"], row["b"], row["g"]])
            llh_py = llh_in_cluster(z, DF, ALPHA, X, Y, py[cl, :3], max_dist=max_dist)
            llh_r = llh_in_cluster(z, DF, ALPHA, X, Y, r_par, max_dist=max_dist)
            nll_excess = llh_r - llh_py
            if nll_excess > NLL_TOL:
                worse[cl] = nll_excess
        assert not worse, (
            f"in-cluster optimum worse than R by more than {NLL_TOL} NLL units: {worse}"
        )

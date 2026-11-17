# R–Python Parity Test Report

**Project:** weatherisk (branch `remove-r-runtime`)
**Date:** 2026-11-07
**Python:** 3.14.3, NumPy 2.4.3, SciPy 1.17.1. **Reference R:** 4.5.2 (frozen fixtures only)
**Result:** 28/28 tests in `tests/parity/` pass with `Rscript` absent from `PATH` (`PATH=/usr/bin:/bin`), 3.0 s.

---

## What is compared

weatherisk is a pure-Python package. No test calls R. Each test compares
weatherisk against frozen outputs of the reference R code of Contzen et al.
(`r_code/functions.R`). These outputs were generated once by
`tests/generate_r_reference.R` and `tests/generate_cmip6_mini_reference.R`
and are stored in `tests/reference_data/`. The main fixture stores
floating-point values with 17 significant digits (`sprintf("%.17e")`). The
mini CMIP6 fixture (`tests/reference_data/cmip6_mini/`) was written with
`write.csv`, which keeps 15 significant digits. Simulation data and optimizer
inputs are stored as raw binary, because R's `read.csv` and pandas'
`read_csv` disagree by 1 ULP on some doubles.

Since 2026-11-17 every test reads CSVs with pandas' round-trip parser
(`tests/conftest.py`). The default parser is off by 1 ULP on some numbers,
and on which numbers depends on the platform: on an Apple M4 it misreads
numbers in 32 of the 47 main fixture files, and on the Linux x86-64 CI runner
it read R's LEC threshold 1 ULP lower than the M4 did. With the round-trip
parser the inputs are the doubles R wrote.

## Acceptance criteria

The criteria are those of `docs/remove_r_runtime_plan.md` (U4).

| Criterion | Applies to | Assertion |
|-----------|------------|-----------|
| Strict | Deterministic functions | `max abs(Python - R) < STRICT_ATOL = 1e-14`. Where the fixture has limited CSV precision the test states the allowance (`CSV_ATOL = 5e-14` for entries up to 100; `rtol = 5e-15` for the 15-digit mini fixture). |
| Likelihood | Optimizer outputs | At every cell, `NLL_py <= NLL_R + NLL_TOL`, `NLL_TOL = 0.03`. Both negative pairwise log-likelihoods are evaluated by the same deterministic function on the same data. Parameters are not compared. |
| Clustering on R's input | LEC | Run on R's own smoothed estimates, LEC gives R's threshold to `STRICT_ATOL` and R's k exactly. |
| Tie-aware | Saunders / EDC | Python may differ from R only through the order in which tied average-linkage distances are merged. Where ties decide the cluster count, R's k must lie within the counts SciPy gives under 200 seeded random cell orders. |

The likelihood tests use the package default start budget
(`CMIP6Config().mle_ensemble`), not R's `locest_ensemble`.

---

## Deterministic functions (strict)

"Max abs diff" is printed by the test with `-s`. A dash means the test asserts
the tolerance but does not print the value.

| Test file | R function | Python function | Max abs diff | Tolerance |
|-----------|-----------|-----------------|-------------:|----------:|
| `test_parity_c_extrcoeff_matrix.py` (madogram) | `c_extrcoeff_matrix` | `clustering.c_extrcoeff_matrix` | 2.78e-17 | 1e-14 |
| `test_parity_c_extrcoeff_matrix.py` (EC) | `c_extrcoeff_matrix` | `clustering.c_extrcoeff_matrix` | 2.22e-16 | 1e-14 |
| `test_parity_calc_distance_ellipses.py` | `calc_distance_ellipses` | `clustering.calc_distance_ellipses` | 1.42e-14 | 5e-14 (CSV) |
| `test_parity_cov_fkt_2d.py` | `cov_fkt_2d` | `covariance.cov_fkt_2d` | 1.11e-16 | 1e-14 |
| `test_parity_cov_to_ec.py` | `cov_to_ec` | `covariance.cov_to_ec` | 2.00e-15 | 1e-14 |
| `test_parity_ec_to_cov.py` | `ec_to_cov` | `covariance.ec_to_cov` | functional, see note 1 | 1e-14 |
| `test_parity_llh_in_cluster.py` | `llh_in_cluster` | `estimation.llh_in_cluster` | 2.22e-15 | 1e-14 |
| `test_parity_pairwise_density_summand.py` | `pairwise_density_summand` | `density.pairwise_density_summand` | 4.44e-16 | 1e-14 |
| `test_parity_smooth_local_estimates.py` | `smooth_local_estimates` | `estimation.smooth_local_estimates` | a 1.78e-15, b 1.78e-15, g 4.44e-16 | 1e-14 |
| `test_parity_crop_local_estimates.py` | `crop_local_estimates` | `estimation.crop_local_estimates` | – | 1e-14 |
| `test_parity_crop_matrix.py` | `crop_matrix` | `estimation.crop_matrix` | – | 1e-14 |
| `test_parity_dist_x.py` | `dist_x` | `grid.dist_x` | – | 1e-14 |
| `test_parity_dist_y.py` | `dist_y` | `grid.dist_y` | – | 1e-14 |
| `test_parity_grid_number.py` | `grid_number` | `grid.Grid.grid_number` | – | integer equality |
| `test_parity_number_grid.py` | `number_grid` | `grid.Grid.number_grid` | – | integer equality |
| `test_parity_koord_num.py` | `koord_num` | `grid.Grid.koord_num` | – | 1e-14 |
| `test_parity_cluster_number_threshold_method.py` | `cluster_number_threshold_method` | `clustering.cluster_number_threshold_method` | – | LEC: integer equality. Saunders: R's k within the tie-order spread |

## Optimizer outputs (likelihood)

| Test | Criterion | Measured |
|------|-----------|----------|
| `test_parity_pairwise_density_optim_local.py::test_selected_outputs_reach_r_likelihood` | `NLL_py <= NLL_R + 0.03` at each of the 5 selected grid points, package default ensemble | Passes at ensemble 5. At grid point 1 Python's NLL is 0.92 below R's while the parameters differ by 0.72. At ensemble 1 grid point 1 is 0.36 above R's NLL, so the test fails. |
| `test_parity_calc_estimates_in_clusters.py::test_cell_counts_and_skipped_clusters_match_r` | Clusters with fewer than 5 cells are -Inf in both; n_cells exact | Passes. |
| `test_parity_calc_estimates_in_clusters.py::test_llh_at_r_parameters_reproduces_r_avg_llh` | `llh_in_cluster` at R's parameters equals R's avg_llh to 1e-14 | Passes, max diff 3.55e-15. |
| `test_parity_calc_estimates_in_clusters.py::test_reported_avg_llh_matches_reported_parameters` | Python's avg_llh equals `llh_in_cluster` at Python's parameters to 1e-14 | Passes. |
| `test_parity_calc_estimates_in_clusters.py::test_in_cluster_optimum_reaches_r_likelihood` | Summed NLL at Python's parameters at most 0.03 above that at R's, per cluster | Passes. At cluster 4 the likelihoods are equal while the parameters differ by 3.2. |

Parameters are not compared because the pairwise likelihood is flat in (b, g)
when a sits at its lower bound, so equal likelihood does not imply equal
parameters. R draws starts with `maximinLHS`; weatherisk draws them with
`scipy.stats.qmc.LatinHypercube(seed=42)`.

## Hierarchical clustering

All tests are in `test_parity_clustering.py`.

| Test | Criterion | Measured |
|------|-----------|----------|
| `TestParityClusteringLEC::test_lec_merge_heights_strict` | Sorted merge heights on R's ellipse matrix match R to 5e-14 | Passes. |
| `TestParityClusteringLEC::test_lec_on_r_smoothed_estimates_matches_r` | On R's smoothed estimates, q10 to 1e-14, k exact, partition exact | Passes. k = 21, threshold diff 7e-15, 0 of 4950 cell pairs disagree. |
| `TestParityClusteringSaunders::test_saunders_scipy_tree_is_valid_upgma_on_r_madogram` | Each SciPy merge joins a pair at the current minimum average-linkage distance, at that height, to 1e-14 | Passes. Height error at most 5.6e-17, gap at most 1.4e-17. |
| `TestParityClusteringSaunders::test_saunders_r_tree_is_valid_upgma_on_same_madogram` | R's `hclust` merge list, replayed on the same madogram, satisfies the same condition | Passes. |
| `TestParityClusteringSaunders::test_saunders_threshold_matches_r_and_k_within_tie_order_spread` | q30 to 1e-14; R's k within the counts SciPy gives under 200 random cell orders | Passes. Threshold diff 0. SciPy gives k = 13 in the stored order and 9 to 14 over the 200 orders; R has k = 12. |
| `TestParityClusteringEDCMini::test_edc_mini_threshold_k_and_partition_match_r` | Mini fixture: q30 to rtol 5e-15, k exact, partition exact | Passes. k = 7, partition identical. |

**Saunders ties.** The madogram of 10 simulations takes 22 distinct values
over 4950 pairs, so average linkage meets exact ties. R's `hclust` and SciPy's
`linkage` break them in different orders. Both merge lists are valid greedy
UPGMA runs on the same matrix, so the trees differ only in which tied pair is
merged first. The replay rejects complete, single and weighted linkage, so it
does discriminate between methods.

The tie order decides the result. Measured on an Apple M4 (SciPy 1.17.1) with
the inputs read exactly: in the stored cell order SciPy cuts into 13 clusters
where R has 12, at R's k = 12 the two partitions disagree on 417 of 4950 cell
pairs, and the sorted merge heights differ by up to 6.56e-3. Over 2000 random
cell orders SciPy gives 9 to 14 clusters (12 in 1032 orders), 0 to 764 pairs
disagree with R, and the sorted heights differ from R's by 2.8e-17 to 1.45e-2.
GitHub's runners (run 44194178528; SciPy 1.15.3 under Python 3.10 and 1.17.1
under 3.11) give the same stored-order figures on Linux x86-64 and macOS arm64.
Over the 2000 orders Linux gives 12 clusters in 1036 orders where macOS gives
1032; the range 9 to 14, the pair range and the height range are the same.
A change of 1 ULP in the input has the same effect: with the default pandas
parser the same code gave k = 12 with 21 disagreeing pairs on the M4 and
k = 11 on the Linux CI runner. Earlier versions of this report gave those
default-parser figures (k = 12, 21 pairs, 3.23e-3).
`python tests/parity/report_edc_tie_orders.py` prints these figures. Continuous
integration runs it on Linux and macOS.

## Mini CMIP6 pipeline (`tests/test_cmip6_rparity.py`)

These tests use the mini CMIP6 fixture and are listed here for completeness.

- `test_local_mle_matches_r` applies the likelihood criterion at each of 25
  cells with the package default ensemble. At ensemble 5 it fails: 4 of 25
  cells are worse than R, by up to 0.79 NLL units. The start budget is set by
  U3 of `docs/remove_r_runtime_plan.md`. Measured with temporary copies of the
  test: ensemble 10 fails (2 cells, 0.26), 20 passes, 40 fails (1 cell, 0.32).
  The failures are not monotone in the ensemble size, because Latin hypercube
  designs of different sizes are not nested.
- `test_ellipse_dissimilarity_matches_r_up_to_r_indexing_bug`. When both
  ellipses of a pair are empty, R's `calc_distance_ellipses` writes
  `dist_matrix[j] = 1`. That is linear index j, element (j, 1), instead of
  element (i, j). After `100 * (D + t(D))` the entry (i, j) is 0 and entry
  (1, j) gains 100. For any non-empty pair, 1 - IoU < 1, so a Python value of
  exactly 100 occurs only when both ellipses are empty. The test applies R's
  bug to Python's matrix and compares it with R's to rtol 5e-15. It also
  asserts that the fixture contains at least one such pair. On this fixture 44
  entries differ by 100 because of the bug. After replay the matrices agree to
  3.1e-13. Entries the bug does not touch agree to 3.55e-14.
- `test_lec_on_r_smoothed_estimates_matches_r_q30_and_k`. On R's smoothed
  estimates and the same frechet margins, LEC gives q30 = 80 and k = 11, as R
  does. The partition is not asserted. Because of the indexing bug above,
  Python's ellipse matrix differs from R's, and Python's partition differs
  from R's in 57 cell pairs. SciPy run on R's own matrix reproduces R's
  partition.

---

## Notes

1. **EC to covariance.** Python uses `scipy.optimize.brentq` (xtol 1e-15). R
   uses `uniroot` with its default tolerance, about 1.22e-4. The test checks
   that `cov_to_ec(df, root)` reproduces the target EC. It does not compare the
   roots. Measured: Python max abs(cov_to_ec(df, root) - target) is 4.44e-16.
   R's is 1.12e-4. The roots differ by at most 2.59e-5.
2. **Mutation checks.** The likelihood test fails when the ensemble is 1. The
   main-fixture LEC-on-R-estimates test fails when smoothed estimate row 4 (0-based),
   column a, is shifted by +0.5. That changes q10 from 36.893 to 36.923. Not
   every +0.5 shift changes q10 or k. Measured over 300 single-entry shifts on
   the main fixture, 204 change q10 or k, and 241 change q10, k or the
   partition.

## How to reproduce

```bash
# Run the parity tests with R hidden
PATH=/usr/bin:/bin .venv/bin/python -m pytest tests/parity/ -v -s

# Regenerate the frozen reference data (requires R 4.5.2; maintainers only)
cd tests && Rscript generate_r_reference.R
```

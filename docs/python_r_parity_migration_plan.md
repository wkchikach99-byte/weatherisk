# Python R-Parity Migration Plan

## Goal

Rebuild the scientifically critical methodology from `r_code/functions.R` as an exact Python conformance layer and prove function-level parity against frozen R outputs.

The acceptance standard for this plan is stricter than the legacy Python/R cross-validation tests:

- every deterministic function that is declared passing must match raw R outputs with absolute differences strictly below `1e-14`,
- any test that accepts larger tolerances does not count as a pass under this plan,
- any comparison against a corrected-R output instead of the historical R output does not count as a pass under this plan,
- optimizer wrappers do not pass unless they return the same selected outputs as R for the same frozen inputs.

## Why Change Direction

- The repository already contains substantial Python implementations of the `functions.R` surface.
- Exact parity is now a methodological requirement, not a performance goal.
- The current Python tests mix exact checks with looser numerical tolerances and corrected-R behavior, so a dedicated Python parity plan is needed to separate strict passes from legacy partial coverage.

## Non-Negotiable Rules

- `r_code/functions.R` is the oracle for intended historical behavior unless a behavior is explicitly classified as a bug and intentionally preserved or replaced.
- A function is only passing if the current Python implementation matches raw R outputs under the `< 1e-14` contract or, for integer/index outputs, the exact normalized output contract.
- Optimizer parity means reproducing the selected optimizer result from R, not merely achieving a similar likelihood.
- Tests that check “close”, “similar partition”, “better than R”, or “corrected R” are informative but do not count toward strict acceptance.
- Known R bugs must be documented explicitly as either:
  - preserved for historical parity, or
  - corrected intentionally in a separate non-parity mode.

## Scope

### In Scope

- Python implementations corresponding to the parity-critical `functions.R` surface.
- Frozen R fixtures under `tests/reference_data/`.
- Strict Python-vs-R conformance tests for each accepted function.
- Function-level documentation of pass/fail status under the `< 1e-14` rule.

### Out of Scope

- Rust bindings and Rust fixture tests.
- “Scientifically equivalent” acceptance without exact historical output matching.
- Treating corrected-R references as interchangeable with raw R outputs.

## Target Architecture

```text
R reference code
  -> frozen fixture generator
  -> Python conformance implementation
  -> Python parity tests and pipeline routing
```

For the parity-critical path, Python must behave as an exact conformance implementation of `functions.R`.

## Exact R Functions To Match

The migration target is specified at the R-function level. Each parity-critical R function must have one Python equivalent and one strict fixture-based conformance check proving that the Python function produces the same outputs under the agreed parity mode.

Each migration step should remain atomic:

1. implement one function or one tightly coupled helper group,
2. add or update raw R fixtures,
3. make all relevant Python tests pass under the strict `< 1e-14` contract,
4. update documentation for that function step,
5. keep the repository in a green state after the step.

### Tier 1: Mandatory parity-critical functions

These define the scientific path used by the current clustering workflow and must be accepted first.

1. `cov_fkt_2d`
  passes: true
2. `pairwise_density_summand`
  passes: true
3. `pairwise_density_optim_local`
  passes: false
4. `smooth_local_estimates`
  passes: true
5. `calc_distance_ellipses`
  passes: false
6. `clustering`
  passes: false
7. `cluster_number_threshold_method`
  passes: false
8. `llh_in_cluster`
  passes: true
9. `calc_estimates_in_clusters`
  passes: false
10. `c_extrcoeff_matrix`
  passes: true

### Tier 2: Required helper functions for exact behavior

These helpers control indexing, array construction, and deterministic transforms and therefore also need exact matching if used by Tier 1 functions.

1. `dist_x`
  passes: true
2. `dist_y`
  passes: true
3. `grid_number`
  passes: true
4. `number_grid`
  passes: true
5. `koord_num`
  passes: true
6. `crop_matrix`
  passes: false
7. `crop_local_estimates`
  passes: false
8. `cov_to_ec`
  passes: true
9. `ec_to_cov`
  passes: false

### Tier 3: Optional or secondary for the first strict slice

These exist in `functions.R` but do not yet control the current strict parity acceptance matrix.

1. `pairwise_density_optim`
2. `cov_fkt_2d_nonstat2`
3. `c_real_extrcoeff_matrix`
4. `number_koord`
5. plotting and visualization helpers
6. simulation helpers

## Python Equivalence Contract

For each Tier 1 and Tier 2 R function, define:

1. one Python function or method with a stable, documented signature,
2. raw R-generated fixtures covering nominal and edge cases,
3. one strict Python test that compares Python output directly against the frozen R fixture,
4. pipeline-routing tests if the function is used by higher-level code,
5. documentation stating whether the function is historical parity, corrected behavior, or dual-mode.

The default rule is exact output parity under absolute tolerance `< 1e-14`. If exact parity is impossible because the function is an optimizer wrapper, then the full selected-result contract must still match exactly: starting design, retries, branch selection, and returned parameters.

Authoritative conformance checks for this plan live in Python because the target implementation is Python. A function is not accepted if its only evidence is:

1. a loose-tolerance test,
2. an end-to-end cluster count or likelihood comparison,
3. a comparison against corrected-R output,
4. a claim that Python is “as good as” or “better than” R.

## Per-Function Acceptance Criteria

The following acceptance criteria apply to each function before its `passes` flag may be changed to `true`.

Status booleans in the tracking matrix below are methodology-scoped:

- `true` means the criterion is satisfied under the strict Python methodology in this file and explicitly verified.
- `false` means it is not yet satisfied, or it only has legacy/partial coverage that does not count for strict acceptance.

Use the short IDs `A1` through `A8` consistently when updating the matrix.

1. `A1` A Python equivalent exists and is reachable from the Python API.
2. `A2` Raw R fixtures exist for nominal and edge cases relevant to that function.
3. `A3` Strict Python conformance tests pass against raw R fixtures with absolute differences `< 1e-14`, or exact normalized integer/index equality where applicable.
4. `A4` The pipeline or higher-level Python code path calls the parity implementation where that step requires it.
5. `A5` No looser acceptance path is required for the function to pass under this plan.
6. `A6` Documentation is updated for the migrated function and any historical-parity caveats.
7. `A7` The relevant strict test subset is green before and after the step.
8. `A8` The function has no unresolved parity caveat that changes the raw R output contract.

`overall_passes: true` is allowed only if all of `A1` through `A8` are `true`.

## Acceptance Tracking Matrix

Legend:

- `A1` Python equivalent exists
- `A2` Raw R fixtures exist
- `A3` Strict Python-vs-R tests pass under `< 1e-14`
- `A4` Python path routed to parity implementation
- `A5` No looser fallback required
- `A6` Documentation updated here
- `A7` Relevant strict subset green
- `A8` No unresolved raw-R parity caveat

| Function | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | overall_passes | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| `cov_fkt_2d` | true | true | true | true | true | true | true | true | true | Observed max abs diff against `cov_fkt_2d_test_cases.csv`: `1.44e-15`. |
| `pairwise_density_summand` | true | true | true | true | true | true | true | true | true | Observed max abs diff against `pairwise_density_test_cases.csv`: `5.33e-15`. |
| `pairwise_density_optim_local` | true | true | false | true | false | true | false | false | false | Selected-point optimizer outputs still differ substantially from R. Observed max abs diff on `local_estimates_selected_points.csv`: `7.24e-01`; exact returned parameters do not match. |
| `smooth_local_estimates` | true | true | true | true | true | true | true | true | true | Observed max abs diffs against `local_estimates_smoothed.csv`: `a=6.22e-15`, `b=7.99e-15`, `g=4.22e-15`. |
| `calc_distance_ellipses` | true | true | false | true | false | true | false | false | false | Raw R parity is not yet met. Observed max abs diff versus raw `ellipse_dissimilarity_matrix.csv`: `5.68e-14`. Existing tests also rely on corrected-R behavior in `tests/test_cmip6_rparity.py`, which does not count for strict historical parity. |
| `clustering` | true | true | false | true | false | true | false | false | false | LEC linkage heights differ from raw R by `6.39e-14`; Saunders heights differ by `1.04e-02`. Existing tests allow looser tolerances and approximate partition agreement. |
| `cluster_number_threshold_method` | true | true | false | true | false | true | false | false | false | LEC cluster count matches raw R on the frozen matrix, but Saunders count differs by `+1` on the raw madogram matrix. Function cannot be accepted while the clustering contract still relies on approximate acceptance. |
| `llh_in_cluster` | true | true | true | true | true | true | true | true | true | Observed max abs diff against `llh_in_clusters_lec.csv`: `5.33e-15`. |
| `calc_estimates_in_clusters` | true | false | false | true | false | true | false | false | false | Python implementation exists, but there are no raw R fixtures or strict conformance tests for returned cluster estimates. |
| `c_extrcoeff_matrix` | true | true | true | true | true | true | true | true | true | Observed max abs diff against raw R fixtures: madogram `4.72e-16`, extremal coefficient `4.99e-16`. |
| `dist_x` | true | true | true | true | true | true | true | true | true | Exact arithmetic helper. Observed max abs diff: `0`. |
| `dist_y` | true | true | true | true | true | true | true | true | true | Exact arithmetic helper. Observed max abs diff: `0`. |
| `grid_number` | true | true | true | true | true | true | true | true | true | Integer mapping matches R exactly after the documented 0-based normalization. |
| `number_grid` | true | true | true | true | true | true | true | true | true | Integer inverse mapping matches R exactly after the documented 0-based normalization. |
| `koord_num` | true | true | true | true | true | true | true | true | true | Coordinate-equivalent fixture confirms the same selected grid point as R under normalized indexing conventions. |
| `crop_matrix` | false | false | false | false | false | true | false | false | false | No Python equivalent was found in the current codebase. |
| `crop_local_estimates` | false | false | false | false | false | true | false | false | false | No Python equivalent was found in the current codebase. |
| `cov_to_ec` | true | true | true | true | true | true | true | true | true | Observed max abs diff against `cov_to_ec_test_cases.csv`: `5.11e-15`. |
| `ec_to_cov` | true | true | false | true | false | true | false | false | false | Observed max abs diff against raw R fixtures: `2.59e-05`. Current tests explicitly allow `1e-4`, so legacy coverage does not count for strict parity. |

## Function Mapping Table

| R function | Required Python equivalent | Exactness requirement | Current status |
|---|---|---|---|
| `dist_x`, `dist_y` | `weatherisk.grid.dist_x`, `weatherisk.grid.dist_y` | exact | accepted |
| `grid_number`, `number_grid`, `koord_num` | `weatherisk.grid.Grid` methods | exact after normalized indexing convention | accepted |
| `cov_fkt_2d` | `weatherisk.covariance.cov_fkt_2d` | exact to `< 1e-14` | accepted |
| `pairwise_density_summand` | `weatherisk.density.pairwise_density_summand` | exact to `< 1e-14` | accepted |
| `pairwise_density_optim_local` | `weatherisk.density.pairwise_density_optim_local` | exact selected output | failing |
| `smooth_local_estimates` | `weatherisk.estimation.smooth_local_estimates` | exact to `< 1e-14` | accepted |
| `calc_distance_ellipses` | `weatherisk.clustering.calc_distance_ellipses` | exact raw-R output | failing |
| `clustering` | `weatherisk.clustering.clustering` | exact linkage output | failing |
| `cluster_number_threshold_method` | `weatherisk.clustering.cluster_number_threshold_method` | exact raw-R cluster count | failing |
| `llh_in_cluster` | `weatherisk.estimation.llh_in_cluster` | exact to `< 1e-14` | accepted |
| `calc_estimates_in_clusters` | `weatherisk.estimation.calc_estimates_in_clusters` | exact selected output | unverified |
| `c_extrcoeff_matrix` | `weatherisk.clustering.c_extrcoeff_matrix` | exact to `< 1e-14` | accepted |
| `crop_matrix`, `crop_local_estimates` | Python equivalents required | exact | missing |
| `cov_to_ec`, `ec_to_cov` | `weatherisk.covariance` helpers | exact | `cov_to_ec` accepted, `ec_to_cov` failing |

## Migration Units

Port and validate these Python function groups separately, in this order:

1. `dist_x`, `dist_y`, `grid_number`, `number_grid`, `koord_num`.
2. `cov_fkt_2d`, `cov_to_ec`, `ec_to_cov`, `pairwise_density_summand`.
3. `c_extrcoeff_matrix`.
4. `pairwise_density_optim_local` including exact start generation and optimizer contract.
5. `smooth_local_estimates`.
6. `calc_distance_ellipses` with an explicit decision on whether strict parity preserves the historical R bug.
7. `clustering` and `cluster_number_threshold_method`.
8. `llh_in_cluster`.
9. `crop_matrix` and `crop_local_estimates`.
10. `calc_estimates_in_clusters`.

This order isolates the current failure points while preserving the exact acceptance methodology.

## Phase Checklist

### Phase 1: Freeze the strict specification and fixtures
passes: false

Required in this phase:

1. every accepted function must have raw R fixtures,
2. the strict acceptance threshold must be `< 1e-14`,
3. corrected-R comparisons must be explicitly labeled non-accepting,
4. optimizer fixtures must include exact starts and exact returned outputs.

### Phase 2: Accept deterministic helpers and kernels
passes: true

Accepted in this phase:

1. `dist_x`
2. `dist_y`
3. `grid_number`
4. `number_grid`
5. `koord_num`
6. `cov_fkt_2d`
7. `cov_to_ec`
8. `pairwise_density_summand`
9. `c_extrcoeff_matrix`
10. `smooth_local_estimates`
11. `llh_in_cluster`

Still failing in this phase:

1. `ec_to_cov`

### Phase 3: Accept optimizer contracts
passes: false

Current state:

1. raw selected-point fixtures exist for `pairwise_density_optim_local`,
2. current Python outputs still differ materially from R,
3. existing “better than R” and loose-tolerance optimizer tests do not count under this plan.

### Phase 4: Accept distance and clustering chain
passes: false

Current state:

1. `calc_distance_ellipses` is close but not inside the `< 1e-14` contract against raw R,
2. existing corrected-R tests prove a different contract than historical parity,
3. LEC and Saunders clustering still rely on approximate acceptance rather than exact raw-R linkage output.

### Phase 5: Accept missing helper and cluster re-estimation functions
passes: false

Current state:

1. `crop_matrix` and `crop_local_estimates` are not implemented in Python,
2. `calc_estimates_in_clusters` has no raw R fixture set and no strict conformance test.

## Verified Evidence Snapshot

The following observed maximum absolute differences were computed directly against frozen raw R fixtures in the current workspace.

| Function | Observed max abs diff |
|---|---|
| `dist_x` | `0` |
| `dist_y` | `0` |
| `cov_fkt_2d` | `1.44e-15` |
| `cov_fkt_2d_nonstat2` | `5.55e-16` |
| `cov_to_ec` | `5.11e-15` |
| `ec_to_cov` | `2.59e-05` |
| `dtdiff` | `2.50e-16` |
| `pairwise_density_summand` | `5.33e-15` |
| `smooth_local_estimates(a)` | `6.22e-15` |
| `smooth_local_estimates(b)` | `7.99e-15` |
| `smooth_local_estimates(g)` | `4.22e-15` |
| `c_extrcoeff_matrix` madogram mode | `4.72e-16` |
| `c_extrcoeff_matrix` EC mode | `4.99e-16` |
| `calc_distance_ellipses` vs raw R | `5.68e-14` |
| `clustering` LEC merge heights vs raw R | `6.39e-14` |
| `clustering` Saunders merge heights vs raw R | `1.04e-02` |
| `llh_in_cluster` | `5.33e-15` |
| `pairwise_density_optim_local` selected outputs | `7.24e-01` |

## Immediate Next Steps

1. Decide whether strict Python parity should preserve the historical `calc_distance_ellipses` linear-indexing bug by default, or expose a separate historical-parity mode and a corrected-science mode.
2. Replace `ec_to_cov` with an implementation that reproduces R’s inverse mapping to `< 1e-14` on the frozen fixture set.
3. Reproduce R’s exact `pairwise_density_optim_local` contract, including start generation and optimizer branch behavior, until selected-point outputs match exactly.
4. Add Python implementations and raw-R fixtures for `crop_matrix`, `crop_local_estimates`, and `calc_estimates_in_clusters`.
5. Tighten the existing Python parity tests so that any acceptance above `< 1e-14` is marked legacy and kept separate from the strict parity gate.
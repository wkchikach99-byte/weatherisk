# Rust R-Parity Migration Plan

## Scope Note

This document defines a Rust migration strategy.
It is not the authority for strict constituent-function Python-to-R parity.

For the current Python strict-acceptance work, the only authoritative document
is:

- `docs/python_r_parity_migration_plan.md`

Any stronger acceptance claim for Python constituent functions must be tracked
there, not here.

## Goal

Rebuild the scientifically critical R methodology in Rust as the authoritative implementation layer, expose it through Python bindings, and prove conformance against R with fixture-based tests at each methodological step.

## Why Change Direction

- The current hybrid path is fast in parts, but exact R parity still fails in the local-MLE and LEC chain.
- Kernel-level math already matches to machine precision; the remaining drift is in algorithmic behavior, especially optimizer starts, retries, and branch selection.
- Treating R as the specification and Rust as the conformance implementation is clearer scientifically than mixing Python reimplementation and Rust acceleration.

## Non-Negotiable Rules

- R is the oracle for intended methodology unless a behavior is confirmed to be a bug.
- Deterministic functions must match R exactly up to fixture precision.
- Optimizer parity means reproducing the full optimizer contract, not just the objective.
- Python becomes orchestration and bindings, not an alternative scientific implementation for the critical Figure 9 path.
- Known R bugs must be documented explicitly as either:
  - preserved temporarily for strict historical parity, or
  - corrected intentionally with a written justification.

## Scope

### In Scope

- CMIP6 Figure 9 parity-critical path.
- Rust implementations for the R methodological steps.
- R-generated fixtures and Rust conformance tests.
- Python bindings to call Rust implementations.

### Out of Scope

- Blind line-by-line port of all R code.
- Reproducing known R bugs silently.
- Premature optimization outside the parity-critical path.

## Target Architecture

```text
R reference code
  -> frozen fixture generator
  -> Rust conformance implementation
  -> Python bindings / pipeline orchestration
```

For the Figure 9 path, the scientific implementation should live in Rust once parity is proven.

## Exact R Functions To Port

The migration target should be specified at the R-function level. Each parity-critical R function must have one Rust equivalent and one fixture-based conformance test proving the Rust function produces the same outputs under the agreed parity mode.

Each function migration step must be atomic:

1. implement one function or one tightly coupled helper group,
2. add or update fixtures,
3. make all relevant Rust and Python tests pass,
4. update documentation for that function step,
5. create one clean commit for that completed step.

Do not batch multiple partially finished functions into one step. The repository should remain in a green state after each atomic addition.

### Tier 1: Mandatory parity-critical functions

These define the Figure 9 scientific path and must be ported first.

1. `cov_fkt_2d`
  passes: false
2. `pairwise_density_summand`
  passes: false
3. `pairwise_density_optim_local`
  passes: false
4. `smooth_local_estimates`
  passes: false
5. `calc_distance_ellipses`
  passes: false
6. `clustering`
  passes: false
7. `cluster_number_threshold_method`
  passes: false
8. `llh_in_cluster`
  passes: false
9. `calc_estimates_in_clusters`
  passes: false

### Tier 2: Required helper functions for exact behavior

These helpers control indexing and array construction and therefore also need exact matching if used by Tier 1 functions.

1. `dist_x`
  passes: false
2. `dist_y`
  passes: false
3. `grid_number`
  passes: false
4. `number_grid`
  passes: false
5. `koord_num`
  passes: false
6. `crop_matrix`
  passes: false
7. `crop_local_estimates`
  passes: false
8. `cov_to_ec`
  passes: false
9. `ec_to_cov`
  passes: false

### Tier 3: Optional or non-critical for the first migration slice

These are not the current source of Figure 9 drift and can wait.

1. `pairwise_density_optim`
2. `cov_fkt_2d_nonstat2`
3. `c_extrcoeff_matrix`
4. `c_real_extrcoeff_matrix`
5. plotting and visualization helpers
6. simulation helpers

## Rust Equivalence Contract

For each Tier 1 and Tier 2 R function, define:

1. one Rust function with a stable, documented signature,
2. one R fixture generator covering nominal and edge cases,
3. one Rust test that compares Rust output against the R fixture,
4. one Python binding test if the function is exposed through PyO3.

The default rule is exact output parity. If exact parity is impossible because the function is an optimizer wrapper, then the full decision contract must still match exactly: starts, retries, branch selection, and selected output.

Authoritative equivalence tests belong in Rust, not only in Python. The canonical conformance check for a migrated function is:

1. frozen R outputs stored under `tests/reference_data/`,
2. Rust fixture tests in the crate reading those frozen R outputs directly,
3. optional Python tests for bindings, routing, and end-to-end orchestration.

Python tests are still required where bindings or pipeline routing are part of the migration step, but they are secondary evidence. A function is not scientifically accepted if its only R-equivalence test lives in Python.

## Per-Function Acceptance Criteria

The following acceptance criteria apply to each function before its `passes` flag may be changed to `true`.

Status booleans in the tracking matrix below are methodology-scoped:

- `true` means the criterion is satisfied under the current migration methodology and explicitly verified.
- `false` means it is not yet satisfied, or it has only legacy/partial coverage that does not yet count for strict acceptance.

Use the short IDs `A1` through `A8` consistently when updating the matrix.

1. `A1` A Rust equivalent exists and is reachable from the crate API.
2. `A2` R-generated fixtures exist for nominal and edge cases relevant to that function.
3. `A3` Rust conformance tests pass against the R fixtures.
4. `A4` Python binding tests pass if the function is exposed through PyO3.
5. `A5` The Python code path calls the Rust implementation where that migration step requires it.
6. `A6` Documentation is updated for the migrated function and any parity-mode caveats.
7. `A7` The relevant test subset is green before and after the step.
8. `A8` The step is committed as one atomic change.

`overall_passes: true` is allowed only if all of `A1` through `A8` are `true`.

## Acceptance Tracking Matrix

Legend:

- `A1` Rust equivalent exists
- `A2` R fixtures exist
- `A3` Rust tests from frozen R outputs pass
- `A4` Python binding tests pass
- `A5` Python path routed to Rust
- `A6` Documentation updated
- `A7` Relevant subset green
- `A8` Atomic step committed

| Function | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | overall_passes | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| `cov_fkt_2d` | true | true | true | true | true | true | true | true | true | Accepted in unit 2: Rust fixture tests are authoritative; Python tests cover binding and routed covariance calls. |
| `pairwise_density_summand` | true | true | true | true | true | true | true | true | true | Accepted in unit 2: Rust fixture tests are authoritative; Python tests cover binding and routed density calls. |
| `pairwise_density_optim_local` | true | true | false | false | false | true | false | false | false | R fixtures are now frozen, including the exact `maximinLHS(10,3)` start matrix and selected-point outputs, but exact selected-result parity is still unresolved in both the Rust binding and the current Python translation. |
| `smooth_local_estimates` | false | false | false | false | false | false | false | false | false | Not yet migrated into Rust. |
| `calc_distance_ellipses` | true | true | false | false | false | false | false | false | false | Rust implementation exists, but parity mode and strict Rust-vs-R fixture acceptance are not finalized. |
| `clustering` | false | false | false | false | false | false | false | false | false | Rust equivalent not yet implemented. |
| `cluster_number_threshold_method` | false | false | false | false | false | false | false | false | false | Rust equivalent not yet implemented. |
| `llh_in_cluster` | false | false | false | false | false | false | false | false | false | Rust equivalent not yet implemented. |
| `calc_estimates_in_clusters` | false | false | false | false | false | false | false | false | false | Rust equivalent not yet implemented. |
| `dist_x` | true | true | true | true | true | true | true | true | true | Helper unit accepted: Rust fixture tests are authoritative; Python tests cover binding and routing. |
| `dist_y` | true | true | true | true | true | true | true | true | true | Helper unit accepted: Rust fixture tests are authoritative; Python tests cover binding and routing. |
| `grid_number` | true | true | true | true | true | true | true | true | true | Helper unit accepted: Rust fixture tests are authoritative; Python tests cover binding and routing. |
| `number_grid` | true | true | true | true | true | true | true | true | true | Helper unit accepted: Rust fixture tests are authoritative; Python tests cover binding and routing. |
| `koord_num` | true | true | true | true | true | true | true | true | true | Helper unit accepted: Rust fixture tests are authoritative; Python tests cover binding and routing. |
| `crop_matrix` | false | false | false | false | false | false | false | false | false | Not yet migrated into Rust. |
| `crop_local_estimates` | false | false | false | false | false | false | false | false | false | Not yet migrated into Rust. |
| `cov_to_ec` | true | true | true | true | true | true | true | true | true | Accepted in unit 2: Rust implementation added and routed into Python covariance helpers. |
| `ec_to_cov` | true | true | true | true | true | true | true | true | true | Accepted in unit 2: Rust implementation added and routed into Python covariance helpers. |

## Function Mapping Table

| R function | Required Rust equivalent | Exactness requirement | Notes |
|---|---|---|---|
| `dist_x`, `dist_y` | `dist_x`, `dist_y` | exact | trivial but needed for pair construction |
| `grid_number`, `number_grid`, `koord_num` | same logical mapping | exact | indexing errors here poison all downstream parity |
| `cov_fkt_2d` | `cov_fkt_2d` | exact to machine precision | deterministic kernel |
| `pairwise_density_summand` | `pairwise_density_summand` | exact to machine precision | already close; keep fixture tests |
| `pairwise_density_optim_local` | `pairwise_density_optim_local` | exact selected output | highest priority; current failure source |
| `smooth_local_estimates` | `smooth_local_estimates` | exact | includes angle-centering behavior |
| `calc_distance_ellipses` | `calc_distance_ellipses` | exact under chosen parity mode | known R bug must be classified explicitly |
| `clustering` | `clustering_average` | exact linkage / merge output | depends on distance input parity |
| `cluster_number_threshold_method` | same logical function | exact | trivial once linkage matches |
| `llh_in_cluster` | `llh_in_cluster` | exact | deterministic aggregator |
| `calc_estimates_in_clusters` | `calc_estimates_in_clusters` | exact selected output | second-phase optimizer port |

## Migration Units

Port and validate these function groups separately, in this order:

1. `dist_x`, `dist_y`, `grid_number`, `number_grid`, `koord_num`.
2. `cov_fkt_2d`, `cov_to_ec`, `ec_to_cov`, `pairwise_density_summand`.
3. `pairwise_density_optim_local` including start generation and optimizer contract.
4. `smooth_local_estimates`.
5. `calc_distance_ellipses`.
6. `clustering` and `cluster_number_threshold_method`.
7. `llh_in_cluster`.
8. `calc_estimates_in_clusters`.

This order isolates the current failure point before touching unrelated parts.

## Phase Checklist

### Phase 1: Freeze the specification and fixtures
passes: false

### Phase 2: Migrate helper and indexing functions
passes: false

Accepted in this phase so far:

1. frozen R fixtures exist for the helper group,
2. authoritative Rust-vs-R fixture tests now live in the Rust crate,
3. Python binding and routing checks are green,
4. the helper acceptance step has been committed cleanly.

1. `dist_x`
2. `dist_y`
3. `grid_number`
4. `number_grid`
5. `koord_num`

### Phase 3: Migrate deterministic Tier 1 kernels
passes: false

Accepted in this phase so far:

1. `cov_fkt_2d`
2. `cov_to_ec`
3. `ec_to_cov`
4. `pairwise_density_summand`

### Phase 4: Migrate `pairwise_density_optim_local`
passes: false

Current state:

- Rust local optimizer now matches more of the R control contract: `b >= 0.01`, scaled-space bounds, gamma-wrap retry, and boundary-triggered extra starts.
- Frozen R artifacts now include the exact `maximinLHS(10,3)` start matrix used by the parity-critical `set.seed(42)` local-estimation path.
- A Rust selected-point fixture check exists as a diagnostic, but remains non-accepting because exact selected-output parity has not yet been achieved.
- The same selected-point mismatch is still observable in the current Python translation, so the remaining blocker is not Rust kernel math alone.

### Phase 5: Migrate post-MLE functions
passes: false

### Phase 6: Bind the full parity-critical path into Python
passes: false

### Phase 7: End-to-end mini Figure 9 validation and HPC revalidation
passes: false

## Testing Contract

### Fixture Strategy

- Generate fixtures from R for each migration unit.
- Store both inputs and expected outputs.
- Use fixed seeds and frozen mini datasets.
- Include edge cases: boundary hits, isotropic cells, flat-likelihood cells, empty ellipses.

### Acceptance Levels

- Tier 2 helpers and deterministic Tier 1 functions: exact or machine-precision agreement on every fixture.
- `pairwise_density_optim_local` and `calc_estimates_in_clusters`: exact agreement on the selected result under the full R decision contract.
- End-to-end mini Figure 9 path: exact agreement on raw local estimates, smoothed estimates, LEC/EDC inputs, `q30`, and `k`.

### Test Layout

- Rust unit tests: authoritative function-level fixture checks against frozen R outputs.
- Python integration tests: binding correctness, routing correctness, and end-to-end mini-pipeline checks.
- R fixture generator tests: ensure fixtures are reproducible and versioned.

## Optimizer Plan

The local-MLE optimizer is the critical migration item.

Rust must reproduce the effective R contract:

1. Same neighborhood/pair construction.
2. Same parameter bounds and parameter scaling.
3. Same start-set generation semantics as `maximinLHS`.
4. Same multi-start execution order.
5. Same boundary-triggered rerun logic.
6. Same gamma-wrap retry behavior at $\pm \pi/2$.
7. Same best-solution selection rule.
8. Deterministic handling of flat-likelihood ties.

If exact R `optim` behavior cannot be reproduced directly, the fallback standard is to emulate the full observed decision contract from R fixtures, not merely to minimize the same objective.

## Known Methodological Decision Needed

At least one R behavior is already identified as a bug in ellipse dissimilarity indexing. Before full migration, classify each such case as:

1. historical-parity mode,
2. corrected-method mode, or
3. dual-mode with explicit documentation.

This decision must be made before claiming exact parity in papers.

## Execution Phases

### Phase 1: Freeze the Specification

- Enumerate the exact R functions listed in this document.
- Create fixture generators for each Tier 1 and Tier 2 function.
- Version the fixtures and document seeds, inputs, and assumptions.

### Phase 2: Rebuild `pairwise_density_optim_local` in Rust

- Implement its helper/indexing contract exactly.
- Implement R-equivalent pair construction.
- Implement R-equivalent start generation.
- Implement R-equivalent multi-start optimizer contract.
- Validate cell-by-cell against mini CMIP6 fixtures.

### Phase 3: Rebuild Post-MLE Steps in Rust

- Port smoothing.
- Port ellipse dissimilarity.
- Validate corrected-vs-historical behavior explicitly where needed.

### Phase 4: Bind into Python

- Replace parity-critical Python implementations with Rust calls.
- Keep Python as orchestration and I/O only.
- Preserve the existing public pipeline API where possible.

### Phase 5: Re-validate the Full Mini Pipeline

- Run the exact mini Figure 9 chain against R reference outputs.
- Require exact parity on raw estimates, smoothed estimates, and final clustering outputs.

### Phase 6: Performance and HPC Validation

- Benchmark the Rust-native path after parity is proven.
- Re-run on albedo and compare runtime, memory, and outputs.

## Deliverables

These are deliverables of the migration program, not pass/fail status markers.

1. Rust conformance modules for the parity-critical Figure 9 path.
2. R fixture generator scripts and frozen reference artifacts.
3. Rust unit tests for each migrated function.
4. Python integration tests for binding and mini-pipeline parity.
5. Updated validation documentation reflecting the new stricter parity standard.

## Risks

- Exact parameter parity is not guaranteed by matching only the objective; optimizer semantics must also match.
- Existing validation notes in `docs/software_x/` are now too optimistic for the tightened parity standard and must not be treated as final evidence.
- Some R behaviors may be accidental or buggy, so parity claims require explicit methodological classification.

## Recommended First Implementation Slice

Start with `pairwise_density_optim_local` and its direct dependencies only.

Success criterion for the first slice:

1. Rust reproduces the R outputs for `pairwise_density_optim_local` on the frozen mini CMIP6 fixtures.
2. The Rust path reproduces the R-selected branch in flat-likelihood cells.
3. Python calls the Rust implementation directly for local MLE in the Figure 9 path.

## Definition of Done

The migration is complete for the Figure 9 path when:

1. The parity-critical steps are implemented in Rust.
2. Rust passes fixture-based conformance tests against R at each step.
3. The Python Figure 9 pipeline uses the Rust implementations for those steps.
4. The mini CMIP6 pipeline matches R exactly under the agreed parity mode.
5. Performance validation is rerun only after scientific parity is demonstrated.
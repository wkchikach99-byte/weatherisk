# Benchmark Results

This file records benchmark runs for the CMIP6 Figure 9 pipeline.
Each run uses the same `run_cmip6_pipeline` orchestration path used by the HPC script,
with deterministic synthetic data injected only at the data-loading boundary.

The obsolete hot-path-only benchmark has been removed. The entries below track the
actual optimization sequence on the full pipeline path.

Note: only comparison runs that support an actual optimization decision are kept below.

## Decision Benchmark Protocol

- Default benchmark class from now on: `decision benchmark`
- Medium case: `48 years`, `16x16` grid, `4 workers`
- Sampling plan: `1 warmup + 5 measured runs`
- Warmup policy: warmup runs are discarded from the reported summary
- Reported statistics: total-time `mean/min/max/std`, per-step `mean/min/max/std`, and peak RSS `mean/max`
- Stability rule: benchmark summaries should also record whether the numerical check dictionary stayed identical across measured runs

## Run 2026-03-09T15:36:19.463268+00:00

- Git revision: `0bc6035`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Variant: `pre-improvement baseline`
- Total time: `5.755s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 1, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 11, 'k_edc': 10}`

| Step | Seconds |
| --- | ---: |
| _stl_detrend_grid | 2.662 |
| _monthly_annual_maxima | 0.000 |
| _compute_frechet_global | 0.459 |
| _run_local_estimation_cmip6 | 1.549 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.004 |
| _incluster_reestimate_cmip6 | 1.079 |

- Checks: `{'frechet_min': 0.17650048124730916, 'frechet_max': 190.1910742661081, 'edc_trace': 0.0, 'est_mean_a': 0.03226419358541949, 'est_mean_b': 1.4137247359882381, 'est_mean_gamma': 0.11852993719742141}`

## Run 2026-03-09T15:36:16.294966+00:00

- Git revision: `b575df2`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Variant: `post low-hanging fruits, serial`
- Total time: `2.543s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 1, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 13, 'k_edc': 10}`

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.000 |
| _monthly_annual_maxima | 0.194 |
| _compute_frechet_global | 0.437 |
| _run_local_estimation_cmip6 | 0.964 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.006 |
| _incluster_reestimate_cmip6 | 0.941 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'edc_trace': 0.0, 'est_mean_a': 0.023881369901431172, 'est_mean_b': 1.283395966261962, 'est_mean_gamma': 0.36477793020266147}`

## Run 2026-03-09T15:42:58.173978+00:00

- Git revision: `b575df2`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Variant: `pre-Step-6 cluster parallelization`
- Total time: `5.988s`
- Config: `{'seed': 12345, 'n_years': 36, 'n_lat': 12, 'n_lon': 12, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 432, 'n_cells': 144, 'n_valid_cells': 144, 'n_years_complete': 36, 'k_lec': 13, 'k_edc': 17}`

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.001 |
| _monthly_annual_maxima | 0.121 |
| _compute_frechet_global | 0.775 |
| _run_local_estimation_cmip6 | 1.313 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.014 |
| _incluster_reestimate_cmip6 | 3.761 |

- Checks: `{'frechet_min': 0.15341623656416714, 'frechet_max': 395.99565803595567, 'edc_trace': 0.0, 'est_mean_a': 0.0172096948404699, 'est_mean_b': 1.2176352663612096, 'est_mean_gamma': -0.15018490051108158}`

## Run 2026-03-09T15:42:24.101710+00:00

- Git revision: `b575df2`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Variant: `post-Step-6 cluster parallelization`
- Total time: `5.121s`
- Config: `{'seed': 12345, 'n_years': 36, 'n_lat': 12, 'n_lon': 12, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 432, 'n_cells': 144, 'n_valid_cells': 144, 'n_years_complete': 36, 'k_lec': 13, 'k_edc': 17}`

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.001 |
| _monthly_annual_maxima | 0.219 |
| _compute_frechet_global | 0.846 |
| _run_local_estimation_cmip6 | 1.356 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.016 |
| _incluster_reestimate_cmip6 | 2.680 |

- Checks: `{'frechet_min': 0.15341623656416714, 'frechet_max': 395.99565803595567, 'edc_trace': 0.0, 'est_mean_a': 0.0172096948404699, 'est_mean_b': 1.2176352663612096, 'est_mean_gamma': -0.15018490051108158}`

## Summary

### Low-Hanging Fruits

- Comparable serial benchmark: `5.755s` baseline → `2.543s` after the first three optimizations.
- Net improvement on the same full pipeline path: `2.26x` faster overall.
- Largest gain: Step 1 de-trending, `2.662s` → `0.0004s`.
- Local estimation also improved: `1.549s` → `0.964s` on the serial comparison.
- The old `4-worker` path could not be benchmarked because the pre-improvement code fails with `Can't pickle local object _stl_detrend_grid.<locals>._worker`.

### Step 6 Parallelization

- Representative `36-year, 12x12, 4-worker` benchmark: `5.988s` before Step 6 parallelization → `5.121s` after.
- Net improvement on the same full pipeline path: `1.17x` faster overall.
- In-cluster re-estimation improved from `3.761s` to `2.680s`, a `1.40x` speedup for the optimized step.

### Clustering Memory Retention

- Representative `36-year, 12x12, 4-worker` benchmark after the clustering-memory refactor: `5.400s` total, peak RSS `0.583 GiB`.
- The numerical checks remained stable (`k_LEC = 13`, `k_EDC = 17`, same estimate means at benchmark precision).
- Runtime stayed within the same band as the pre-refactor optimized pipeline, which is expected because this change targets peak memory and duplicate allocations rather than raw arithmetic cost.

### Vectorized Pair Assembly

- Representative `36-year, 12x12, 4-worker` benchmark after vectorizing pair-array setup: `5.309s` total, peak RSS `0.583 GiB`.
- The benchmark stayed numerically stable at the same check values.
- The runtime gain on this benchmark size was modest, which suggests Python pair-list construction was not yet the dominant cost relative to optimizer evaluations.

### Condensed EDC Path

- Representative `36-year, 12x12, 4-worker` benchmark after switching the no-artifact EDC path to condensed storage: `5.316s` total, peak RSS `0.579 GiB`.
- The clustering step itself became slightly cheaper while the end-to-end runtime remained essentially flat.
- This is consistent with a memory-saving refactor applied to a step that is not the main wall-clock bottleneck on the benchmark grid.

## Run 2026-03-09T19:42:19.806555+00:00

- Git revision: `b575df2`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `5.400s`
- Config: `{'seed': 12345, 'n_years': 36, 'n_lat': 12, 'n_lon': 12, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 432, 'n_cells': 144, 'n_valid_cells': 144, 'n_years_complete': 36, 'k_lec': 13, 'k_edc': 17}`
- Max memory: `0.583 GiB` (`625999872` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.001 |
| _monthly_annual_maxima | 0.213 |
| _compute_frechet_global | 0.933 |
| _run_local_estimation_cmip6 | 1.345 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.050 |
| _incluster_reestimate_cmip6 | 2.824 |

- Checks: `{'frechet_min': 0.15341623656416714, 'frechet_max': 395.99565803595567, 'labels_edc_sum': 1304, 'est_mean_a': 0.0172096948404699, 'est_mean_b': 1.2176352663612096, 'est_mean_gamma': -0.15018490051108158}`
## Run 2026-03-09T19:44:11.245125+00:00

- Git revision: `23acf77`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `5.309s`
- Config: `{'seed': 12345, 'n_years': 36, 'n_lat': 12, 'n_lon': 12, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 432, 'n_cells': 144, 'n_valid_cells': 144, 'n_years_complete': 36, 'k_lec': 13, 'k_edc': 17}`
- Max memory: `0.583 GiB` (`625721344` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.001 |
| _monthly_annual_maxima | 0.140 |
| _compute_frechet_global | 0.845 |
| _run_local_estimation_cmip6 | 1.419 |
| _smooth_estimates_cmip6 | 0.002 |
| _run_clustering_cmip6 | 0.060 |
| _incluster_reestimate_cmip6 | 2.813 |

- Checks: `{'frechet_min': 0.15341623656416714, 'frechet_max': 395.99565803595567, 'labels_edc_sum': 1304, 'est_mean_a': 0.0172096948404699, 'est_mean_b': 1.2176352663612096, 'est_mean_gamma': -0.15018490051108158}`
## Run 2026-03-09T19:45:21.083073+00:00

- Git revision: `23acf77`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `5.316s`
- Config: `{'seed': 12345, 'n_years': 36, 'n_lat': 12, 'n_lon': 12, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 432, 'n_cells': 144, 'n_valid_cells': 144, 'n_years_complete': 36, 'k_lec': 13, 'k_edc': 17}`
- Max memory: `0.579 GiB` (`622067712` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.001 |
| _monthly_annual_maxima | 0.141 |
| _compute_frechet_global | 0.872 |
| _run_local_estimation_cmip6 | 1.353 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.047 |
| _incluster_reestimate_cmip6 | 2.868 |

- Checks: `{'frechet_min': 0.15341623656416714, 'frechet_max': 395.99565803595567, 'labels_edc_sum': 1304, 'est_mean_a': 0.0172096948404699, 'est_mean_b': 1.2176352663612096, 'est_mean_gamma': -0.15018490051108158}`
## Decision Benchmark 2026-03-09T19:54:49.439181+00:00

- Git revision: `23b6ceb`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Method: `1 warmup + 5 measured runs` (warmups excluded from summary)
- Benchmark case: `medium`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 38, 'k_edc': 20}`
- Total time summary: mean `9.176s`, min `9.116s`, max `9.227s`, std `0.038s`
- Peak memory summary: mean `0.757 GiB`, max `0.772 GiB` (`peak_process_tree_rss`, Δt=0.05s)
- Checks stable across measured runs: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 0.002 | 0.001 | 0.002 | 0.000 |
| _monthly_annual_maxima | 0.001 | 0.000 | 0.001 | 0.000 |
| _compute_frechet_global | 1.022 | 0.999 | 1.066 | 0.024 |
| _run_local_estimation_cmip6 | 2.304 | 2.286 | 2.338 | 0.018 |
| _smooth_estimates_cmip6 | 0.003 | 0.002 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.076 | 0.073 | 0.080 | 0.003 |
| _incluster_reestimate_cmip6 | 5.721 | 5.681 | 5.767 | 0.035 |

- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.01500266795292907, 'est_mean_b': 0.8624371652491263, 'est_mean_gamma': 0.02906508271908983}`

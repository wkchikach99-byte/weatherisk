# Benchmark Results

This file records benchmark evidence for the CMIP6 Figure 9 pipeline.
Each run uses the real `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
orchestration path, with deterministic synthetic data injected only at the
data-loading boundary.

The benchmark policy changed on 2026-03-09. From that point onward, only the
medium repeated-run benchmark is treated as decision-grade evidence. Earlier
one-off runs are kept as historical optimization checkpoints.

## Current Decision Benchmark Baseline

- Benchmark class: `decision benchmark`
- Protocol: `1 warmup + 5 measured runs`
- Warmup policy: warmup runs are excluded from the reported summary
- Medium case: `48 years`, `16x16` grid, `4 workers`
- Reported statistics: total-time `mean/min/max/std`, per-step `mean/min/max/std`, and peak RSS `mean/max`
- Stability rule: the numerical check dictionary must stay identical across measured runs

## Decision Benchmark 2026-03-09T19:54:49.439181+00:00

- Git revision: `23b6ceb`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 38, 'k_edc': 20}`
- Total time summary: mean `9.176s`, min `9.116s`, max `9.227s`, std `0.038s`
- Peak memory summary: mean `0.757 GiB`, max `0.772 GiB` (`peak_process_tree_rss`, `Δt=0.05s`)
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

## Historical Optimization Checkpoints

The runs below were useful during implementation, but they are one-off
measurements rather than repeated-run decision benchmarks.

| Revision | Scenario | Config | Total | Peak RSS | Main observation |
| --- | --- | --- | ---: | ---: | --- |
| `0bc6035` | pre-improvement baseline | `24y, 8x8, 1 worker` | `5.755s` | n/a | serial baseline before the optimization pass |
| `b575df2` | post low-hanging fruits, serial | `24y, 8x8, 1 worker` | `2.543s` | n/a | `2.26x` faster than the serial baseline |
| `b575df2` | pre-Step-6 cluster parallelization | `36y, 12x12, 4 workers` | `5.988s` | n/a | reference point before parallel cluster re-estimation |
| `b575df2` | post-Step-6 cluster parallelization | `36y, 12x12, 4 workers` | `5.121s` | n/a | in-cluster step dropped from `3.761s` to `2.680s` |
| `b575df2` | clustering memory retention | `36y, 12x12, 4 workers` | `5.400s` | `0.583 GiB` | memory-saving refactor with stable outputs |
| `23acf77` | vectorized pair assembly | `36y, 12x12, 4 workers` | `5.309s` | `0.583 GiB` | modest runtime gain; optimizer cost still dominates |
| `23acf77` | condensed EDC path | `36y, 12x12, 4 workers` | `5.316s` | `0.579 GiB` | small memory win; total runtime essentially flat |

## Historical Notes

### Low-Hanging Fruits

- Comparable serial benchmark: `5.755s` baseline to `2.543s` after the first optimization pass.
- Net improvement on the same full pipeline path: `2.26x` faster overall.
- Largest gain: de-trending, from `2.662s` to effectively zero at this benchmark scale.
- Local estimation also improved, from `1.549s` to `0.964s` on the serial comparison.
- The old `4-worker` path could not be benchmarked before the de-trending fix because the pre-improvement code fails with `Can't pickle local object _stl_detrend_grid.<locals>._worker`.

### Step 6 Parallelization

- Representative `36-year, 12x12, 4-worker` benchmark: `5.988s` before cluster parallelization and `5.121s` after.
- Net improvement on the same full pipeline path: `1.17x` faster overall.
- In-cluster re-estimation improved from `3.761s` to `2.680s`, a `1.40x` speedup for that step.

### Memory-Oriented Refactors

- Clustering memory retention refactor: `5.400s` total, peak RSS `0.583 GiB`.
- Vectorized pair assembly: `5.309s` total, peak RSS `0.583 GiB`.
- Condensed EDC path: `5.316s` total, peak RSS `0.579 GiB`.
- These one-off runs stayed numerically stable at benchmark precision, but their runtime differences were small enough that they should not be treated as decision-grade evidence on their own.

## Run 2026-03-09T20:33:51.999721+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `3.161s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 14, 'k_edc': 10}`
- Max memory: `0.539 GiB` (`578420736` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.000 |
| _monthly_annual_maxima | 0.150 |
| _compute_frechet_global | 0.714 |
| _run_local_estimation_cmip6 | 0.757 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.034 |
| _incluster_reestimate_cmip6 | 1.476 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'labels_edc_sum': 404, 'est_mean_a': 0.023654652789610237, 'est_mean_b': 1.2909360326913053, 'est_mean_gamma': 0.36473074153383245}`
## Run 2026-03-09T20:34:50.085619+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `3.149s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 14, 'k_edc': 10}`
- Max memory: `0.538 GiB` (`577601536` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.000 |
| _monthly_annual_maxima | 0.138 |
| _compute_frechet_global | 0.700 |
| _run_local_estimation_cmip6 | 0.755 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.036 |
| _incluster_reestimate_cmip6 | 1.489 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'labels_edc_sum': 404, 'est_mean_a': 0.023654652789610237, 'est_mean_b': 1.2909360326913053, 'est_mean_gamma': 0.36473074153383245}`
## Decision Benchmark 2026-03-09T20:35:37.857473+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Method: `1 warmup + 5 measured runs` (warmups excluded from summary)
- Benchmark case: `medium`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 40, 'k_edc': 20}`
- Total time summary: mean `7.729s`, min `7.572s`, max `7.900s`, std `0.119s`
- Peak memory summary: mean `0.795 GiB`, max `0.811 GiB` (`peak_process_tree_rss`, Δt=0.05s)
- Checks stable across measured runs: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 0.002 | 0.001 | 0.002 | 0.000 |
| _monthly_annual_maxima | 0.000 | 0.000 | 0.000 | 0.000 |
| _compute_frechet_global | 0.957 | 0.922 | 1.001 | 0.031 |
| _run_local_estimation_cmip6 | 1.889 | 1.850 | 1.946 | 0.041 |
| _smooth_estimates_cmip6 | 0.003 | 0.002 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.068 | 0.065 | 0.071 | 0.002 |
| _incluster_reestimate_cmip6 | 4.771 | 4.692 | 4.852 | 0.052 |

- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.014427466088444751, 'est_mean_b': 0.8624030479245355, 'est_mean_gamma': 0.029206774817643926}`
## Run 2026-03-09T20:37:22.771154+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `3.211s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 14, 'k_edc': 10}`
- Max memory: `0.539 GiB` (`578568192` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.000 |
| _monthly_annual_maxima | 0.137 |
| _compute_frechet_global | 0.752 |
| _run_local_estimation_cmip6 | 0.752 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.040 |
| _incluster_reestimate_cmip6 | 1.502 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'labels_edc_sum': 404, 'est_mean_a': 0.023654652789610237, 'est_mean_b': 1.2909360326913053, 'est_mean_gamma': 0.36473074153383245}`
## Decision Benchmark 2026-03-09T20:38:12.243390+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Method: `1 warmup + 5 measured runs` (warmups excluded from summary)
- Benchmark case: `medium`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 40, 'k_edc': 20}`
- Total time summary: mean `8.119s`, min `8.060s`, max `8.180s`, std `0.043s`
- Peak memory summary: mean `0.794 GiB`, max `0.797 GiB` (`peak_process_tree_rss`, Δt=0.05s)
- Checks stable across measured runs: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 0.002 | 0.002 | 0.002 | 0.000 |
| _monthly_annual_maxima | 0.001 | 0.000 | 0.001 | 0.000 |
| _compute_frechet_global | 1.046 | 1.012 | 1.098 | 0.037 |
| _run_local_estimation_cmip6 | 1.971 | 1.966 | 1.983 | 0.006 |
| _smooth_estimates_cmip6 | 0.003 | 0.003 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.077 | 0.075 | 0.078 | 0.001 |
| _incluster_reestimate_cmip6 | 4.976 | 4.942 | 5.028 | 0.031 |

- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.014427466088444751, 'est_mean_b': 0.8624030479245355, 'est_mean_gamma': 0.029206774817643926}`
## Run 2026-03-09T20:44:54.827979+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `3.244s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24, 'k_lec': 13, 'k_edc': 10}`
- Max memory: `0.529 GiB` (`568426496` bytes; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.000 |
| _monthly_annual_maxima | 0.139 |
| _compute_frechet_global | 0.709 |
| _run_local_estimation_cmip6 | 0.807 |
| _smooth_estimates_cmip6 | 0.001 |
| _run_clustering_cmip6 | 0.032 |
| _incluster_reestimate_cmip6 | 1.528 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'labels_edc_sum': 404, 'est_mean_a': 0.023881369901431172, 'est_mean_b': 1.283395966261962, 'est_mean_gamma': 0.36477793020266147}`
## Decision Benchmark 2026-03-09T20:45:48.639764+00:00

- Git revision: `aa3b738`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Method: `1 warmup + 5 measured runs` (warmups excluded from summary)
- Benchmark case: `medium`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 38, 'k_edc': 20}`
- Total time summary: mean `8.828s`, min `8.780s`, max `8.892s`, std `0.043s`
- Peak memory summary: mean `0.784 GiB`, max `0.801 GiB` (`peak_process_tree_rss`, Δt=0.05s)
- Checks stable across measured runs: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 0.001 | 0.001 | 0.001 | 0.000 |
| _monthly_annual_maxima | 0.000 | 0.000 | 0.001 | 0.000 |
| _compute_frechet_global | 0.964 | 0.957 | 0.975 | 0.006 |
| _run_local_estimation_cmip6 | 2.219 | 2.212 | 2.232 | 0.007 |
| _smooth_estimates_cmip6 | 0.003 | 0.003 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.068 | 0.067 | 0.070 | 0.001 |
| _incluster_reestimate_cmip6 | 5.533 | 5.488 | 5.605 | 0.049 |

- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.01500266795292907, 'est_mean_b': 0.8624371652491263, 'est_mean_gamma': 0.02906508271908983}`

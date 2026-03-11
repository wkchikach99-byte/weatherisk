# Benchmark Results

This file is the canonical benchmark log for the CMIP6 Figure 9 path.
Each benchmark runs the real `scripts/reproduce_fig9.py` entrypoint used by the
SLURM Figure 9 job, with deterministic synthetic data injected only at the
data-loading boundary so the reduced case finishes quickly while preserving the
production call chain.

## Canonical Workflow

- Canonical runner: `python scripts/run_fig9_benchmark.py`
- Script entrypoint under test: `scripts.reproduce_fig9.main`
- Pipeline entrypoint under test: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Default intent: reduced Figure 9 case that preserves the production call chain while finishing in roughly a minute or less
- Recorded metrics: total wall time, per-step timings, peak process-tree RSS, and generated figure count
- Memory is reported in `bytes`, `KiB`, `MiB`, and `GiB`

## Notes

- A real script entrypoint is required for multiprocessing benchmarks; heredoc or stdin entrypoints can give misleading results or fail under macOS spawn mode.
- For the CMIP6 Figure 9 path, Step 3 local MLE remains the primary bottleneck; optimizer cost dominates pair-array assembly once the production likelihood is active.
- This file now tracks only the active Python-only Figure 9 benchmark workflow.

| **Peak RSS** | **0.555 GiB** | **0.596 GiB** | **1.07×** | ~4× if O(n²) |
| Total time | 11.790s | 43.974s | 3.73× | 4× |
| _detrend_grid_fast | 6.061s | 24.952s | 4.12× | 4× (O(n)) |
| _compute_frechet_global | 1.001s | 3.523s | 3.52× | 4× (O(n)) |
| _run_local_estimation_cmip6 | 2.570s | 10.546s | 4.10× | 4× (O(n)) |
| _run_clustering_cmip6 | 0.044s | 0.098s | 2.23× | 4–16× (O(n²)) |
| _incluster_reestimate_cmip6 | 2.068s | 4.775s | 2.31× | varies |

#### Scaling analysis

- **Memory scales sub-linearly**: 4× cells → only 1.07× memory. The
  condensed form and block-based NLL prevent quadratic memory growth.
  The dominant memory consumer is now the Python interpreter and worker
  processes (~0.5 GiB baseline), not the data arrays.
- **Compute scales linearly**: Steps 1–4 (per-cell work) show ~4× as
  expected. Clustering (Step 5) shows ~2.2× because the condensed
  vector is O(n²) in size but the Rust backend streams pair-by-pair
  without materialising large intermediates. In-cluster re-estimation
  (Step 6) shows 2.3× rather than 4× because larger clusters hit the
  `max_dist` filter, reducing effective pair count.
- **Extrapolation to 192×96 (18,432 cells)**: With sub-linear memory
  scaling and the 0.6 GiB baseline, the full grid should peak at
  approximately 1.5–2.0 GiB — well within the 8 GiB target.

#### 16×16 detail (Rust backend, post-optimisation)

- Total: mean `11.790s`, min `11.513s`, max `12.068s`, std `0.189s`
- Peak RSS: mean `0.555 GiB`, max `0.559 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 6.061 | 5.836 | 6.216 | 0.137 |
| _monthly_annual_maxima | 0.002 | 0.000 | 0.005 | 0.002 |
| _compute_frechet_global | 1.001 | 0.975 | 1.091 | 0.045 |
| _run_local_estimation_cmip6 | 2.570 | 2.564 | 2.588 | 0.009 |
| _smooth_estimates_cmip6 | 0.003 | 0.002 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.044 | 0.037 | 0.051 | 0.004 |
| _incluster_reestimate_cmip6 | 2.068 | 2.035 | 2.100 | 0.023 |

#### 32×32 detail (Rust backend, post-optimisation)

- Total: mean `43.974s`, min `43.381s`, max `44.720s`, std `0.430s`
- Peak RSS: mean `0.596 GiB`, max `0.605 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 24.952 | 24.479 | 25.320 | 0.275 |
| _monthly_annual_maxima | 0.001 | 0.001 | 0.001 | 0.000 |
| _compute_frechet_global | 3.523 | 3.408 | 3.836 | 0.159 |
| _run_local_estimation_cmip6 | 10.546 | 10.486 | 10.616 | 0.051 |
| _smooth_estimates_cmip6 | 0.011 | 0.011 | 0.011 | 0.000 |
| _run_clustering_cmip6 | 0.098 | 0.093 | 0.100 | 0.002 |
| _incluster_reestimate_cmip6 | 4.775 | 4.766 | 4.793 | 0.009 |

### Previous historical checkpoints

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
## Decision Benchmark 2026-03-10T16:51:00.918262+00:00

- Git revision: `4cf18b2`
- Entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Method: `1 warmup + 5 measured runs` (warmups excluded from summary)
- Benchmark case: `medium`
- Config: `{'seed': 12345, 'n_years': 48, 'n_lat': 16, 'n_lon': 16, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 576, 'n_cells': 256, 'n_valid_cells': 256, 'n_years_complete': 48, 'k_lec': 33, 'k_edc': 32}`
- Total time summary: mean `11.051s`, min `10.902s`, max `11.221s`, std `0.132s`
- Peak memory summary: mean `0.613 GiB`, max `0.617 GiB` (`peak_process_tree_rss`, Δt=0.05s)
- Checks stable across measured runs: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _detrend_grid_fast | 5.225 | 5.160 | 5.369 | 0.076 |
| _monthly_annual_maxima | 0.001 | 0.000 | 0.001 | 0.000 |
| _compute_frechet_global | 1.013 | 0.976 | 1.101 | 0.049 |
| _run_local_estimation_cmip6 | 2.594 | 2.553 | 2.691 | 0.052 |
| _smooth_estimates_cmip6 | 0.003 | 0.003 | 0.003 | 0.000 |
| _run_clustering_cmip6 | 0.040 | 0.036 | 0.043 | 0.003 |
| _incluster_reestimate_cmip6 | 2.133 | 2.098 | 2.164 | 0.026 |

- Reference checks: `{'frechet_min': 0.15719810969238618, 'frechet_max': 517.3014620256827, 'labels_edc_sum': 4090, 'est_mean_a': 0.07405666957724144, 'est_mean_b': 6.947762730241841, 'est_mean_gamma': 0.09303462842023175}`
## Run 2026-03-11T20:50:52.086251+00:00

- Git revision: `225b588`
- Script entrypoint: `scripts.reproduce_fig9.main`
- Pipeline entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `4.587s`
- Config: `{'seed': 12345, 'n_years': 4, 'n_lat': 3, 'n_lon': 3, 'n_workers': 1, 'year_start': 1980, 'dpi': 300, 'generate_plots': False, 'backend': 'python', 'suppress_script_output': True}`
- Derived: `{'n_months': 48, 'n_cells': 9, 'n_valid_cells': 9, 'n_years_complete': 4, 'k_lec': 3, 'k_edc': 2, 'saved_figure_count': 0}`
- Max memory: `236666880` bytes (`231120.0 KiB`, `225.703 MiB`, `0.220 GiB`; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.131 |
| _monthly_annual_maxima | 0.155 |
| _compute_frechet_global | 1.081 |
| _run_local_estimation_cmip6 | 2.253 |
| _smooth_estimates_cmip6 | 0.000 |
| _run_clustering_cmip6 | 0.255 |
| _incluster_reestimate_cmip6 | 0.678 |

- Checks: `{'frechet_min': 0.4039738356016117, 'frechet_max': 16.0, 'labels_edc_sum': 15, 'est_mean_a': 0.11202639737889099, 'est_mean_b': 9.294280877432646, 'est_mean_gamma': 0.2644610756069325}`
## Run 2026-03-11T20:53:08.449560+00:00

- Git revision: `225b588`
- Script entrypoint: `scripts.reproduce_fig9.main`
- Pipeline entrypoint: `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
- Total time: `18.160s`
- Config: `{'seed': 12345, 'n_years': 12, 'n_lat': 6, 'n_lon': 6, 'n_workers': 4, 'year_start': 1980, 'dpi': 300, 'generate_plots': False, 'backend': 'python', 'suppress_script_output': True}`
- Derived: `{'n_months': 144, 'n_cells': 36, 'n_valid_cells': 36, 'n_years_complete': 12, 'k_lec': 7, 'k_edc': 9, 'saved_figure_count': 0}`
- Max memory: `927809536` bytes (`906064.0 KiB`, `884.828 MiB`, `0.864 GiB`; peak_process_tree_rss, Δt=0.05s)

| Step | Seconds |
| --- | ---: |
| _detrend_grid_fast | 0.798 |
| _monthly_annual_maxima | 0.144 |
| _compute_frechet_global | 3.990 |
| _run_local_estimation_cmip6 | 10.388 |
| _smooth_estimates_cmip6 | 0.000 |
| _run_clustering_cmip6 | 0.275 |
| _incluster_reestimate_cmip6 | 2.534 |

- Checks: `{'frechet_min': 0.2296057314788301, 'frechet_max': 144.0, 'labels_edc_sum': 162, 'est_mean_a': 0.11799152428727161, 'est_mean_b': 9.744587339675194, 'est_mean_gamma': -0.1480983950616838}`

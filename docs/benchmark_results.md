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

## Rust Backend A/B Comparison — 2026-03-09

Revision `9f24784`. First Rust integration: PyO3 bindings for `neg_log_likelihood_sum` and LEC matrix.
The Rust backend is called per-evaluation from Python's L-BFGS-B optimizer.

### Python backend (`WEATHERISK_BACKEND=python`)

- Total time: mean `9.169s`, min `8.834s`, max `9.317s`, std `0.181s`
- Peak RSS: mean `0.769 GiB`, max `0.783 GiB`
- Checks stable: `True`
- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.01500266795292907, 'est_mean_b': 0.8624371652491263, 'est_mean_gamma': 0.02906508271908983}`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _run_local_estimation_cmip6 | 2.288 | 2.225 | 2.322 | 0.035 |
| _incluster_reestimate_cmip6 | 5.736 | 5.525 | 5.852 | 0.130 |

### Rust backend (`WEATHERISK_BACKEND=auto`, Rust detected)

- Total time: mean `12.320s`, min `12.106s`, max `12.516s`, std `0.148s`
- Peak RSS: mean `0.606 GiB`, max `0.618 GiB`
- Checks stable: `True`
- Reference checks: `{'frechet_min': 0.14875993551681177, 'frechet_max': 771.6627392983846, 'labels_edc_sum': 2485, 'est_mean_a': 0.015160195138307478, 'est_mean_b': 0.8625245785429787, 'est_mean_gamma': 0.0289673269060139}`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _run_local_estimation_cmip6 | 2.248 | 2.168 | 2.311 | 0.049 |
| _incluster_reestimate_cmip6 | 8.940 | 8.792 | 9.026 | 0.084 |

### Analysis

- **Memory**: Rust uses 21% less memory (0.606 vs 0.769 GiB) — the LEC matrix avoids materializing a 3D boolean tensor.
- **Speed regression**: Rust is 34% slower overall. The bottleneck is `_incluster_reestimate_cmip6` (+56%).
- **Root cause**: PyO3 boundary crossing overhead. The Python L-BFGS-B optimizer calls into Rust ~350 times per cell (7 finite-difference gradient evaluations × ~50 iterations). Each call converts NumPy arrays to Rust slices and back. The per-call FFI cost (~20µs) dominates compared to the ~5µs Python-native NumPy evaluation.
- **Next step**: Move the entire optimizer loop into Rust to eliminate FFI overhead. The L-BFGS-B solver, objective function, finite-difference gradient, and bounds checking should all run compiled, crossing PyO3 only once per cell.

## Rust nll_with_gradient A/B Comparison — 2026-03-09

Approach: Keep SciPy's Fortran L-BFGS-B but compute NLL + forward-difference
gradient in a single Rust FFI call (`nll_with_gradient`). This reduces FFI
crossings from ~4/iteration (1 f + 3 approx_fprime) to 1/iteration. The pure
Rust L-BFGS-B (`lbfgsb-rs-pure`) was 14× slower than SciPy's Fortran and was
abandoned.

### Python backend (`WEATHERISK_BACKEND=python`)

- Total time: mean `9.201s`, min `9.144s`, max `9.304s`, std `0.061s`
- Peak RSS: mean `0.755 GiB`, max `0.763 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _run_local_estimation_cmip6 | 2.311 | 2.291 | 2.380 | 0.035 |
| _incluster_reestimate_cmip6 | 5.731 | 5.701 | 5.769 | 0.027 |

### Rust backend (`WEATHERISK_BACKEND=rust`, nll_with_gradient)

- Total time: mean `8.320s`, min `8.274s`, max `8.417s`, std `0.050s`
- Peak RSS: mean `0.610 GiB`, max `0.616 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _run_local_estimation_cmip6 | 2.144 | 2.117 | 2.234 | 0.045 |
| _incluster_reestimate_cmip6 | 5.019 | 5.007 | 5.030 | 0.008 |

### Analysis

- **Total**: Rust is **9.6% faster** (8.320s vs 9.201s).
- **Incluster**: 12.4% faster (5.019s vs 5.731s) — the biggest bottleneck.
- **Local MLE**: 7.2% faster (2.144s vs 2.311s).
- **Memory**: 19.2% less (0.610 vs 0.755 GiB).
- **vs previous Rust attempt**: The per-call Rust NLL was 34% *slower*; `nll_with_gradient` is 9.6% *faster*. The `lbfgsb-rs-pure` crate was 14× slower than SciPy's Fortran L-BFGS-B in micro-benchmarks.

## Rust Kernel Optimisation — 2026-03-10

Deep optimisation of the inner `pairwise_density_summand` kernel in Rust.
All changes preserve mathematical equivalence (parity tests pass to 1e-12
element-wise).

**Optimisations applied:**

1. **CSE on `powf`**: 14 `powf()` calls per element → 2 base powers
   `z^(1/df)`, all other fractional exponents derived via mul/div.
2. **Closed-form `t_cdf`**: For integer `df/2` (covers default df=5),
   replaced iterative `beta_reg` with exact polynomial
   `I_z(a, 1/2) = 1 − √(1−z) Σ c_j z^j`. 5–8× faster than statrs.
3. **Fast `t_pdf`**: Precomputed log-coefficient (avoids 2× `ln_gamma`
   per call). Integer-power decomposition for half-integer exponents
   (e.g., `u^{-3.5} = 1/(u³√u)`) replaces one `powf` per call.
4. **Inlined `dtdiff`**: Reuses `t_pdf` values already computed for
   `dt_m1`/`dt_m2`, eliminating 2 redundant `t_pdf` evaluations.
5. **Covariance caching**: In `neg_log_likelihood_sum`, detects contiguous
   (x,y) blocks (from `np.repeat` in pair-building) and computes
   `cov_fkt_2d` once per spatial pair — 48× reduction for n_sim=48.
6. **`alpha=1.0` fast path**: Skips `powf(alpha)` in `cov_fkt_2d`
   when `alpha=1.0` (the pipeline default).

### Micro-benchmark: kernel ns/element

| Metric | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| `neg_log_likelihood_sum` µs/call | 3366 | 471 | **7.1×** |
| Per-element ns | 369 | 52 | **7.1×** |
| `nll_with_gradient` µs/call | 13700 | 1858 | **7.4×** |
| NLL/trivial ratio | 742× | 118× | 6.3× more efficient |

### Scaling: full optimization time by cluster size

| Cells | Before | After | Speedup |
| ---: | ---: | ---: | ---: |
| 10 | 0.184s | 0.023s | **8.0×** |
| 20 | 0.496s | 0.087s | **5.7×** |
| 30 | 1.243s | 0.275s | **4.5×** |
| 50 | 4.702s | 0.876s | **5.4×** |

### Python backend (`WEATHERISK_BACKEND=python`)

- Total time: mean `9.281s`, min `8.955s`, max `9.506s`, std `0.191s`
- Peak RSS: mean `0.769 GiB`, max `0.791 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _run_local_estimation_cmip6 | 2.341 | 2.224 | 2.376 | 0.059 |
| _incluster_reestimate_cmip6 | 5.740 | 5.571 | 5.858 | 0.113 |

### Rust backend (`WEATHERISK_BACKEND=rust`, optimised kernel)

- Total time: mean `3.270s`, min `3.214s`, max `3.400s`, std `0.070s`
- Peak RSS: mean `0.600 GiB`, max `0.610 GiB`
- Checks stable: `True`

| Step | Mean (s) | Min (s) | Max (s) | Std (s) |
| --- | ---: | ---: | ---: | ---: |
| _compute_frechet_global | 1.059 | 1.025 | 1.134 | 0.041 |
| _run_local_estimation_cmip6 | 0.702 | 0.690 | 0.728 | 0.014 |
| _incluster_reestimate_cmip6 | 1.377 | 1.365 | 1.403 | 0.014 |

### Analysis

- **Total**: Rust is **2.84× faster** than Python (3.270s vs 9.281s, 64.8% reduction).
- **Incluster**: **4.17× faster** (1.377s vs 5.740s, 76.0% reduction).
- **Local MLE**: **3.34× faster** (0.702s vs 2.341s, 70.0% reduction).
- **Memory**: 22% less (0.600 vs 0.769 GiB).
- **vs previous Rust**: **2.54× faster** than the pre-optimisation Rust (3.270s vs 8.320s).
- **Bottleneck shift**: `_compute_frechet_global` (GEV fitting, pure SciPy) is now 32% of total time. The NLL kernel is no longer the dominant cost.

## Historical Optimization Checkpoints

The runs below were useful during implementation, but they are one-off
measurements rather than repeated-run decision benchmarks.

### Memory Optimisation & Scaling Test — 2026-03-10

Applied three memory optimisations to bring peak RSS under 8 GiB for
the full 18,432-cell CMIP6 grid:

1. **Condensed LEC dissimilarity**: compute the upper-triangle vector
   directly (both Python and Rust backends), skipping the full n×n
   matrix allocation.
2. **Pre-filtered pair expansion**: filter pairs by `max_dist` BEFORE
   repeating them across `n_sim` years.
3. **Block-based NLL evaluation**: evaluate the negative log-likelihood
   in blocks of 5,000 pairs, capping per-call memory regardless of
   cluster size.

#### Scaling comparison: 16×16 vs 32×32

Both runs: 48 years, 4 workers, Rust backend, 1 warmup + 5 measured.

| Metric | 16×16 (256 cells) | 32×32 (1,024 cells) | Ratio | Expected (linear) |
| --- | ---: | ---: | ---: | ---: |
| n_cells | 256 | 1,024 | 4.0× | 4.0× |
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

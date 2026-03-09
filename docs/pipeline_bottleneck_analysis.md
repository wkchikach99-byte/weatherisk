# Pipeline Data Flow and Bottleneck Analysis

This note replaces the older pre-optimization bottleneck write-up. It is based on the current code in `weatherisk.cpc_pipeline`, `weatherisk.cmip6_pipeline`, `weatherisk.density`, `weatherisk.clustering`, and the benchmark evidence recorded in `docs/benchmark_results.md`.

Two distinctions matter:

1. The CMIP6 Figure 9 path has received several performance improvements and now behaves differently from the older notes.
2. The CPC real-data path still contains slower serial implementations for some of the same conceptual steps.

## Executive Summary

| Pipeline | Main bottleneck now | Secondary bottleneck | Important qualification |
| --- | --- | --- | --- |
| CMIP6 (`weatherisk cmip6`) | Step 3 local MLE | Step 6 in-cluster re-estimation | GEV fitting, EDC matrix, and cluster re-estimation are no longer purely serial in this path |
| CPC (`weatherisk maps`) | Step 3 local MLE | Step 6 in-cluster re-estimation and Step 5 EDC matrix | CPC has not yet inherited the same optimizations as CMIP6 |

In both pipelines, the dominant computational cost is still repeated optimization of the pairwise composite likelihood. The core reason is unchanged: `pairwise_density_summand()` is evaluated many times inside multi-start L-BFGS-B runs.

## What Changed Relative to the Older Note

The following statements are no longer true for the current CMIP6 implementation:

- GEV fitting is not always serial anymore. `cmip6_pipeline._compute_frechet_global()` uses a `multiprocessing.Pool` when `n_workers > 1`.
- The CMIP6 EDC matrix is no longer built with the older Python pair loop. `cmip6_pipeline._edc_matrix_flat()` now uses `scipy.spatial.distance.cdist`.
- CMIP6 in-cluster re-estimation is no longer serial-only. `_incluster_reestimate_cmip6()` parallelizes across clusters when `n_workers > 1`.
- The old STL bottleneck is gone. `_detrend_grid_fast()` is vectorized and negligible in benchmark runs.
- The CMIP6 path now supports checkpointing after Steps 2, 3, and 4, reducing rerun cost after failures.

Those improvements do not automatically apply to the CPC pipeline, which still uses separate implementations.

## Current CMIP6 Execution Profile

### Step 0: Data discovery and loading

- Entrypoint: `cmip6_data.load_monthly_precipitation`
- Operational role: find files locally, on HPC, or via ESGF if needed
- Compute role: low relative to later statistical steps
- Bottleneck status: not primary

### Step 1a: De-trending

- Function: `cmip6_pipeline._detrend_grid_fast`
- Method: monthly climatology plus centered running-mean trend removal
- Current status: no longer a bottleneck

Benchmark evidence in `docs/benchmark_results.md` shows the vectorized implementation reduced this step from the earlier multi-second baseline to effectively negligible time on the synthetic full-pipeline benchmark.

### Step 1b: Annual maxima

- Function: `cmip6_pipeline._monthly_annual_maxima`
- Method: keep complete years, then take the maximum over the 12 monthly values
- Bottleneck status: negligible

### Step 2: GEV fit and Fréchet transform

- Function: `cmip6_pipeline._compute_frechet_global`
- Method: fit a GEV at each valid cell, transform to unit Fréchet, clip extreme tails
- Parallelism: optional process pool over valid cells
- Bottleneck status: moderate, not dominant

This step still performs many independent `genextreme.fit` calls, so it remains a meaningful cost center at scale, but it is no longer correctly described as a serial bottleneck in the CMIP6 path.

### Step 3: Local MLE

- Functions: `cmip6_pipeline._run_local_estimation_cmip6`, `_local_mle_one_cmip6`
- Core kernel: `density.pairwise_density_summand`
- Parallelism: optional process pool over grid cells
- Bottleneck status: dominant

This remains the main runtime bottleneck because each valid cell triggers:

1. Neighborhood construction in grid-point-distance space.
2. Assembly of repeated pair arrays `zi`, `zj`, `xl`, `yl`.
3. Multiple L-BFGS-B restarts.
4. Many evaluations of `pairwise_density_summand` per optimizer iteration.

The cost is dominated less by orchestration and more by the repeated statistical objective evaluation. Even after replacing SciPy's higher-overhead Student-`t` wrappers inside `density.py`, this step still scales with the number of cells, neighbors, years, optimizer iterations, and starts.

### Step 4: Spatial smoothing

- Function: `cmip6_pipeline._smooth_estimates_cmip6`
- Method: neighborhood mean for `a` and `b`, wrapped-angle mean for `γ`
- Bottleneck status: small

This is a serial loop, but the work per cell is small relative to the optimization steps.

### Step 5: LEC and EDC clustering

#### Step 5a: LEC matrix

- Function: `clustering.calc_distance_ellipses`
- Method: rasterized ellipse-overlap dissimilarity with chunked vectorization
- Bottleneck status: moderate and memory-sensitive

The LEC matrix still has an `O(n^2)` footprint and uses large intermediate boolean arrays. On small synthetic benchmarks it is minor, but on the full T63 grid it remains the most memory-sensitive part of clustering.

#### Step 5b: EDC matrix

- Function: `cmip6_pipeline._edc_matrix_flat`
- Method: rank each cell, then use `cdist(..., metric="cityblock")`
- Bottleneck status: improved; no longer a primary concern relative to Steps 3 and 6

The old note correctly identified the Python-loop version as wasteful, but that statement is now historical for CMIP6.

#### Step 5c: Hierarchical clustering

- Function: `clustering.clustering`
- Method: SciPy average-linkage clustering on the condensed distance vector
- Bottleneck status: secondary

### Step 6: In-cluster re-estimation

- Function: `cmip6_pipeline._incluster_reestimate_cmip6`
- Core kernel: `density.pairwise_density_optim` → `pairwise_density_summand`
- Parallelism: optional process pool over clusters
- Bottleneck status: second-largest hotspot

This step still repeats expensive composite-likelihood optimization, but it now parallelizes across clusters. That makes the old “serial over clusters” diagnosis obsolete for the CMIP6 path. The step remains expensive because cluster-level optimizations can still involve many cell pairs and repeated optimizer calls.

### Step 7: Plotting and writes

- Functions: `plot_figure9`, `np.savez_compressed`
- Bottleneck status: minor

### Operational bottleneck: restart cost after failures

Historically, a failed CMIP6 run could lose all progress before `pipeline_results.npz` was written. The current checkpoint system mitigates this by saving:

- `step2.npz` after Fréchet transformation
- `step3.npz` after local MLE
- `step4.npz` after smoothing

This is not a speed optimization in the strict sense, but it materially reduces wasted cluster time on long HPC runs.

## Current CPC Execution Profile

The CPC path implements the same scientific workflow but with different performance characteristics.

### Step 1: Load and subset CPC data

- Function: `cpc_pipeline._load_subregion`
- Bottleneck status: modest I/O, usually not dominant

### Step 2: GEV and Fréchet transform

- Function: `cpc_pipeline._compute_frechet`
- Parallelism: serial loop over retained land cells
- Bottleneck status: meaningful, but usually still below local MLE on realistic runs

Unlike CMIP6, the CPC path does not currently fan this step out across worker processes.

### Step 3: Local MLE

- Function: `cpc_pipeline._run_local_estimation`
- Parallelism: serial loop over retained cells
- Bottleneck status: dominant

This is the single clearest remaining performance gap between the CPC and CMIP6 implementations. The conceptual workload is the same as in CMIP6 local MLE, but the CPC path still evaluates it in a single Python loop over cells.

### Step 4: Smoothing

- Function: `cpc_pipeline._smooth_estimates`
- Bottleneck status: small

### Step 5: Clustering

#### LEC

- Function: `clustering.calc_distance_ellipses`
- Bottleneck status: moderate for larger retained grids

#### EDC

- Function: `cpc_pipeline._edc_matrix`
- Parallelism: serial Python loop over rows
- Bottleneck status: avoidable secondary hotspot

This is another place where the CPC path is behind the CMIP6 path: it still uses a row loop rather than the `cdist`-based implementation.

### Step 6: In-cluster re-estimation

- Function: `cpc_pipeline._incluster_reestimate`
- Parallelism: serial loop over clusters
- Bottleneck status: major secondary hotspot

### Step 7 and 8: Risk metrics, GDP weighting, map generation

- Functions: `_cluster_risk`, GDP logic inside `run_cpc_pipeline`, `generate_maps`
- Bottleneck status: not primary

## Root Cause of the Remaining Hotspots

Across both pipelines, the deepest cost is still the same statistical inner loop.

### `pairwise_density_summand()` remains the core hot path

Each optimizer call repeatedly evaluates:

- `cov_fkt_2d`
- Student-`t` PDF and CDF terms
- derivative terms via `_dtdiff`
- many power and ratio transforms on large pair arrays

The current implementation has already reduced wrapper overhead in `density.py`, but the function is still called many times, and it is called inside optimization, which multiplies the cost.

### Why local and in-cluster MLE dominate

- They scale with the number of neighbors or cluster pairs.
- They rebuild pair arrays for each optimization target.
- They use multiple optimizer restarts.
- They allow up to `maxiter=10000`, which is conservative.

That combination dominates both CPU time and, for large problems, wall-clock time even after the lower-level Student-`t` optimization work.

## Benchmark-Backed Interpretation

The synthetic full-pipeline benchmarks in `docs/benchmark_results.md` support the current qualitative ranking for CMIP6.

- After the low-hanging-fruit optimization pass, the benchmark dropped from `5.755s` to `2.543s` on the comparable serial setup.
- On the optimized serial benchmark, `_run_local_estimation_cmip6` and `_incluster_reestimate_cmip6` are still among the largest steps.
- After Step 6 parallelization, `_incluster_reestimate_cmip6` dropped from `3.761s` to `2.680s` in the representative multi-worker benchmark, confirming that the step remains important even after improvement.

These are synthetic-scale comparisons rather than full-HPC production timings, so they should be used to compare code paths, not to infer exact wall-clock runtime on the full 96×192 grid.

## Current Documentation-Safe Conclusions

The following statements are precise and up to date.

- For CMIP6, the dominant computational cost is still local MLE, with in-cluster re-estimation second.
- For CMIP6, GEV fitting and in-cluster re-estimation are parallel-capable, so they should not be described as purely serial bottlenecks anymore.
- For CMIP6, EDC matrix construction has already been moved to a `cdist`-based implementation.
- For CPC, the main remaining bottlenecks are still the serial local-MLE loop, the serial EDC-matrix loop, and serial in-cluster re-estimation.
- The CPC and CMIP6 code paths should be documented separately when discussing performance.

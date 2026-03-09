# Code Execution Flow

This document maps the public entrypoints in `weatherisk` to the concrete
functions that run, the intermediate objects they produce, and the files
they write. It is meant to answer three practical questions:

1. Which command or script should be used for which workflow?
2. Which functions actually perform the computation?
3. Where do the CPC and CMIP6 implementations currently differ?

## Public Entry Points

| Entry point | Function path | Purpose |
| --- | --- | --- |
| `weatherisk validate` | `weatherisk.cli.validate` → `weatherisk.pipeline.run_pipeline` | Lightweight synthetic sanity check |
| Python API | `weatherisk.pipeline.run_nonstationary_pipeline` | Full non-stationary synthetic study |
| `weatherisk maps` | `weatherisk.cli.maps` → `weatherisk.cpc_pipeline.run_cpc_pipeline` | Main CPC real-data analysis |
| `weatherisk cmip6` | `weatherisk.cli.cmip6` → `weatherisk.cmip6_pipeline.run_cmip6_pipeline` | Figure 9 CMIP6 reproduction |
| `weatherisk risk-pipeline` | `weatherisk.cli.risk_pipeline` | Post-process a precomputed risk CSV |
| `scripts/reproduce_fig9.py` | Thin wrapper around `run_cmip6_pipeline` and `plot_figure9` | Scriptable Figure 9 run outside the package entrypoint |
| `weatherisk.benchmarks.run_hotpath_benchmark` | Wraps `run_cmip6_pipeline` with synthetic input and timing hooks | Benchmark the actual CMIP6 orchestration path |

## Shared Computational Kernels

Several pipelines call the same lower-level functions.

- `weatherisk.extremes.block_maxima`: block aggregation for daily data
- `weatherisk.extremes.fit_gev`: per-cell GEV fit using `scipy.stats.genextreme.fit`
- `weatherisk.extremes.to_frechet`: GEV-to-unit-Fréchet marginal transform
- `weatherisk.density.pairwise_density_summand`: pairwise composite log-density kernel
- `weatherisk.density.pairwise_density_optim`: global multi-start L-BFGS-B for one cluster
- `weatherisk.density.pairwise_density_optim_local`: local version used in the synthetic grid API
- `weatherisk.clustering.calc_distance_ellipses`: LEC dissimilarity matrix
- `weatherisk.clustering.clustering`: SciPy hierarchical clustering using average linkage
- `weatherisk.risk.compute_var` and `weatherisk.risk.compute_es`: cluster-level tail-risk summaries on the Fréchet scale

## Flow 1: `weatherisk validate`

The validation command is intentionally lightweight.

### Call chain

`weatherisk.cli.validate`
→ `weatherisk.parameters.get_preset`
→ `weatherisk.pipeline.run_pipeline`

### What `run_pipeline` actually does

1. Build a regular `Grid` on `[-5, 5] × [-5, 5]`.
2. Simulate a stationary max-stable field with `simulation.sim_expt_2d`.
3. Create synthetic local estimates by adding small noise around the known input `(a, b, γ)`.
4. Compute the ellipse-overlap dissimilarity matrix with `calc_distance_ellipses`.
5. Run average-linkage clustering and cut the dendrogram at a simple size-based `k`.
6. Optionally save `clusters.npy` and `estimates.npy`.

### Important implication

This command does not run local MLE, smoothing, EDC clustering, or in-cluster re-estimation. It is best viewed as a fast structural check of the clustering machinery.

## Flow 2: `run_nonstationary_pipeline`

This is the full synthetic methodology path and is closer to the Extremes-paper workflow.

### Call chain

`weatherisk.pipeline.run_nonstationary_pipeline`
→ `simulation.sim_expt_2d_nonstat`
→ `density.run_local_mle_parallel`
→ `estimation.smooth_local_estimates`
→ `clustering.calc_distance_ellipses`
→ `clustering.c_extrcoeff_matrix`
→ `estimation.calc_estimates_in_clusters`

### Step-by-step data flow

1. Evaluate preset functions `a_func`, `b_func`, and `g_func` on the grid.
2. Simulate a non-stationary max-stable process.
3. Run local MLE at every grid cell.
4. Smooth the local estimates.
5. Build the LEC dissimilarity matrix and cluster it.
6. Build the EDC dissimilarity matrix and cluster it.
7. Re-estimate `(a, b, γ)` inside each LEC and EDC cluster.
8. Return all intermediate arrays in a single result dictionary.

## Flow 3: `weatherisk maps`

This is the current main real-data climate-risk pipeline.

### Call chain

`weatherisk.cli.maps`
→ `weatherisk.cpc_pipeline.PipelineConfig`
→ `weatherisk.cpc_pipeline.run_cpc_pipeline`
→ optional `weatherisk.cpc_pipeline.generate_maps`

### Internal execution order

1. `_load_subregion`
   Loads yearly CPC NetCDF files, subsets the requested lat-lon box, and coarsens the grid.
2. `_compute_frechet`
   Flattens the grid, applies a land mask, fills remaining NaNs columnwise, computes annual block maxima, fits a GEV at each retained cell, transforms to unit Fréchet, and normalizes coordinates to `[-5, 5] × [-5, 5]`.
3. `_run_local_estimation`
   Runs `_local_mle_one` for every retained cell.
4. `_smooth_estimates`
   Moving-average smoothing for `a` and `b`, wrapped-angle smoothing for `γ`.
5. `_run_clustering`
   Builds the LEC matrix with `calc_distance_ellipses` and the EDC matrix with the CPC-specific `_edc_matrix`, then applies average-linkage clustering and the 30% quantile threshold rule.
6. `_incluster_reestimate`
   Re-estimates one `(a, b, γ)` triple per cluster using `pairwise_density_optim`.
7. `_cluster_risk`
   For each cluster, forms the spatial maximum over cells and computes VaR and ES on the Fréchet scale.
8. Optional GDP weighting
   If `gdp_path` is provided, each cell inherits its cluster ES and multiplies it by regridded GDP.
9. `generate_maps`
   Writes the standard PNG outputs to the configured output directory.

### Core returned objects

- `frechet`: shape `(n_blocks, n_land)`
- `estimates`, `smoothed`: shape `(n_land, 3)`
- `labels_lec`, `labels_edc`: one label per retained land cell
- `risk_lec`, `risk_edc`: list-of-dict cluster summaries
- `gdp_per_cell`, `risk_gdp_lec`, `risk_gdp_edc`: optional exposure outputs

## Flow 4: `weatherisk cmip6`

This path is specialized to the AWI-ESM-1-1-LR historical precipitation setup used for Figure 9.

### Call chain

`weatherisk.cli.cmip6`
→ `weatherisk.cmip6_pipeline.CMIP6Config`
→ `weatherisk.cmip6_pipeline.run_cmip6_pipeline`
→ optional `weatherisk.cmip6_pipeline.plot_figure9`

### Internal execution order

1. `load_monthly_precipitation`
   Loads or discovers the monthly CMIP6 NetCDF files.
2. `_detrend_grid_fast`
   Removes the long-term component with a vectorized monthly-climatology plus running-mean approach.
3. `_monthly_annual_maxima`
   Keeps only complete years and takes the maximum over the 12 monthly values in each year.
4. `_compute_frechet_global`
   Fits a GEV at each valid cell, optionally in parallel, and transforms to unit Fréchet.
5. `_grid_coords`
   Converts valid flat indices into row-column grid coordinates for grid-point-distance neighborhoods.
6. `_run_local_estimation_cmip6`
   Runs `_local_mle_one_cmip6` at every valid cell, optionally in parallel.
7. `_smooth_estimates_cmip6`
   Smooths local estimates across neighboring grid points.
8. `_run_clustering_cmip6`
   Computes LEC with `calc_distance_ellipses`, computes EDC with `_edc_matrix_flat`, clusters both, and selects `k` from the 30% quantile threshold.
9. `_incluster_reestimate_cmip6`
   Re-estimates cluster parameters for LEC and EDC clusters, optionally in parallel across clusters.
10. Save and plot
   Writes `pipeline_results.npz` and optionally generates the Figure 9 maps.

### Checkpoint and resume behavior

The current CMIP6 implementation can checkpoint intermediate stages.

- `step2.npz`: after GEV and Fréchet transformation
- `step3.npz`: after local MLE
- `step4.npz`: after smoothing

These checkpoints live under `output/cmip6_fig9/checkpoints/` unless `checkpoint_dir` is overridden in `CMIP6Config`. This does not accelerate successful runs directly, but it prevents long reruns after failures or OOM events before the final `pipeline_results.npz` write.

## Flow 5: `weatherisk risk-pipeline`

This path assumes the risk field already exists and only performs spatial post-processing.

### Call chain

`weatherisk.cli.risk_pipeline`
→ `weatherisk.risk_pipeline.load_and_grid`
→ `smooth_field`
→ `quantile_bands`
→ `connected_patches`
→ `merge_tiny_regions`
→ `remap_ids_to_sequential`
→ `compute_cluster_stats`

### Data flow

1. Load a CSV with at least `lat`, `lon`, `VaR_95`, and `ES_95`.
2. Reshape the table to 2-D ES and VaR grids.
3. Build a land mask from exposure.
4. Gaussian-smooth the risk field without leaking zeros into masked cells.
5. Discretize the field into quantile bands.
6. Find connected regions inside each band.
7. Merge tiny regions into the nearest larger region.
8. Relabel the final regions to `0..K-1` and compute regional summary statistics.

## CPC vs CMIP6: Current Implementation Split

The CPC and CMIP6 pipelines implement the same broad methodology, but they are not performance-equivalent.

- The CMIP6 path has explicit multiprocessing for GEV fitting, local MLE, and in-cluster re-estimation.
- The CMIP6 path uses the faster `_edc_matrix_flat` implementation backed by `scipy.spatial.distance.cdist`.
- The CMIP6 path now has checkpointing for Steps 2–4.
- The CPC path still uses a serial local-MLE loop, a CPC-specific serial `_edc_matrix`, and serial in-cluster re-estimation.

That split matters for both performance and documentation: optimization claims that are true for the CMIP6 path are not automatically true for the CPC path.
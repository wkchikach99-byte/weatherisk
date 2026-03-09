# weatherisk

## Dir structure

```
docs/ - the papers we are writing for a PhD, all in latex format
  phd_thesis/ - the main PhD thesis, in progress
  softwarex/ - a paper for the SoftwareX publication about the code we're writing for the PhD thesis, in progress
refernce_papers/ - papers we use as reference in our PhD work, with important methodology that we reuse and build on
weatherisk/ - the main python code base, which is the target code
r_code/ - the historical code base used in the reference papers by Justus, that we have refactored to python and built upon in python
scripts/ - scripts to produce some figures for our papers
data/ - data we've downloaded and using in our modelling. precipitation data, heat data, and more.
meeting_notes/ - meeting notes with the supervisors with feedback on the papers and the code

```

**Climate risk analysis via max-stable process clustering.**

`weatherisk` is a Python package that models the spatial dependence structure of climate extremes (heatwaves, heavy precipitation) using max-stable processes. It estimates local anisotropy parameters, clusters grid cells by dependence similarity, and computes tail-risk metrics (VaR, ES) per cluster.

The package consolidates an R + Shell/SLURM + Python workflow into a single, tested, pip-installable library with a CLI.

## Installation

```bash
pip install .
```

For development (includes pytest, ruff, mypy):

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.10.

## Quick start

```python
from weatherisk.pipeline import run_pipeline

result = run_pipeline(resolution=10, n_sim=20, seed=42)
print(f"Found {len(set(result['clusters']))} clusters")
```

Or via the CLI:

```bash
weatherisk validate --resolution 10 --n-sim 20 --seed 42
```

---

## Command Status

The current CLI has five subcommands, but they are not all at the same maturity level.

| Command | Status | Actual entrypoint | Notes |
|--------|--------|-------------------|-------|
| `weatherisk validate` | Implemented | `weatherisk.pipeline.run_pipeline` | Lightweight synthetic sanity check |
| `weatherisk maps` | Implemented | `weatherisk.cpc_pipeline.run_cpc_pipeline` | Main CPC real-data pipeline with risk maps |
| `weatherisk cmip6` | Implemented | `weatherisk.cmip6_pipeline.run_cmip6_pipeline` | Figure 9 reproduction path |
| `weatherisk risk-pipeline` | Implemented | `weatherisk.risk_pipeline` helpers | Post-process precomputed risk CSVs |
| `weatherisk risk` | Placeholder | `weatherisk.cli.risk` | Accepts arguments but exits with “Not yet implemented for real data” |

## Documentation Guide

- `docs/code_execution_flow.md`: function-level execution flow for each implemented pipeline
- `docs/pipeline_bottleneck_analysis.md`: current bottlenecks, split by CMIP6 and CPC paths
- `docs/methodology_notes.md`: methodological interpretation notes for the CPC maps pipeline
- `docs/benchmark_results.md`: benchmark evidence for the optimized CMIP6 Figure 9 path

## Pipeline Flows

### Flow 1: Lightweight Validation (`weatherisk validate`)

`weatherisk validate` is a lightweight synthetic sanity check. It does not run the full paper-style local-MLE pipeline.

```
Parameter preset
  →  Simulate stationary max-stable field
  →  Inject small noise around known (a, b, γ)
  →  Ellipse-overlap dissimilarity matrix
  →  Average-linkage clustering
  →  Save clusters / estimates
```

**Usage:**

```bash
weatherisk validate --params stripes --resolution 10 --n-sim 20 --seed 42
weatherisk validate --output-dir output/validate
```

**Python API:**

```python
from weatherisk.parameters import get_preset
from weatherisk.pipeline import run_pipeline

p = get_preset("stripes")
result = run_pipeline(
    resolution=p.resolution,
    n_sim=p.n_sim,
    df=p.df,
    alpha=p.alpha,
    seed=42,
)
clusters = result["clusters"]
linkage = result["linkage"]
```

For the full synthetic methodology from the paper, use `run_nonstationary_pipeline` from Python instead of `weatherisk validate`.

### Flow 2: Full Synthetic Study (`run_nonstationary_pipeline`)

This is the paper-style synthetic workflow used to reproduce the non-stationary simulation study.

```
Preset parameter fields a(x, y), b(x, y), γ(x, y)
  →  Simulate non-stationary max-stable process
  →  Local MLE at each grid cell
  →  Spatial smoothing
  →  LEC clustering (D₂)
  →  EDC clustering (D₁)
  →  In-cluster re-estimation
```

**Python API:**

```python
from weatherisk.pipeline import run_nonstationary_pipeline

result = run_nonstationary_pipeline(
    preset="paper_stripes",
    resolution=51,
    n_sim=100,
    n_workers=4,
)
```

### Flow 3: CPC Real-Data Maps (`weatherisk maps`)

This is the main real-data risk-analysis pipeline currently implemented in the package.

```
CPC daily NetCDF files
  →  Sub-region selection and coarsening
  →  Annual block maxima
  →  GEV fit per retained cell
  →  Unit-Fréchet transform
  →  Local pairwise CL MLE of (a, b, γ)
  →  Spatial smoothing
  →  LEC clustering (ellipse overlap D₂)
  →  EDC clustering (madogram D₁)
  →  In-cluster re-estimation
  →  Cluster-level VaR / ES on the Fréchet scale
  →  Optional GDP-weighted cell-level risk
  →  Cartopy map generation
```

**Usage:**

```bash
weatherisk maps --variable precip --gdp-path data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc
weatherisk maps --variable tmax --file-prefix tmax
weatherisk maps --variable precip --lat-range 30 65 --lon-range 5 55 --coarsen 4
```

**Output:** 11 PNG maps in `docs/figures/`, including cluster maps, parameter fields, ES maps, GDP exposure, and GDP-weighted risk maps.

### Flow 4: CMIP6 Figure 9 Reproduction (`weatherisk cmip6`)

This is the dedicated AWI-ESM-1-1-LR Figure 9 pipeline.

```
Monthly CMIP6 precipitation
  →  Vectorized de-trending
  →  Annual maxima of monthly data
  →  GEV fit per valid grid cell
  →  Unit-Fréchet transform
  →  Local MLE in grid-point distance units
  →  Spatial smoothing
  →  LEC / EDC clustering
  →  In-cluster re-estimation
  →  Save pipeline_results.npz
  →  Plot Figure 9 maps
```

The CMIP6 path also supports checkpoint-and-resume for Steps 2–4 through `output/cmip6_fig9/checkpoints/`.

**Usage:**

```bash
weatherisk cmip6
weatherisk cmip6 --workers 16
weatherisk cmip6 --data-dir /pool/data/CMIP6/.../pr/gn/
```

### Flow 5: Risk-Map Post-Processing (`weatherisk risk-pipeline`)

This path does not fit max-stable models. It post-processes an already computed risk grid.

```
CSV with lat / lon / VaR_95 / ES_95
  →  Reshape to 2-D grid
  →  Gaussian smoothing
  →  Quantile banding
  →  Connected-component labeling
  →  Merge tiny regions
  →  Sequential relabeling
  →  Per-region summary statistics
```

**Usage:**

```bash
weatherisk risk-pipeline --csv data/risk_map_grid.csv --bands 6 --sigma 0.8
```

## Library-Only Scalable Helpers

`weatherisk.scalable` contains coarse-grid proxy helpers such as `downsample_estimates()` and `propagate_cluster_labels()`. These are available for experimentation, but the current CLI pipelines do not call them by default.

---

## Module Reference

| Module | Purpose |
|--------|---------|
| `grid` | 2-D regular grid, index conversions (`grid_number`, `number_grid`, `koord_num`, `number_koord`), column-major (Fortran-order) flattened coordinates (`X_flat`, `Y_flat`), `rad()`/`deg()` |
| `covariance` | Stationary (`cov_fkt_2d`) and non-stationary (`cov_fkt_2d_nonstat2`) anisotropic covariance; extremal-coefficient conversions (`cov_to_ec`, `ec_to_cov`) |
| `simulation` | Max-stable process simulation via spectral representation (`sim_expt_2d`, `sim_expt_2d_nonstat`) |
| `density` | Pairwise composite likelihood (`pairwise_density_summand`); global and local MLE optimisation with R-matching heuristics (parscale, gamma-wrapping retry, boundary-proximity retry) |
| `estimation` | Spatial smoothing of local estimates (angular wrapping for γ); in-cluster re-estimation with average log-likelihood computation |
| `clustering` | Ellipse-overlap dissimilarity (`calc_distance_ellipses`); madogram-based Saunders method (`c_extrcoeff_matrix`, column-major flatten); hierarchical clustering; threshold-based *k* selection |
| `extremes` | Block-maxima extraction, GEV fitting, Fréchet transform |
| `risk` | VaR, ES, per-cluster risk aggregation |
| `risk_pipeline` | Risk-map loading, quantile banding, connected-component clustering, region statistics |
| `cpc_pipeline` | End-to-end CPC real-data pipeline: load CPC files → GEV → Fréchet → LEC/EDC clustering → risk → Cartopy maps |
| `map_plotting` | Filled-region Cartopy maps (pcolormesh): cluster maps, parameter fields, risk choropleths, summary panels |
| `gdp` | Gridded GDP exposure loading (Kummu et al. 2018), regridding to pipeline grid, cell-level extraction |
| `netcdf` | NetCDF climate data ingestion (CPC, ERA5), longitude wrapping |
| `scalable` | Library-only helpers for coarse-grid proxy clustering and label propagation |
| `parameters` | Named parameter presets (`stripes`, `bigsmall`, `rotate`) as dataclasses |
| `plotting` | Heatmaps, cluster maps, dendrograms, choropleth, bar charts (synthetic data) |
| `pipeline` | Synthetic orchestration: lightweight validation (`run_pipeline`) and full non-stationary study (`run_nonstationary_pipeline`) |
| `cli` | Click-based CLI entry point with `validate`, `maps`, `cmip6`, `risk-pipeline`, and placeholder `risk` |
| `io` | CSV, NumPy, and RDS I/O helpers |

## Parameter Presets

Three presets reproduce scenarios from the original paper:

| Preset | Resolution | Description |
|--------|-----------|-------------|
| `stripes` | 10 | Semi-major axis `b` varies linearly with the x-coordinate, creating stripe-like anisotropy regions |
| `bigsmall` | 51 | Semi-minor axis `a` varies radially from the grid centre (large centre, small edges) |
| `rotate` | 51 | Rotation angle γ varies linearly with x, creating a smooth rotation field |

```python
from weatherisk.parameters import get_preset

p = get_preset("stripes")
print(p.resolution, p.df, p.alpha)
```

## Mathematical Background

### Covariance Model

The stationary anisotropic covariance function is:

$$C(\mathbf{h}) = \exp\!\bigl(-\|\Sigma^{-1}\mathbf{h}\|^{\alpha}\bigr)$$

where $\mathbf{h} = (x, y)^T$ is the spatial lag, $\alpha \in (0, 2]$ is the smoothness exponent, and $\Sigma^{-1}$ encodes ellipse axes $a$, $a + b$ and rotation $\gamma$:

$$\Sigma^{-1} = R(\gamma) \begin{pmatrix} 1/a & 0 \\ 0 & 1/(a+b) \end{pmatrix} R(\gamma)^T$$

The non-stationary variant blends two anisotropy matrices $M_1, M_2$ (the inverse covariance matrices at two locations) via the harmonic mean:

$$\Omega = \bigl(\tfrac{1}{2}(M_1^{-1} + M_2^{-1})\bigr)^{-1}$$

$$C_{\text{ns}}(\mathbf{h}) = \min\!\Bigl(1,\; \sqrt{\det(\Omega)\,a_1(a_1+b_1)\,a_2(a_2+b_2)}\; \exp\!\bigl(-\sqrt{\mathbf{h}^T\Omega\,\mathbf{h}}^{\,\alpha}\bigr)\Bigr)$$

### Extremal Coefficient

The extremal coefficient $\theta \in [1, 2]$ links pairwise dependence to the covariance $\rho$ via the Student-$t$ CDF:

$$\theta(\rho) = 2\, T_{df+1}\!\left(\sqrt{\frac{(df+1)(1-\rho)}{1+\rho}}\right)$$

At $\rho = 0$ with finite $df$, $\theta < 2$ (residual tail dependence). Complete independence ($\theta = 2$) is reached only as $df \to \infty$.

### Clustering

Two methods are implemented:

1. **Ellipse dissimilarity** — Jaccard-like overlap between estimated anisotropy ellipses, followed by hierarchical agglomerative clustering (average linkage).
2. **Saunders method** — rank-based madogram estimation of pairwise extremal coefficients, converted to a dissimilarity matrix.

## Testing

```bash
# Run tests
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=weatherisk --cov-report=term-missing
```

105 tests covering all modules:

- **55 unit tests** for all modules using small synthetic grids (5×5 and 10×10) with fixed seeds.
- **28 cross-validation tests** (`test_r_cross_validation.py`) that load reference CSV data generated by the original R code and assert numerical equivalence at every pipeline stage (grid, covariance, density, madogram, clustering, log-likelihood).
- **22 R-parity fix tests** (`test_r_parity_fixes.py`) verifying the five algorithmic fixes that align Python with R's optimisation heuristics:
  - Column-major (Fortran-order) grid indexing matching R's `grid_number`/`number_grid`
  - `X_flat`/`Y_flat` properties using Fortran-order flattening
  - Parameter rescaling (`parscale = (upper - lower) / 100`) in both local and global optimisers
  - Gamma-wrapping retry when γ hits the ±π/2 boundary
  - Boundary-proximity retry (up to 5 extra LHS starts when parameters are near bounds)
  - Average log-likelihood computation in `calc_estimates_in_clusters` (column 4)
  - Column-major reshape in `c_extrcoeff_matrix` matching R's `as.vector()`

R reference data is generated by `tests/generate_r_reference.R` (requires R ≥ 4.0 with `lhs` package).

## License

MIT

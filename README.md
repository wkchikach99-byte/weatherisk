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

## Pipeline Flows

### Flow 1: Method Validation (`weatherisk validate`)

Reproduces the clustering method on **synthetic data** to verify correctness.

```
Parameter preset  →  Simulate max-stable process
  →  Local parameter estimation (pairwise composite likelihood)
  →  Spatial smoothing of estimates
  →  Ellipse-shape dissimilarity matrix
  →  Hierarchical agglomerative clustering
  →  In-cluster re-estimation
  →  Comparison plots
```

**Usage:**

```bash
# Using a named parameter preset (stripes, bigsmall, rotate)
weatherisk validate --params stripes --resolution 10 --n-sim 20 --seed 42

# Custom output directory
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
clusters = result["clusters"]   # 1-D array, one label per grid cell
linkage  = result["linkage"]    # scipy linkage matrix for dendrograms
```

**Grid sizes:** 10–51 points per axis (100–2 601 cells). Runs in minutes on a laptop.

---

### Flow 2: Real-Data Risk Analysis (`weatherisk risk`)

Applies the method to observed climate data (CPC, ERA5, CMIP6) and computes risk metrics.

```
NetCDF climate data
  →  Extract block maxima (annual/seasonal)
  →  Fit GEV marginals per grid cell
  →  Transform to unit Fréchet margins
  →  Local pairwise MLE of (a, b, γ)        [embarrassingly parallel]
  →  Spatial smoothing
  →  Scalable clustering (coarse-grid proxy)
  →  In-cluster re-estimation
  →  Compute VaR / ES per cluster
  →  Produce cluster map, risk choropleth, bar charts
```

**Usage:**

```bash
weatherisk risk --netcdf data/cpc_tmax.nc --hazard heat -k 25 --workers 8
```

**Python API (step by step):**

```python
import numpy as np
from weatherisk.netcdf import load_climate_data
from weatherisk.extremes import block_maxima, fit_gev, to_frechet
from weatherisk.grid import Grid
from weatherisk.density import pairwise_density_optim_local
from weatherisk.estimation import smooth_local_estimates
from weatherisk.clustering import calc_distance_ellipses, clustering
from weatherisk.risk import compute_var, compute_es

# 1. Load climate data
ds = load_climate_data("data/cpc_tmax.nc", variable="tmax")

# 2. Extreme value analysis per grid cell
bm = block_maxima(daily_series, block_size=365)
loc, scale, shape = fit_gev(bm)
frechet = to_frechet(bm, loc, scale, shape)

# 3. Local estimation (parallelisable)
estimates = pairwise_density_optim_local(
    sim_data=frechet, df=5.0, alpha=1.0,
    x=0, y=0, grid=grid, neighbourhood=3,
)

# 4. Smooth and cluster
smoothed = smooth_local_estimates(estimates, smoothing_dist=2, grid=grid)
D = calc_distance_ellipses(smoothed, res=grid.resolution)
hc = clustering(D, method="average")

# 5. Risk metrics
var_95 = compute_var(cluster_data, p=0.95)
es_95  = compute_es(cluster_data, p=0.95)
```

**Grid sizes:** up to 360 × 720 = 259 200 cells (0.5° global). Uses the coarse-grid proxy strategy for scalable clustering (see below).

---

### Flow 3: Risk Pipeline (`weatherisk risk-pipeline`)

Lightweight post-processing of **pre-computed risk maps**. Takes a CSV with lat/lon/VaR/ES columns, smooths the fields, creates contiguous risk regions via quantile banding + connected components, and computes per-region statistics.

**Usage:**

```bash
weatherisk risk-pipeline --csv data/risk_map_grid.csv --bands 6 --sigma 0.8
```

**Python API:**

```python
from weatherisk.risk_pipeline import (
    load_and_grid, smooth_field, quantile_bands,
    connected_patches, merge_tiny_regions,
    remap_ids_to_sequential, compute_cluster_stats,
)

data = load_and_grid("data/risk_map_grid.csv")
ES_s = smooth_field(data["ES"], sigma=0.8, land_mask=data["land_mask"])
bands, thresholds = quantile_bands(ES_s, n_bands=6)
clusters = connected_patches(bands, min_patch=30)
clusters = merge_tiny_regions(clusters, data["lon_grid"], data["lat_grid"])
clusters, K = remap_ids_to_sequential(clusters)
stats = compute_cluster_stats(clusters_df, data["df"])
```

---

### Flow 4: CPC Maps (`weatherisk maps`)

The primary pipeline for real-data analysis. Loads NOAA CPC daily NetCDF files (precipitation or temperature), computes annual block maxima, fits GEV marginals, transforms to unit Fréchet, estimates local anisotropy parameters, runs LEC and EDC clustering, and outputs filled-region Cartopy maps.

Optionally weights risk by GDP exposure using the Kummu et al. (2018) gridded GDP PPP dataset.

```
CPC daily NetCDF files (20 years)
  →  Sub-region selection & coarsening (~2° grid)
  →  Annual block maxima → GEV fit → Fréchet transform
  →  Local pairwise composite-likelihood MLE of (a, b, γ)
  →  Spatial smoothing
  →  LEC clustering (ellipse-overlap D₂) + EDC clustering (madogram D₁)
  →  30%-quantile threshold on dissimilarity matrix → k clusters
  →  In-cluster re-estimation
  →  Tail-risk: ES₉₅ of Fréchet spatial block max per cluster
  →  [Optional] GDP exposure × ES₉₅ per cell
  →  Filled-region Cartopy maps (PNG)
```

**Usage:**

```bash
# Precipitation with GDP exposure weighting
weatherisk maps --variable precip --gdp-path data/gdp/GDP_PPP_1990_2015_5arcmin_v2.nc

# Temperature without GDP
weatherisk maps --variable tmax --file-prefix tmax

# Custom region and resolution
weatherisk maps --variable precip --lat-range 30 65 --lon-range 5 55 --coarsen 4
```

**Output:** 11 PNG maps in `docs/figures/` including cluster maps, parameter fields, tail-risk intensity (ES₉₅), GDP exposure, and exposure-weighted risk.

**Performance:** ~100 seconds on Apple M2 (16 GB RAM) for 384 land cells × 20 years.

---

## Scalable Clustering Strategy

The full dissimilarity matrix for a 0.5° global grid (259K cells) requires ~500 GB.
`weatherisk` provides a **coarse-grid proxy** approach:

1. Estimate (a, b, γ) at full resolution (embarrassingly parallel).
2. Downsample estimates to ~2° (16 200 cells) via `downsample_estimates()`.
3. Compute the full dissimilarity matrix at coarse resolution (fits in memory).
4. Hierarchical clustering on the coarse grid → *k* clusters.
5. Propagate labels to fine grid via `propagate_cluster_labels()`.
6. Re-estimate parameters within each cluster.

```python
from weatherisk.scalable import downsample_estimates, propagate_cluster_labels

coarse_est = downsample_estimates(fine_estimates, fine_shape=(360, 720), coarse_shape=(90, 180))
# ... cluster on coarse_est ...
fine_labels = propagate_cluster_labels(coarse_labels, coarse_shape=(90, 180), fine_shape=(360, 720))
```

For HPC, local estimation can be distributed across SLURM job arrays using `chunk_indices()`, `save_chunk()`, and `load_chunk()`.

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
| `cpc_pipeline` | End-to-end CPC real-data pipeline: load NetCDF → GEV → Fréchet → LEC/EDC clustering → risk → Cartopy maps |
| `map_plotting` | Filled-region Cartopy maps (pcolormesh): cluster maps, parameter fields, risk choropleths, summary panels |
| `gdp` | Gridded GDP exposure loading (Kummu et al. 2018), regridding to pipeline grid, cell-level extraction |
| `netcdf` | NetCDF climate data ingestion (CPC, ERA5), longitude wrapping |
| `scalable` | Coarse-grid proxy clustering, chunked parallel estimation, checkpoint/resume |
| `parameters` | Named parameter presets (`stripes`, `bigsmall`, `rotate`) as dataclasses |
| `plotting` | Heatmaps, cluster maps, dendrograms, choropleth, bar charts (synthetic data) |
| `pipeline` | End-to-end orchestration for synthetic validation (`run_pipeline`) |
| `cli` | Click-based CLI entry point with four subcommands |
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

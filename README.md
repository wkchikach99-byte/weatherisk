# weatherisk

[![Tests](https://github.com/wkchikach99-byte/weatherisk/actions/workflows/tests.yml/badge.svg)](https://github.com/wkchikach99-byte/weatherisk/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20108563.svg)](https://doi.org/10.5281/zenodo.20108563)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

## Dir structure

```
docs/ - papers and thesis chapters (LaTeX)
  phd_thesis/ - the main PhD thesis, in progress
  software_x/ - the SoftwareX paper about this package
reference_papers/ - reference papers used in the PhD work
weatherisk/ - the main Python package source code
r_code/ - the original R reference implementation (Contzen et al.)
scripts/ - figure-reproduction and benchmarking scripts
data/ - downloaded climate data (precipitation, heat)
tests/ - pytest test suite
notebooks/ - quickstart Jupyter notebook
meeting_notes/ - supervisor meeting notes (not in public repo)
```

**Climate risk analysis via max-stable process clustering.**

`weatherisk` is a Python package that models the spatial dependence structure of climate extremes (heatwaves, heavy precipitation) using max-stable processes. It estimates local anisotropy parameters, clusters grid cells by dependence similarity, and computes tail-risk metrics (VaR, ES) per cluster.

The package replaces the R, shell and SLURM scripts of the original workflow of Contzen et al. (2025) with a single pip-installable Python library and a CLI. R is not needed to install, run, or test the package. The tests compare against outputs that the reference R implementation produced once and that are stored in `tests/reference_data/`.

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

The CLI provides three subcommands.

| Command | Status | Actual entrypoint | Notes |
|--------|--------|-------------------|-------|
| `weatherisk validate` | Implemented | `weatherisk.pipeline.run_pipeline` | Lightweight synthetic sanity check |
| `weatherisk maps` | Implemented | `weatherisk.cpc_pipeline.run_cpc_pipeline` | Main CPC real-data pipeline with risk maps |
| `weatherisk cmip6` | Implemented | `weatherisk.cmip6_pipeline.run_cmip6_pipeline` | Figure 9 reproduction path |

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

The CMIP6 path also supports checkpoint-and-resume through `output/cmip6_fig9/checkpoints/`: full checkpoints after Steps 2–5 and incremental per-cluster checkpoints during Step 6 for both LEC and EDC re-estimation.

**Usage:**

```bash
weatherisk cmip6
weatherisk cmip6 --workers 16
weatherisk cmip6 --data-dir /pool/data/CMIP6/.../pr/gn/
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
| `density` | Pairwise composite likelihood (`pairwise_density_summand`); global and local MLE optimisation with the heuristics of the R reference (parscale, gamma-wrapping retry, boundary-proximity retry) |
| `backend` | NumPy and SciPy kernels shared by the pipelines: summed pairwise negative log-likelihood, ellipse dissimilarity wrappers, multistart L-BFGS-B optimisers |
| `estimation` | Spatial smoothing of local estimates (angular wrapping for γ); in-cluster re-estimation with average log-likelihood computation |
| `clustering` | Ellipse-overlap dissimilarity (`calc_distance_ellipses`); madogram-based Saunders method (`c_extrcoeff_matrix`, column-major flatten); hierarchical clustering; threshold-based *k* selection |
| `extremes` | Block-maxima extraction, GEV fitting, Fréchet transform |
| `cmip6_pipeline` | AWI-ESM-1-1-LR Figure 9 pipeline (`CMIP6Config`, `run_cmip6_pipeline`) with checkpoint-and-resume |
| `cmip6_data` | Discovery, download, and loading of AWI-ESM-1-1-LR monthly precipitation files |
| `risk` | Value at Risk and Expected Shortfall per cluster |
| `cpc_pipeline` | End-to-end CPC real-data pipeline: load CPC files → GEV → Fréchet → LEC/EDC clustering → risk → Cartopy maps |
| `map_plotting` | Filled-region Cartopy maps (pcolormesh): cluster maps, parameter fields, risk choropleths, summary panels |
| `gdp` | Gridded GDP exposure loading (Kummu et al. 2018), regridding to pipeline grid, cell-level extraction |
| `netcdf` | NetCDF climate data ingestion (CPC, ERA5), longitude wrapping |
| `scalable` | Library-only helpers for coarse-grid proxy clustering and label propagation |
| `parameters` | Named parameter presets (`stripes`, `bigsmall`, `rotate`) as dataclasses |
| `plotting` | Heatmaps, cluster maps, dendrograms, choropleth, bar charts (synthetic data) |
| `pipeline` | Synthetic orchestration: lightweight validation (`run_pipeline`) and full non-stationary study (`run_nonstationary_pipeline`) |
| `cli` | Click-based CLI entry point with `validate`, `maps`, and `cmip6` |
| `io` | CSV and NumPy I/O helpers; optional reading of R `.rds` files through `pyreadr` (no R installation needed) |

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

## HPC Deployment (AWI *albedo* cluster)

For large-scale publication workloads—such as the full AWI-ESM-1-1-LR T63 grid
(192 × 96 = 18,432 cells, 156 years: 1850–2005)—the pipeline should be run on
HPC infrastructure. All SLURM scripts and environment setup files are in `hpc/`.

### Data volume context

A single monthly precipitation variable on the T63 grid spans
192 × 96 × 156 × 12 = **34,504,704 values** (≈10⁷ per field). Multi-variable
workflows easily reach 10⁸–10⁹ values. Scaling to a global 0.25° grid increases
cell count by ~56×; switching from monthly to daily output adds another ~30×—over
**three orders of magnitude** more data than the CMIP6-class example here.

The pairwise-comparison stage grows **quadratically** in cell count: moving from
~18,000 cells (T63) to ~10⁶ cells (km-scale) increases the leading-order cost by
~10⁴. GPU acceleration, sparse graph-based strategies, or approximate nearest-
neighbour methods will be required for next-generation Earth system simulations.

### SLURM configuration (pure-Python runs, November 2026)

The Figure 9 runs reported in the SoftwareX paper used `hpc/run_fig9_pyonly.slurm`, which refuses to start if `Rscript` is on the search path.

| Parameter | Value |
|---|---|
| Partition | `smp` (shared memory) |
| Nodes / tasks | 1 / 1 |
| CPUs per task | 64 |
| Memory limit | 240 GiB |
| Wall-time limit | 48 h (`--qos=48h`) |
| Python | 3.11.7 |
| Thread control | `OMP/MKL/OPENBLAS_NUM_THREADS=1` |

Thread libraries are forced to single-threaded mode to prevent over-subscription
when 64 worker processes are active simultaneously.

### Running on albedo

```bash
# Environment setup (run once)
bash hpc/setup_albedo.sh

# Submit the full Figure 9 pipeline (STL detrending)
sbatch hpc/run_fig9_pyonly.slurm

# Without detrending
FIG9_EXTRA_ARGS=--no-detrend sbatch hpc/run_fig9_pyonly.slurm

# Check job status
squeue -u $USER

# Resume after an interruption: resubmit with the same output directory.
# The pipeline loads the checkpoints written there.
FIG9_OUTPUT_DIR=/path/to/previous/output sbatch hpc/run_fig9_pyonly.slurm
```

### Results of jobs 48932578 and 48932579

Both runs covered the T63 grid (18,432 cells) and 156 years on 64 cores and resumed from a checkpoint written after stage 1.

| | STL detrending (48932578) | No detrending (48932579) | Contzen et al. (2025) |
|---|---|---|---|
| Wall time (`sacct`) | 36:02:48 | 25:58:03 | |
| Job MaxRSS (`sacct`) | 95 GiB | 85 GiB | |
| Largest single process (`time -v`) | 36.4 GiB | 38.0 GiB | |
| k_LEC | 14 | 23 | 24 |
| k_EDC | 142 | 104 | 104 |

The job MaxRSS sums resident memory over all processes and counts shared pages more than once.
Stage 6 runs one worker per cluster, so the last large clusters run alone while the other cores idle.

---

## Testing

```bash
# Run tests
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=weatherisk --cov-report=term-missing
```

The suite has 235 tests: 116 test individual functions against R, 16 compare the stages of the climate pipeline with R on synthetic data, and 103 test structure and behaviour. None of the tests calls R. The R comparisons read frozen outputs in `tests/reference_data/`, which were generated once with the R functions of Contzen et al., R's `stl()` and the scripts in `tests/`. `python -m pytest tests/` runs all 235 tests, including one slow end-to-end run; continuous integration runs `pytest -m "not slow"`. Each stage is held to a criterion that an independent implementation can meet:

- **Deterministic functions** (covariance, extremal-coefficient conversions, pairwise likelihood summands, madogram, ellipse-overlap distance, LEC merge heights) must match the R outputs on identical inputs to absolute tolerances between 1e-14 and 1e-10, depending on the function (the largest difference observed is 1.4e-14). See `tests/parity/` and `tests/test_r_cross_validation.py`.
- **Local optimiser.** At every cell of the synthetic monthly precipitation generated by `tests/generate_cmip6_mini_reference.R` (5×5 grid, 20 years), the summed pairwise negative log-likelihood at the Python optimum may exceed its value at R's optimum by at most 0.03. Parameter values are not compared. With 20 annual maxima per cell the likelihood is nearly flat in `b` and γ where `a` is at its lower bound, so very different parameters reach the same likelihood. See `tests/test_cmip6_rparity.py`.
- **LEC clustering** receives R's smoothed local estimates and must reproduce R's quantile threshold and cluster count *k* exactly.
- **EDC (Saunders) clustering.** SciPy and R's `hclust` break ties between equal distances in a different order, so later merge heights may differ. Each dendrogram must be a valid average-linkage run on the same distance matrix, and the quantile threshold and cluster count must equal R's exactly.
- **Structural tests** cover simulation properties, checkpoint-and-resume, serial and parallel equivalence, and parameter presets.

To regenerate the fixtures, run `tests/generate_r_reference.R` and `tests/generate_cmip6_mini_reference.R`. This step needs R ≥ 4.0 with the `lhs` package. Running the tests does not.

## License

MIT

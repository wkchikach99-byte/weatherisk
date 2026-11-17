# Changelog

All notable changes to `weatherisk` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.2.0] — 2026-11-17

### Removed
- R is no longer called at run time. The package does not call `Rscript` for local optimisation, hierarchical clustering, GEV fitting, or STL detrending.
- The optional Rust backend (`crates/weatherisk_core`), the vendored `lbfgsb-rs-pure` crate, `hpc/setup_rust.sh`, and `tests/test_rust_parity.py`.
- R and Rust setup steps from the HPC scripts in `hpc/`.

### Changed
- Local maximum-likelihood estimation uses the SciPy L-BFGS-B multistart (Latin hypercube starts, parameter scaling, gamma wrapping, boundary retries).
- Hierarchical clustering uses `scipy.cluster.hierarchy.linkage`.
- GEV fitting and serial STL detrending use the Python implementations.
- A failed optimiser start is no longer silently ignored.
- Validation criteria. Tests read frozen outputs of the reference R implementation in `tests/reference_data/` and do not need R.
  - Deterministic computations (grid, distances, covariances, densities, madogram, extremal coefficients, smoothing, annual maxima, STL) match R within tolerances of 5e-15 to 1e-9. The ellipse dissimilarities match once R's indexing bug for pairs of empty ellipses is replayed. The inverse map from extremal coefficient to covariance is checked by evaluating the forward map at the Python root, to 1e-14, because R's `uniroot` stops at its default tolerance of about 1.2e-4. One STL series with mostly zero values is held to R's own sensitivity to a one-ulp change of the input, with trends to 1e-9 relative and robustness weights to 2e-5.
  - GEV fits match R within 0.001 in the parameters and 0.2% in the Fréchet values.
  - The local optimiser is judged by likelihood. At each cell of the mini CMIP6 reference, the summed pairwise negative log-likelihood at the Python optimum may exceed its value at R's optimum by at most 0.03, with the default start budget. Parameter values are not compared, because the start designs and L-BFGS-B implementations differ and the likelihood is nearly flat in `b` and γ where `a` is close to 0.
  - LEC on R's smoothed local estimates reproduces R's threshold to 1e-14 and R's *k* exactly.
  - EDC on the mini CMIP6 reference reproduces R's threshold to 1e-12 and R's *k* exactly. On the simulated reference, where many madogram distances are tied, merge heights may differ by up to 0.015 and *k* by 1, because SciPy and R merge tied distances in different orders.
- README, `CITATION.cff`, and the SoftwareX manuscript describe the package as independent of R at run time.

---

## [0.1.1] — 2026-05-10

### Changed
- Updated Zenodo DOI in paper metadata table to v0.1.1 archive (`10.5281/zenodo.20108563`).
- Expanded Contzen et al. (2025) reference to full citation with DOI in cover letter.
- Restored clean final v20 paper source; moved tracked-changes variant to separate file.

### Removed
- Risk-pipeline work-in-progress removed from public repository (separate research branch).
- Non-essential internal files stripped from public release.

---

## [0.1.0] — 2026-05-07

### Added
- Initial public release of `weatherisk`.
- Complete Python pipeline for spatial clustering of climate extremes via max-stable processes.
- Seasonal-trend decomposition (STL), GEV fitting, unit-Fréchet transform.
- Local pairwise composite-likelihood MLE for anisotropy parameters `(a, b, γ)`.
- Two clustering schemes: Local Estimates Clustering (LEC) and Extremal Dependence Clustering (EDC).
- In-cluster re-estimation with average log-likelihood computation.
- Three CLI commands: `weatherisk validate`, `weatherisk maps`, `weatherisk cmip6`.
- Three synthetic parameter presets: `stripes`, `bigsmall`, `rotate`.
- Full R-parity test suite validating numerical agreement with the reference R implementation.
- CMIP6 Figure 9 reproduction pipeline (AWI-ESM-1-1-LR, T63 grid, 1850–2005).
- Checkpoint-and-resume support for long HPC runs.
- `CITATION.cff` and Zenodo archive (`10.5281/zenodo.20069739`).

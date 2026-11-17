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
<!-- TODO(U4): confirm the criteria below against the final test inventory -->
- Validation criteria. Tests read frozen outputs of the reference R implementation in `tests/reference_data/` and do not need R. Deterministic functions must match R to near machine precision. The local optimiser is judged by likelihood: at each cell, the summed pairwise negative log-likelihood at the Python optimum may exceed its value at R's optimum by at most 0.03. Parameter values are not compared, because on the mini CMIP6 reference the likelihood is nearly flat in `b` and γ. LEC clustering is tested on R's smoothed local estimates and must reproduce R's threshold and *k* exactly. EDC clustering may differ from R only in how tied distances are merged.
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

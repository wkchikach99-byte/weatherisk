# Remove R from the weatherisk runtime plan

weatherisk is presented as a Python package validated against the R code of Contzen, Dickhaus and Lohmann (2025). Since commit `0a79109` (2026-03-11) four steps call `Rscript` at runtime instead of running Python. This plan removes every runtime call to R. The R code keeps one role. It is the frozen reference that the tests compare against. The rule the program enforces is that `weatherisk/` contains no call to `Rscript` and the full test suite passes on a machine without R.

Units in order. U1 removes dead Rust code. U2 restores the Python implementations. U3 picks the optimizer start budget from measurement. U4 rewrites the parity tests around criteria an independent implementation can meet. U5 adds the structural guard and CI. U6 regenerates Figure 2 and Figure 9 in pure Python. U7 updates the paper, README, metadata and release.

## How to read this

One box is one unit of work. Every box names the evidence that checks it. A nested box is a sub-step of the box above it. Check a box only when its evidence exists, a file, a log line, a figure, a test run, or a SHA. The body is a how-to. The appendices explain and record.

The program runs `playbooks/autopilot-stack.md` from the pstack poteto-mode skill, adapted to a single-author research repository. The agent builds and verifies each unit as one commit on a local branch. Widad Kchikach reviews and lands every unit. Nothing is pushed without her explicit go.

Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

## Program checklist

### Arm the program

- [ ] State the protocol and this plan to the operator, then stop. Start execution only on the operator's explicit go.
- [ ] Create branch `remove-r-runtime` from the current working tree after the operator commits or stashes the pending ruff and version fixes.
- [ ] Record the baseline. Run `PATH=/usr/bin:/bin .venv/bin/python -m pytest -q` and save the summary line. Expected today is 53 failures caused by `RuntimeError: Rscript not found on PATH`.

### Hold the boundaries

- [ ] U1 touches only `crates/`, `vendor/`, `tests/test_rust_parity.py` and Rust branches in `weatherisk/`.
- [ ] U2 touches only `weatherisk/backend.py`, `weatherisk/clustering.py`, `weatherisk/extremes.py`, `weatherisk/cmip6_pipeline.py`, `weatherisk/density.py`.
- [ ] U4 touches only `tests/`.
- [ ] U7 touches only `docs/`, `README.md`, `CHANGELOG.md`, `CITATION.cff`, `pyproject.toml`, `weatherisk/__init__.py`.
- [ ] Never edit `r_code/` or the frozen CSVs in `tests/reference_data/`. They are the oracle.

### Every unit

- [ ] Run `.venv/bin/ruff check weatherisk/` before the commit.
- [ ] Run the suite with R hidden, `PATH=/usr/bin:/bin .venv/bin/python -m pytest -q -m "not slow"`, and paste the summary line into the commit body.
- [ ] Show the operator the diff and the summary line. Land only on her go.

## Delete the dead Rust code (U1)

**Depends on.** None.

**Files.**

- [ ] Delete `tests/test_rust_parity.py`. All 40 of its tests are skipped today.
- [ ] Delete `crates/` and `vendor/lbfgsb-rs-pure/` after `rg -n "crates/|lbfgsb" --glob '!crates/**' --glob '!vendor/**'` shows no live reference.
- [ ] Edit `weatherisk/backend.py`, `weatherisk/density.py`, `weatherisk/covariance.py`, `weatherisk/grid.py` to drop every `_HAS_RUST_*` and `_USE_RUST` flag.

**Build.**

- [ ] Remove the Rust branches only. Python behaviour is unchanged.

**You see.**

- [ ] `rg -n "_HAS_RUST|_USE_RUST|_rc\." weatherisk/` prints nothing.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] The suite summary with R on PATH equals the baseline minus the 40 deleted skips. Run `.venv/bin/python -m pytest -q -m "not slow"`.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Run `weatherisk --help` and the README quickstart. Pass when both exit 0.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] None. U1 deletes unreachable code and changes no executed path.

**Review gate.** None. U1 is not review-gated.

## Restore the Python implementations (U2)

**Depends on.** U1.

**Files.**

- [ ] Edit `weatherisk/backend.py`.
- [ ] Edit `weatherisk/clustering.py`.
- [ ] Edit `weatherisk/extremes.py`.
- [ ] Edit `weatherisk/cmip6_pipeline.py`.
- [ ] Edit `weatherisk/density.py` docstrings that say the optimizer delegates to R.

**Build.**

- [ ] Replace the body of `backend.optimize_local_mle` with the pre-`0a79109` SciPy L-BFGS-B multistart (`git show 0a79109^:weatherisk/backend.py`). Keep the `lower_bounds` and `upper_bounds` parameters added later.
- [ ] Replace `clustering.clustering` with SciPy `linkage(condensed, method=method)`.
- [ ] Delete the R branch at the top of `extremes.fit_gev`. The Python MLE below it is already the fallback.
- [ ] Delete the `n_workers == 1` R branch in the CMIP6 detrending step. The serial and parallel paths then both use Python STL.
- [ ] Delete `_optimize_local_mle_via_r`, `clustering_via_r`, `_fit_gev_via_r`, `_detrend_grid_via_r`, the three `_rscript_path` helpers, the R-only cache code, and their `subprocess` and `shutil` imports.
- [ ] Keep everything else `0a79109` added in Python, such as `crop_matrix`, `calc_distance_ellipses_condensed` and the `brentq` tolerance.

**You see.**

- [ ] `rg -n "Rscript|_via_r|subprocess" weatherisk/ --glob '!cmip6_data.py' --glob '!benchmarks.py'` prints nothing. Those two files call `wget` and a child Python process, not R.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Run the suite with R hidden. The prototype in Appendix A predicts 188 passed and 7 failed. All 7 are strict optimizer or clustering parity tests that U3 and U4 resolve.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Run the mini CMIP6 pipeline end to end with R hidden. Pass when it writes `k_edc`, `k_lec` and the cluster maps without error.
- [ ] Run Listing 1 and Appendix A.1 at 11×11 grid and 30 replicates with R hidden through the existing `verify.py`. Pass when both scripts produce identical labels.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Metric. Wall time of the mini CMIP6 pipeline.
- [ ] Probe. Three interleaved runs each of the R-backed tree and U2.
- [ ] Baseline. Record the R-backed time first.
- [ ] Rule. Report the ratio with its spread. No pass threshold, because U3 sets the final cost.

**Review gate.** The operator reviews the diff before landing.

## Choose the optimizer start budget from measurement (U3)

**Depends on.** U2.

**Files.**

- [ ] Edit `weatherisk/backend.py`.
- [ ] Edit `weatherisk/cmip6_pipeline.py`, `weatherisk/cpc_pipeline.py`, `weatherisk/presets.py` or wherever `ensemble` defaults live.
- [ ] Create `scripts/optimizer_start_sweep.py`.

**Build.**

- [ ] Write `scripts/optimizer_start_sweep.py`. For each configuration it reports, per cell of the mini CMIP6 reference, `NLL(Python optimum) − NLL(R optimum)` and the wall time. It extends `nll_compare.py` from Appendix A.
- [ ] Sweep `ensemble` ∈ {5, 10, 15, 20} crossed with the start design ∈ {plain Latin hypercube, `qmc.LatinHypercube(optimization="random-cd")`}.
- [ ] Set the default to the cheapest configuration where no cell is worse than R by more than 0.03. Use one source of truth for the default instead of repeating the number in each config class.

**You see.**

- [ ] The sweep table in the commit body, with a row for each configuration showing cells worse than R, cells better than R, and wall time.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] `test_local_mle_matches_r[python]` passes. Run `PATH=/usr/bin:/bin .venv/bin/python -m pytest -q "tests/test_cmip6_rparity.py::TestMiniCMIP6Pipeline::test_local_mle_matches_r"`.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Run the stripes simulation preset with known true parameters at 21×21. Pass when the recovered semi-major axis tracks the true `b(x, y) = (x + 5)/2` at least as well as the pre-U3 default, measured as RMSE.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Metric. Seconds per local-MLE cell on the mini reference and on a 16×16 CMIP6 subgrid.
- [ ] Probe. Interleave the old default and the chosen default, three runs each.
- [ ] Baseline. Record the 5-start time first.
- [ ] Rule. Report the measured slowdown. Extrapolate the Figure 9 HPC cost from it before U6 and flag it to the operator if the projection exceeds the albedo wall-time limit.

**Review gate.** The operator reviews the sweep table before landing.

## Rewrite the parity tests around criteria an independent implementation can meet (U4)

**Depends on.** U3.

**Files.**

- [ ] Edit `tests/test_cmip6_rparity.py`.
- [ ] Edit `tests/parity/test_parity_pairwise_density_optim_local.py`.
- [ ] Edit `tests/parity/test_parity_calc_estimates_in_clusters.py`.
- [ ] Edit `tests/parity/test_parity_clustering.py`.
- [ ] Edit `tests/parity/PARITY_REPORT.md` and `docs/python_r_parity_migration_plan.md`.
- [ ] Create `tests/test_no_r_runtime.py`.

**Build.**

- [ ] Keep strict tolerances for every deterministic function that already meets them, such as likelihood summands, covariance, madogram, ellipse distance and LEC merge heights. These are the comparisons that show the Python code computes the same quantities as R.
- [ ] Replace the optimizer parameter-equality assertions with a likelihood criterion. At every cell the Python optimum must reach a pairwise log-likelihood at least as high as R's selected optimum, within 0.03. Appendix B explains why parameter equality is not a valid criterion here.
- [ ] Test LEC clustering on R's own local estimates. Feed `local_estimates_smoothed.csv` from the R reference into Python's ellipse distance and clustering, and require R's `q30` and `k_lec` exactly. This isolates the clustering stage from optimizer non-identifiability.
- [ ] Handle Saunders and EDC ties explicitly. Measure whether SciPy and R `hclust` give the same partition at the selected k on the reference madogram. If they match, assert the partition. If they differ only through tied distances, assert that each SciPy merge height equals an R merge height up to tie order and document the difference.
- [ ] Drop the `backend` parametrisation over `python` and `rust`.
- [ ] Every test must pass with R hidden. Tests that need R to regenerate fixtures stay in `tests/generate_*.R` and are never collected by pytest.

**You see.**

- [ ] `PATH=/usr/bin:/bin .venv/bin/python -m pytest -q` reports 0 failed.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Break the behaviour on purpose and confirm the new tests catch it. Set `ensemble=1` and confirm the likelihood test fails. Perturb one smoothed estimate by 0.5 and confirm the LEC-on-R-estimates test fails. Revert both.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Rebuild the paper's validation table counts from `pytest --co -q` output by test file. Pass when every number in the table is generated by a script, not typed.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Metric. Wall time of `pytest -m "not slow"`. Rule. Stay under 10 minutes so CI on GitHub runners remains practical.

**Review gate.** The operator and, ideally, T. Dickhaus review the new acceptance criteria before landing, because they become the paper's validation claim.

## Enforce no R at runtime in structure (U5)

**Depends on.** U4.

**Files.**

- [ ] Create `tests/test_no_r_runtime.py`.
- [ ] Edit `.github/workflows/tests.yml`.

**Build.**

- [ ] `tests/test_no_r_runtime.py` scans every file under `weatherisk/` and fails if it finds `Rscript`.
- [ ] CI already installs no R. Keep it that way, so every CI run proves the package works without R.

**You see.**

- [ ] CI is green on GitHub after the operator pushes.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Add an `Rscript` string to a scratch module, confirm the guard fails, remove it.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] In a clean virtual environment, `pip install .` from a fresh clone and run Listing 1 at small size. Pass when it runs on a machine where `which Rscript` prints nothing.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] None. U5 adds a source scan and changes no runtime path.

**Review gate.** None. U5 is not review-gated.

## Regenerate the paper's results in pure Python (U6)

**Depends on.** U5.

**Files.**

- [ ] Create new outputs under `output/` and `docs/software_x/final/figures/` without overwriting the R-backed files, so the two can be compared.

**Build.**

- [ ] Figure 2. Run the `paper_stripes` preset, 51×51 grid, 250 replicates, seed 42, on albedo together with Figure 9. Record `k_edc`, `k_lec` and the run time. Full-scale runs never run on the workstation. Local checks before this unit stay at 11×11 or smaller, with n_workers ≤ 2 and a 10-minute timeout per command.
- [ ] Figure 9. Run `hpc/run_fig9_clean.slurm` on albedo with the U3 defaults. Record `k_edc`, `k_lec` and the run time.
- [ ] Compare each figure with the R-backed version and with Contzen et al. (2025). Report adjusted Rand index between partitions and the change in k.

**Figure criteria, fixed on 2026-11-07 before any pure-Python full run.** Each criterion is either asserted or reported. Only the deterministic steps are asserted. An asserted criterion that fails is a bug to find. A reported number is never tuned toward a target.

- F1, asserted. The EDC chain contains no random step: STL, block maxima, GEV fit, Fréchet transform, madogram, `hclust` and the selection rule. On the same full Figure 9 input, Thorsten's R code (run offline on albedo for validation only) and weatherisk give the same `k_edc`. The partitions must have ARI = 1, or differ only where tied distances are merged in a different order, and that must be shown.
- F2, asserted. The LEC clustering stage is deterministic once the smoothed estimates are given. R's ellipse distance and clustering, run on weatherisk's full-data smoothed estimates, give weatherisk's `q30` and `k_lec` exactly. The indexing bug in R's `calc_distance_ellipses` is handled the way the paper describes.
- F3, reported. Local MLE uses random starts on a flat likelihood, so `k_lec` is a random quantity. Report `k_lec` and the ARI between partitions over several optimizer seeds if HPC time allows, otherwise for one seed with that limitation stated. Contzen et al.'s `k_lec = 24` is reported next to it, not used as a target.
- F4, reported. Compare the `k` values and maps with Contzen et al. (2025). Describe every difference and attribute it only to a cause that has been tested, otherwise call it untested.
- F5, reported. Figure 2 (stripes) has known true parameters. Report the recovered `k` and the ARI against the true stripe partition, and compare them with the earlier version of the figure.

**You see.**

- [ ] A short comparison note, `docs/software_x/final/python_vs_r_results.md`, with both k values, the ARI and side-by-side maps.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] `tests/test_paper_reproduction.py` passes with R hidden.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Open both new figures and confirm visually that the cluster structure is spatially coherent and comparable to the published one.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Metric. Figure 9 wall time on 16 workers. Baseline. The earlier R-backed run, roughly 4 hours according to earlier notes, not re-measured. Rule. Report the measured time in the paper's performance paragraph.

**Review gate.** The operator and both supervisors review the new figures before the paper changes. If k changes, the text around Figure 2 and Figure 9 changes with it.

## Update the paper and release (U7)

**Depends on.** U6.

**Files.**

- [ ] Edit `docs/software_x/final/softwarex_paper_v23.tex`.
- [ ] Edit `README.md`, `CHANGELOG.md`, `CITATION.cff`, `pyproject.toml`, `weatherisk/__init__.py`.
- [ ] Edit the cover letter.

**Build.**

- [ ] Metadata C6 reads "Python 3.10+; SLURM (optional, HPC runs)". Remove "optional Rust backend" and "R (validation)" from runtime requirements. Say instead that validation fixtures were generated once from the reference R implementation.
- [ ] Rewrite the validation section and both tables from the U4 test inventory, including the likelihood criterion and why parameter-level equality is not used.
- [ ] Remove Rust and PyO3 from the text, the Conclusions and the bibliography.
- [ ] Replace Listing 1 and Appendix A.1 with the scripts already verified in this session.
- [ ] Fill funding. W. Kchikach by the MarDATA Helmholtz School for Marine Data Science (HIDSS-0005). G. Lohmann by ERC Synergy Grant "i2B" (101118519). T. Dickhaus only if he supplies a grant.
- [ ] Fix uncited references, the doubled "Appendix" heading and the doubled periods.
- [ ] Bump to v0.1.2, tag it, and let Zenodo mint a new DOI. Put the new DOI in the paper and the cover letter.

**You see.**

- [ ] `tectonic softwarex_paper_v23.tex` builds with no undefined references.

**Verify, unit.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Every code listing in the PDF is extracted and executed by a script. Pass when all exit 0.

**Verify, live.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] Read the compiled PDF end to end against the repository at the v0.1.2 tag.

**Verify, perf.** Tests alone are not sufficient verification. A unit is verified only when its unit, live, and perf boxes are all checked.

- [ ] None. U7 changes text only.

**Review gate.** Both supervisors approve the final PDF before submission.

## Close the program

- [ ] Every box above is checked with its evidence.
- [ ] Report to the operator. Include the final test summary with R hidden, the k values before and after, and the new DOI.

## Appendix A. Prototype evidence

All runs used a throwaway copy of the working tree, including uncommitted changes, at `/private/tmp/claude-782984363/-Users-wikchi001--claude-spares-spare-yRrH4o/892b6b31-8b03-40e3-bd5c-5dfa47b387cb/scratchpad/proto`, with R removed from `PATH`. The real repository was not modified.

- Current code with R hidden. 53 tests fail with `Rscript not found on PATH`.
- U2 applied, 5 optimizer starts. 188 passed, 7 failed, 40 skipped. All 7 failures are strict optimizer or clustering parity tests.
- Optimizer starts. With 5 starts Python's optimum was worse than R's by more than 0.03 in negative log-likelihood at 4 of 25 mini CMIP6 cells, up to 0.79. With 20 and with 50 starts no cell was worse.
- Likelihood comparison at 20 starts, script `nll_compare.py` in the scratchpad. Python is better than R at 5 cells, ties at 20 and is worse at 0.
- Parameters at tied likelihoods differ a lot. For example cell 1 has Python (0.01, 0.01, −0.04) against R (0.01, 6.67, 0.41) with the same likelihood. The mini reference has 20 annual maxima per cell, so with `a` at its lower bound the likelihood is flat in `b` and γ.
- Downstream effect on the mini grid. LEC gives k = 7 in Python against k = 11 in R. Smoothed `b` differs by up to 3.07. Both follow from the flat likelihood.
- Clustering. LEC merge heights already match R to better than 5e-14. Saunders merge heights differ by up to 3.2e-3, which `tests/parity/test_parity_clustering.py` itself attributes to tie-breaking between R `hclust` and SciPy `linkage`.

Still unproven. Whether k changes on the full Figure 2 and Figure 9 data, where identifiability is better. U6 answers that.

## Appendix B. Alternatives rejected

- **Keep R and disclose it.** Rejected by the operator. It contradicts the paper's purpose, a Python implementation.
- **Port R's `optim` L-BFGS-B and `maximinLHS` bit for bit into Python.** This would chase the March 2026 contract of `< 1e-14` agreement on optimizer outputs, from `docs/python_r_parity_migration_plan.md`. Appendix A shows that contract is the wrong premise. At tied likelihoods the R parameters are one arbitrary point on a flat ridge, chosen by R's random starts. Matching them proves agreement with R's random number generator, not correctness of the estimator. Every independent implementation would fail that test.
- **Keep R as an optional parity mode.** Rejected under the Laziness Protocol. It keeps four subprocess paths and a second code path that CI never runs. The frozen CSVs already preserve R's outputs for comparison.

## Appendix C. Risks

- **k changes in Figure 2 or Figure 9.** Lands in U6. If it does, the paper's numbers change, and the supervisors must agree with the new text.
- **More optimizer starts make Figure 9 slower.** Lands in U3. The owner measures the cost on a subgrid before booking HPC time.
- **SciPy and R disagree on tied EDC distances.** Lands in U4. The owner measures partition agreement on the reference data before choosing the assertion.
- **The validation claim changes.** Lands in U4 and U7. The supervisors must accept the likelihood criterion, because a reviewer will read it as the method's validation.

## Appendix D. Links and reading list

- `git show 0a79109` for the change that introduced the R calls.
- `docs/python_r_parity_migration_plan.md` for the March 2026 strict-parity contract this plan replaces.
- `tests/parity/PARITY_REPORT.md` for the per-function parity status.
- `tests/generate_r_reference.R` and `tests/generate_cmip6_mini_reference.R` for how the frozen fixtures were made.
- Contzen, Dickhaus and Lohmann (2025), Extremes, https://doi.org/10.1007/s10687-025-00519-2.

# Judgment calls on branch `remove-r-runtime` (for Widad and the supervisors)

This list covers what was decided during the unattended run of 2026-11-07 without asking.
Each item says what was decided, why, and what you could tell me to do instead.
Items marked **Supervisors** change a scientific claim in the paper and should be confirmed by T. Dickhaus and G. Lohmann.
The full decision trail with timestamps is `.audit/remove-r-runtime.tsv` (local, not committed).

## A. Validation criteria against the reference R code of Contzen et al. in r_code/ (Supervisors)

1. **Deterministic functions are compared strictly.**
   STL, GEV transform, madogram, ellipse distance, smoothing and the k rule must match R to about 1e-13 to 1e-14 relative.
   Some fixtures are 15-digit CSV files, so those tests use rtol 5e-15 instead of exact equality.
   The paper still says "≤ 1e-10" for Level 1 until the final test inventory measures the actual worst case.

2. **The local optimizer is judged by likelihood, not by parameter equality.**
   Criterion is Python NLL ≤ R NLL + 0.03 per cell, with NLL the summed pairwise negative log-likelihood of that cell.
   Reason. The mini-data likelihood is flat. In 14 of 25 mini cells the range parameter `a` sits on its lower bound 0.01.
   Two optimizers can reach different parameters at practically the same likelihood, so parameter equality is not a well-posed test.
   The tolerance 0.03 was chosen by us. It has no external source.
   Alternative. A tolerance relative to the NLL, or a likelihood-ratio style bound.

3. **Optimizer start count is larger than R's.**
   R uses 5 random starts (`locest_ensemble = 5`). Python uses SciPy L-BFGS-B from Latin-hypercube starts.
   With 5 starts Python was worse than R by more than 0.03 at some cells. The sweep found that 50 plain LHS starts is the smallest tested budget that meets the criterion. Re-measured on 2026-11-09 over 30 cells (25 mini cells, 5 cells of the 10x10 field): 50 starts pass for all 21 seeds 0–9, 42 and 100–109, worst excess +0.0287 (seeds 4 and 103, cell mini:16), only 0.0013 inside the tolerance. The recurring value suggests a secondary local optimum at that cell. Of seeds 0–9, 6 fail with 20 starts (worst +1.13) and 4 with 40 (worst +0.79); 1 of seeds 100–109 and seed 42 fail with 40. The first sweep's "20 fails 5/10" used a seed set that was not recorded and is not reproduced; the paper now gives the re-measured numbers.
   Plain LHS designs are not nested, so 40 starts can fail where 20 pass. The nested design was not adopted.
   This is an optimizer setting, not a model change. It increases run time roughly in proportion.
   Alternative. Keep 5 starts and report the likelihood gap instead of asserting the criterion.

4. **LEC is tested on R's smoothed estimates (and, on HPC, R LEC is run on Python's estimates).**
   Because local estimates differ between optimizers within the likelihood tolerance, the LEC clustering step is tested in isolation.
   Given the same smoothed estimates, Python must reproduce R's q30 threshold and k exactly.

5. **EDC/Saunders clustering may differ only by tie-breaking.**
   Average linkage on a matrix with tied dissimilarities has no unique merge order. On the Saunders data, changing only the cell order moved R's own k from 12 to 9.
   Tests therefore assert agreement up to ties, and report where ties occur.
   Full grid (F1, job 48931982). The 169,860,096 madogram values take only 5,134 distinct values.
   On R's exact doubles SciPy gives k = 142 and R gives k = 140.
   A NumPy transcription of R's Fortran hclust (`hpc/r_reference/linkage_tie_check.py`) reproduces R's partition and k = 140 exactly, with heights within 8.3e-17.
   Permuting only the cell order moves k to 139, 141 in the transcription and to 139, 141, 142 in SciPy.
   So 140 vs 142 is tie-breaking, and F1 is met up to ties.
   weatherisk keeps SciPy. Shipping the R-order linkage would make k match R in R's cell order only. It would not make k better defined.

6. **Figure criteria were fixed before the full runs (F1 to F5).**
   * F1 (asserted). The EDC chain gives the same k in R and Python on the same full input, up to shown ties.
   * F2 (asserted). R LEC on Python's full-data estimates, with the bug fix of item B1, gives exactly Python's q30 and k.
   * F3 (reported). k_lec over several seeds. Not asserted, because the optimizer and simulation are random.
   * F4 (reported). Comparison with Contzen et al. (k_edc 104, k_lec 24). Not asserted.
   * F5 (reported). Stripes recovery. See item A7.
   Reason. Exact figure equality with the published paper is impossible where the pipeline is stochastic. We only assert what is deterministic.

7. **F5 metric was substituted.**
   The plan said "adjusted Rand index against the true stripe partition".
   The truth b = (x+5)/2 is continuous, so there is no true partition and the ARI is undefined.
   Replacement. (i) between-cluster share of the variance of true b, 1 − SS_within/SS_total, and (ii) RMSE of the in-cluster estimate of b against true b per cell (`scripts/stripes_metrics.py`).
   Alternative. Discretise b into a chosen number of stripes and compute ARI against that.

## B. Defects found in the R reference code or tooling (Supervisors for B1)

1. **Indexing bug in the reference R code's `calc_distance_ellipses` (`r_code/functions.R` around line 515).**
   The empty-ellipse branch writes `dist_matrix[j] = 1` instead of `dist_matrix[i, j] = 1`.
   On the mini data 44 entries of the distance matrix differ.
   Python implements the corrected version. F2 compares against a patched copy of R, created at run time by `hpc/r_reference/validation_patch_calc_distance_ellipses_fixed.R`. `r_code/` is untouched.
   The paper describes this. Please confirm with J. Contzen or the supervisors that `[i, j]` is the intended code.

2. **statsmodels STL is not R's `stl(s.window="periodic", robust=TRUE)`.**
   With R's parameters it differs from R by up to 0.12 over the 25 mini cells. The difference appears only when robustness weights and R's `jump` interpolation combine.
   A NumPy port of the netlib STL algorithm (as used by R) replaces statsmodels in the EDC chain.
   License note. The port follows the published Fortran algorithm of Cleveland et al. It does not copy GPL R source.

3. **R's `read.csv` on arm64 macOS is not correctly rounded.**
   Some frozen hclust references in `tests/reference_data/` were produced in R from values that `read.csv` parsed with a 1-ulp error.
   Exact comparisons now use binary files or `%.17e` strings. The frozen CSVs were not edited.

4. **Released v0.1.1 cannot be installed from source.**
   In `pyproject.toml`, `[project.urls]` came before `classifiers` and `dependencies`, so both were parsed into `project.urls` and the build failed.
   Fixed in commit 9e0bcaf. The tag v0.1.1 itself still has the defect. A v0.1.2 release is needed.

## C. Paper and release items that need you

Updated 2026-11-09. Lines marked `% RELEASE BLOCKER` or `% AUTHOR CHECK` in v23.tex and cover_letter.tex are the open items. A search for those two strings lists all of them.

1. **Version bump and Zenodo DOI.** Metadata C1 says v0.2.0, but `__version__` and `pyproject.toml` still say 0.1.1 and C3 still has the v0.1.1 DOI. Tagging v0.2.0, the Zenodo release and the push are yours to do. The cover letter says v0.2.0 with the same blocker.
2. **Example 1 (stripes) is stale.** Its k values, figure and text come from a run made before R was removed. The rerun is albedo job 48947462 (pure Python, seed 42, output `weatherisk_pyonly/output/fig2_stripes_seed42_48947462/`). It was still running at 27 min elapsed on 2026-11-09. Replace the figure and numbers when it finishes. The embedded title "Fig. 3 — EDC vs LEC" in the current image is also stale.
3. **Acknowledgement of J. Contzen (AUTHOR CHECK).** The paper calls `r_code/` "publicly available". Confirm that it is, and confirm the "helpful discussions", before submission.
4. **T. Dickhaus funding.** Not known. He is omitted from the funding sentence.
5. **Cover letter (AUTHOR CHECK).** It states that all authors have read the final version. That must be true on the day of submission.
6. **Module count.** The paper says 22 modules and a CLI entry point. `git ls-files` agrees (23 files under `weatherisk/` besides `__init__.py`, one of them `cli.py`). `koch_risk.py` and `risk_pipeline.py` exist locally but are not tracked.
7. **Test counts.** The committed tree has 235 tests (Level 1: 116, Level 2: 16, Level 3: 103). An earlier count of 238 included `tests/test_risk_pipeline.py` (3 tests), which is gitignored since commit 683d167 and is not in the public repository. With R removed from the search path, `pytest` runs every test (there is no `addopts`) and they pass in about 134 s on an Apple M4 laptop (measured 2026-11-09 by two independent runs under a load of about 3). `pytest -m "not slow"`, which CI runs, skips one end-to-end test that takes 1.5 s. An earlier note said 67 s; that figure was not reproduced, so the paper now says "about 2 min". Recount if tests change before the tag.
8. **Paper structure.** `docs/software_x/AGENTS.md` describes four illustrations (transformation chain, Fig 3, LEC map, GDP-weighted LEC map). The paper keeps its two examples. Not changed.
9. **Motivation section.** Left as Thorsten corrected it.
10. **CI.** Whether GitHub CI has R available was not checked. The R-free suite does not need it.

## D. Process notes (for transparency)

1. The ruff and version fixes made before the plan landed inside the U1+U2 commit, named in its message.
2. The first plan commit accidentally included staged Rust deletions. It was soft-reset (local, unpushed) and redone as two commits.
3. I ran `git stash` / `git stash pop` once while a subagent edited `weatherisk/`. The diff was checked afterwards and nothing was lost.
4. `tests/test_r_cross_validation.py` was left unchanged.
5. The network dropped during the night. The STL agent's last write was at about 05:20 UTC, after which it died on a DNS error (`ENOTFOUND`). At 10:02 UTC albedo did not resolve either, so the VPN was down and no HPC job could be submitted or collected.

## E. Why our k_EDC (140 to 142) differs from the published 104 (Supervisors)

1. **Which code and paper.** The 104 is from Contzen et al., Extremes 28:713–737 (2025), Sect. 5 and the Fig. 9 caption.
   `r_code/main.pdf` is byte-identical to the June 2023 preprint (`reference_papers/regionalization-approaches-for-the-spatial-analysis.pdf`), which has no climate application.
   So `r_code/` is the simulation-study code. It has no script for the 104.
   A grep of `r_code/` for STL, detrending, GEV fitting, NetCDF reading and precipitation finds none of them.
   The STL settings in weatherisk (`s.window = "periodic"`, `robust = TRUE`, subtract the trend) come from our own R scripts of March 2026 (`tests/generate_cmip6_mini_reference.R`, `scripts/replay_fig9_benchmark_in_r.R`), not from J. Contzen.
   His only public code is the Zenodo record 10.5281/zenodo.7413932 for Contzen, Dickhaus, Lohmann, PLOS ONE 2023. Its R files (`main.R`, `functions.R`, `import_*.R`) contain no STL or detrending call; the only grep matches are "coastlines". That paper states that it makes no seasonality adjustment. It is a different application (daily temperature, time-dependent GEV), so this does not show what was done for Sect. 5.
   Our R number, 140, is Contzen et al.'s unchanged madogram, hclust and k functions applied to our own preprocessing (`hpc/r_reference/r_edc_fullgrid.R`).

2. **What the paper leaves open.** The STL settings, whether the seasonal part is removed, the ensemble member and the exact data product.
   The marginal transform cannot matter, because `c_extrcoeff_matrix` uses only ranks.
   Only r1i1p1f1 exists on the albedo CMIP6 pool.

3. **Sweep (job 48932026, `hpc/r_reference/edc_preprocessing_sweep.py`).** Same data, years 1850–2005, 30 % quantile. "R order" is the transcription of R's hclust.

   | Preprocessing before annual maxima | k SciPy | k R order |
   |---|---|---|
   | STL periodic, robust, minus trend (our choice) | 142 | 140 |
   | STL periodic, robust, minus trend and seasonal | 142 | 143 |
   | STL periodic, not robust, minus trend | 144 | 147 |
   | STL s.window 7, robust, minus trend | 138 | 140 |
   | STL s.window 13, robust, minus trend | 146 | 142 |
   | Our choice, years 1850–2014 | 138 | 137 |
   | **No detrending** | **104** | **106** |

   Correction. The no-detrending row was first computed with tied maxima ranked by SciPy `method="max"`, not the reference code's rule (equivalent to `method="min"`, item E6). The recomputation (job 48932523) gives the same k under all three rules: 104 (SciPy) and 106 (R order) for `"min"`, `"average"` and `"max"`, with identical q30 0.160726 and identical k(p) curves. The ties are small. Of the 231 tied cells, the median loses 1 year to ties and the worst loses 6; no tied value is zero; the rank of any year changes by at most 4 of 156 between rules (measured). The detrended rows have no ties.

   Tie-breaking alone moves k by about 2 (item A5), so 104 and 106 are the same result.
   k is steep in the quantile level. For our choice, k = 104 needs a level of about 0.34 instead of 0.30.

4. **Reading.** The published 104 is reproduced when the annual maxima are taken from the monthly series without detrending, and by none of the six detrended variants.
   This is consistent with the published EDC result having been computed without the STL step, or with maxima taken from the raw series. It does not prove it.
   Only J. Contzen's script for Sect. 5 can settle it.
   A further check is the LEC count on the non-detrended input (the paper reports 24).
   Both full Fig. 9 runs were submitted on 2026-11-07 from commit 7561a6e, 50 local-MLE starts, 16 CPUs, qos 48h: job 48932143 with STL, job 48932144 with `--no-detrend`. Measured on the step-2 checkpoint (job 48932524), step 3 with 50 starts costs 110 s per cell (101 to 116 s over 6 cells), about 35 h on 16 workers, and step 3 writes its checkpoint only at the end. Both jobs were cancelled before step 3 could finish and resubmitted from the step-2 checkpoint with 64 CPUs and 240 GB: job 48932578 with STL, job 48932579 with `--no-detrend`. Step 3 took about 11.9 h on 64 workers (step-3 checkpoints written 2026-11-08 01:36 and 01:39 CET), longer than the 9 h guess from linear scaling.
   Step-5 results (measured, logs `fig9_pyonly_48932578.out` and `_48932579.out`):

   | Run | k_LEC | k_EDC |
   |---|---|---|
   | STL detrending (48932578) | 14 | 142 |
   | No detrending (48932579) | 23 | 104 |
   | Published (Contzen et al. 2025) | 24 | 104 |

   Without detrending both counts are at or next to the published ones. With STL, k_LEC is 14, far from 24. k_LEC depends on the optimizer through the local estimates, so 23 versus 24 is within what we would expect from a different optimizer (inference, not tested). The March 2026 STL run with the old R backend gave k_LEC 21 and k_EDC 142. Both jobs completed on 64 CPUs (measured, `sacct`): STL 36:02:48 wall time, no-detrend 25:58:03. Memory needs two numbers. `sacct` MaxRSS for the whole job is 95 GiB (STL) and 85 GiB (no-detrend). The largest single process, from `/usr/bin/time -v`, is 38,178,196 KiB = 36.4 GiB (STL) and 39,809,608 KiB = 38.0 GiB (no-detrend). An earlier version of this note called the per-process figure the job's max RSS, which was wrong. The paper's "approximately 4 hours" was removed on 2026-11-09 and Table 2 now gives the measured times. Step 6 runs one worker per cluster, LEC clusters first, so the last LEC cluster to finish (3247 cells in the STL run, the fourth largest there, and 6230 cells in the no-detrend run) ran alone for hours while the other cores idled (measured with `ps` inside the jobs on 2026-11-08; average CPU use 2214 % and 3016 % of 6400 %). The cluster counts in the final `pipeline_results.npz` equal the step-5 counts above.

   Agreement with the March 2026 STL run (R backend, `output/cmip6_fig9/pipeline_results.npz`), measured locally:

   | Comparison | Value |
   |---|---|
   | Adjusted Rand index, LEC labels, STL new vs March | 0.559 |
   | Adjusted Rand index, EDC labels, STL new vs March | 0.721 |
   | Adjusted Rand index, LEC labels, STL vs no-detrend (new) | 0.341 |
   | Adjusted Rand index, EDC labels, STL vs no-detrend (new) | 0.534 |
   | Median abs. difference of local estimates (a, b, gamma), STL new vs March | 0.008, 0.013, 0.001 |
   | Cells whose annual-maxima year ranks differ, STL new vs March | 18230 of 18432 (median max shift 3 ranks) |
   | Median relative difference of Frechet margins, STL new vs March | 2.1 % (99th percentile 19.7 %) |

   The March run detrended with a different STL implementation than the current NumPy port of R's `stl()`. The year-rank changes show that the detrended series differ in almost every cell, which explains why the March EDC labels are not reproduced exactly (inference). The cause has not been traced further; the full-grid check of the STL port against R's `stl_trend.f64le.bin.gz` is still open.
   Both jobs loaded the code from before the switch to min ranks (item E6). `hpc/r_reference/edc_labels_rule_check.py` recomputed the EDC labels from each run's step-2 checkpoint under both rules. Average ranks reproduce the step-5 checkpoint labels exactly, and min ranks give identical labels in both runs (0 cells differ). So the runs need no restart.

5. **What the SoftwareX paper says (changed 2026-11-09).** The untested sentence "we attribute this primarily to differences in preprocessing" is gone. Widad decided on 2026-11-09 to make the paper edits without waiting for supervisor review and to share this document afterwards. Correction to the proposal that stood here: Contzen et al. do state that the data were detrended with STL. They do not report the STL settings. The paper now says exactly that (v23.tex, Example 2 discussion), reports both runs as equal rows, and makes no claim about which preprocessing produced the published counts. It also gives the sensitivity results of the sweep (137 to 147 across STL variants, 140 with R's merge order, 139 to 142 under tie permutations, and a quantile level of about 0.34 needed for 104 with STL). **Supervisors.** Please check that this description of your 2025 paper is fair.

6. **Rank rule for tied maxima.** `c_extrcoeff_matrix` ranks each year as `sum(x[k] <= x)`, the number of years with a maximum at least as large. The madogram uses only rank differences, so this is equivalent to SciPy `rankdata(method="min")`. A check on tied integer data matched R bitwise with `"min"` and not with `"max"` or `"average"`.
   weatherisk used average ranks in four copies of the ranking and `method="max"` in `clustering.c_extrcoeff_matrix`. The rules differ only when a cell has tied annual maxima. That happens in 231 cells without detrending and in none with STL, so the STL results and the F1 check are unaffected.
   Effect on the full grid (job 48932523). Without detrending, k is 104 (SciPy) and 106 (R order) under `"min"`, `"average"` and `"max"`. The rule does not change k here. Decision (Widad, 2026-11-07). Follow the reference code. All five rank sites now call `clustering.madogram_ranks`, which computes `sum(x[k] <= x)` exactly. `tests/test_clustering.py::TestMadogramTies` compares `c_extrcoeff_matrix`, the CMIP6 EDC and the CPC EDC with a line-by-line transcription of the R loop on tied data; three of its tests fail on the previous code. No number reported so far changes.

7. `docs/maternity_leave_progress_update_june2026.md:37` lists "104 EDC / 24 LEC" as our reproduction. That was never our result. The file is your June report, so it was left unchanged.

## F. Paper edits of 2026-11-09 and the defaults behind them

Made without supervisor review at Widad's request. Each default can be reversed.

1. **Validation section rewritten.** Three levels with counts, criteria and the largest measured deviations (tables `tab:validation` and `tab:pipeline-parity`). Tolerances in the text are the ones the tests assert after commit 1acca38.
2. **Optimizer claim.** The paper states the criterion (NLL at most R's + 0.03 per cell) and the measurement (lower at 12 of 25 cells, equal within 1e-6 at 13, never higher by more than 7.3e-7). It does not say "at least as good as R" without the tolerance.
3. **Figure 9.** Both runs appear as a 2x2 figure and a table of counts and run times. The comparison with the March run (ARI 0.559 / 0.721) is not in the paper, because that run used a different STL implementation whose cause is untraced (E4).
4. **Speed claims.** The comparison with R speed, the 21x21 timing and the "4 hours" were removed. The paper gives the 6x6 laptop time (median 16.25 s) and the measured full-grid times.
5. **Checkpoint claim.** The paper says a checkpoint is written after each step and after each cluster within the final step. It no longer says the runs fit batch-queue limits, because step 3 has no partial checkpoint and took about 12 h.
6. **Abstract and conclusions.** Detrending is described as optional; both runs' counts are given.

## G. Second check of 2026-11-09 (all corrections made in the commit that adds this section)

Each item was measured or read in the repository on 2026-11-09. "Paper" is `softwarex_paper_v23.tex`.

1. **Test count 238 → 235.** See C7. Changed in the contributions list, Level 3 heading and text (106 → 103, other 53 → 50), Table `tab:validation`, Impact and Conclusions, and the README.
2. **"Serial vs. parallel (5)".** Only 2 of the 5 tests in `test_cmip6_performance_equivalence.py` compare serial with parallel runs. The other 3 compare a fast implementation with a reference loop. The table row and text now say both; the tolerance is 1e-10 or tighter, not "identical".
3. **STL "14 series".** `tests/reference_data/stl/cases.csv` has 14 configurations holding 78 series; the near-zero configuration (`dry240`) holds 15 series. Text corrected.
4. **STL row of `tab:pipeline-parity`.** The criterion is relative to max|x| but the deviation was given in absolute units (5.8e-13). Now 2.3e-14 max|x|.
5. **Start-count numbers.** See A3. The paper now names the 21 seeds and the 30 cells and gives 6 of 10 (20 starts) and 4 of 10 (40 starts).
6. **Vacuous test.** `test_r_parity_fixes.py::TestGammaWrappingRetry::test_gamma_wrap_triggers_at_boundary` ended in `assert ... or True`, which cannot fail. The `or True` is removed and the test passes. It still does not reach the boundary (gamma = -0.733 at that cell), so it does not exercise the wrap; the other gamma test does. Writing a case that hits the bound is open.
7. **Benchmark command.** `scripts/run_fig9_benchmark.py` with default arguments crashes in plotting on the coarse 6x6 grid (cartopy 0.25.0, shapely 2.1.2: `'GeometryCollection' object is not subscriptable` in cartopy's antimeridian wrap of a 60-degree cell centred on 180°). The global T63 maps plot without error with the same libraries (checked with `scripts/plot_fig9_comparison.py`). The paper's timings are without plotting, so the paper now gives the command with `--skip-plots` and the caption says "without plotting". Fixing the plot for coarse grids is open.
8. **Level 1 composition.** Of the 116 tests, 90 read stored R output. The other 26 (22 in `test_r_parity_fixes.py`, 4 madogram-tie tests) check translation details and tie handling against formulas or a transcription of the R loop. The paper text describes them that way; the table row label "Individual functions" covers both.
9. **Level 2 composition.** 10 of the 16 tests compare with stored R output on the synthetic data. 3 compare the pipeline's EDC functions with the function-level implementation, 1 checks pipeline vs. entry point, 1 checks the configuration, and 1 checks that the detrending is closer to STL than to a running mean. The table footnotes say so. In one of these tests R's k_saunders (12) on the 10x10 field is loaded but not asserted; Python gives 13 because of tie order.
10. **Synthetic, not AWI.** The Level 2 data are generated by `tests/generate_cmip6_mini_reference.R` (seasonal sine, weak trend, Gaussian noise, seed 42), not a subset of AWI-ESM-1-1-LR. Paper and README now say synthetic.
11. **Peak memory units.** `time -v` reports KiB: 38,178,196 and 39,809,608 KiB = 36.4 and 38.0 GiB, as the paper says. An earlier note in E4 wrote "38.2 GB / 39.8 GB"; corrected.
12. **Cover letter.** "seven-stage" → "six-stage"; trend removal marked optional; tolerance sentence matches the paper; "In an earlier study with J. Contzen, two of us …" makes the co-authorship of Contzen et al. (2025) explicit; a "not under consideration elsewhere" sentence was added. An AUTHOR CHECK comment asks you to confirm the word limit and section structure before claiming compliance.
13. **README.** HPC section now describes the pure-Python runs (jobs 48932578, 48932579; 64 CPUs, 240 GiB, 48 h) instead of the March R-backend run; TODO markers removed; EDC bullet states the exact threshold and k check.
14. **SLURM and docstrings.** `hpc/run_fig9_pyonly.slurm` header now matches the runs (64 CPUs, 240G, 48 h, qos 48h). `scripts/reproduce_fig9.py` example uses 64 workers. The `weatherisk cmip6` docstring gives both runs' counts.
15. **References.** Saunders et al. (2021) corrected to Extremes 24(2), 215–240, with A.G. Stephenson. Opitz (2013) added for the extremal-t model. Bauer et al. (2021), Nat. Clim. Change 11, 80–83, added for Destination Earth.
16. **Test tolerance tightening.** In `tests/test_r_cross_validation.py`, 18 `assert_allclose` calls had only `atol`, so numpy's default `rtol=1e-7` also applied. `rtol=0` was added so the stated absolute tolerances are the ones enforced; all 29 tests pass.

Open items for you (not changed):

- **Word count and structure.** The paper is about 6,000 words. The Elsevier SoftwareX page mentions a 3,000-word limit; the Guide for Authors was not read. The top-level "Model validation" section is not one of the five required sections. Decide whether to shorten and move validation under Software description or Impact.
- **CRediT and AI-use statements.** SoftwareX normally asks for a CRediT author statement and a declaration on generative AI use. Neither is in the paper.
- **CMIP6 data citation.** The AWI-ESM-1-1-LR historical data (ESGF) should be cited with its dataset DOI and the CMIP6 acknowledgement text.
- **Example 1 (stripes, Fig. 2).** RELEASE BLOCKER: rerun job 48947462 on albedo; the k values, figure and description in the paper date from a run before R was removed.
- **Version and DOI.** RELEASE BLOCKER: v0.2.0 tag and Zenodo DOI.


## H. Example 1 rerun and the status of the R code (2026-11-11)

1. **Example 1 is now the pure-Python run.** Albedo job 48947462 ran `scripts/reproduce_fig3.py --seed 42 --workers 32` on the `paper_stripes` preset. Wall time 8:24:36; largest process 11.6 GiB. It gives k_edc = 7 and k_lec = 5, the same counts as the old figure. The EDC partition is identical to the old one; the LEC partition differs in detail but has the same form.
2. **Measured on the run's saved arrays** (`labels_*.npy`, `inclusters_*.npy`, `b_matrix.npy`). Assigning each cell its cluster's estimate of b gives a mean absolute error of 0.373 (LEC) and 0.858 (EDC) over the 2601 cells. LEC in-cluster b, left to right: 0.37, 1.68, 3.42, 4.60, plus a 28-cell cluster in the upper-right corner with b = 5.07. EDC in-cluster b ranges from 0.31 to 3.73. 24 of 51 grid columns lie in a single LEC cluster; none lie in a single EDC cluster.
3. **Figure.** `figures/fig3_stripes.png` is replotted from those arrays with `plot_cluster_comparison`, without the "Fig. 3 —" title. The old `fig3_regenerated.png` is kept in the repository but no longer used by the paper.
4. **Text.** The paragraph now gives these numbers instead of qualitative claims. It also states that Contzen et al. obtain six clusters with both methods on their data, and that the counts need not agree because the replicates come from a different random number generator. Their Figure 3 discussion (EDC clusters span a wide x-range, LEC clusters resemble vertical stripes, small LEC clusters at the edge of the map) matches what we see, and the paper says so.
5. **The R code is not published by Contzen.** Contzen et al. (2025, Extremes 28:713–737) have only a data availability statement (for the Section 5 data, Danek et al. 2020) and no code availability statement. The arXiv preprint 2306.07780 has no code link. Zenodo has one Contzen software record, the R code for a different paper ("Long-term temporal evolution of extreme temperature in a warming Earth", 10.5281/zenodo.7413932, 2022). The Acknowledgements claim that he made the code "publicly available" is therefore removed. The sentence now reads "We thank J. Contzen for the original R implementation of the methods, on which weatherisk is based." The unverified "helpful discussions" was dropped; add it back if it is true.
6. **But `r_code/` is public through our repository.** The public `main` branch of github.com/wkchikach99-byte/weatherisk has contained `r_code/` (34 files) since the initial commit of 2026-03-09. The paper links to `r_code/functions.R#L515`. The repository's MIT licence names only Kchikach, Dickhaus and Lohmann as copyright holders, and `r_code/` has no licence of its own. `r_code/main.pdf` is the author manuscript of Contzen et al. Decision for you: get J. Contzen's written agreement to redistribute `r_code/` and say under which licence (for example a `r_code/LICENSE` or README line naming him as author), or remove it from the public tree and cite the R behaviour without a link. The published article is open access under CC BY 4.0, so redistributing the article is allowed with attribution, but `main.pdf` is a pre-publication manuscript and is not needed in a software repository; consider removing it. Nothing in the repository was changed.
7. **SoftwareX rules, from the official template (Version 6, March 2026).** The Guide for Authors page on ScienceDirect returns 403 to automated access, so it was not read. The DOCX and TEX templates at legacyfileshare.elsevier.com/promis_misc/softwarex-osp-template.{docx,tex} were read. They say:
   - "Our word limit is 4,000 words", and at most six pages "excluding metadata, tables, figures, references", with priority on the word limit. Exceptions "should be discussed with the Managing Editor or publisher prior to submission."
   - The five sections (Motivation and significance; Software description; Illustrative examples; Impact; Conclusions) are mandatory, and a paper without them "will be returned or rejected". Extra sections are neither allowed nor forbidden; "Any submission that deviates widely will likely face desk rejection." The words testing and validation do not appear in the template.
   - Abstract "ca. 100 words"; at most six keywords; at most six figures; Code metadata table C1 to C8 mandatory; a GitHub repository with README.md and "a Licence.txt file" is mandatory.
8. **Measured against those rules** (word count of the LaTeX source with figures, tables, listings, math and commands removed):

   | Part | Words |
   |---|---|
   | Abstract | 217 |
   | 1 Motivation and significance | 876 |
   | 2 Software description | 1963 |
   | 3 Illustrative examples | 1041 |
   | 4 Model validation (extra section) | 959 |
   | 5 Impact | 491 |
   | 6 Conclusions | 589 |
   | Main body total | 5919 |

   The main body is about 1900 words over the limit and the abstract about twice the guideline. Figures (3) and keywords (5) are within the limits. Because "Model validation" is a numbered section, Impact and Conclusions are numbered 5 and 6 instead of 4 and 5. The repository has `LICENSE`, not `Licence.txt`.
9. **Statements added.** A CRediT statement above the Acknowledgements replaces the "Author contributions" paragraph and keeps its content: Kchikach (Methodology, Software, Validation, Investigation, Visualization, Writing – original draft); Dickhaus and Lohmann (Conceptualization, Supervision, Writing – review & editing). A Data availability section cites the AWI-ESM-1-1-LR historical dataset, Danek et al. (2020), doi 10.22033/ESGF/CMIP6.9328 (checked: the DOI resolves to the DKRZ CERA record). A generative-AI declaration uses Elsevier's required wording and names Claude (Anthropic) for code and text editing; an AUTHOR CHECK asks you to confirm the tools and purposes. Funding and Competing Interests were already present.

## I. Length rule: Guide for Authors versus template (2026-11-12)

1. **The current Guide for Authors could not be read.** sciencedirect.com/journal/softwarex/publish/guide-for-authors returns HTTP 403 to both the automated fetcher and the in-app browser (an Elsevier bot check). The Wayback Machine holds no successful capture of that URL. The latest successful capture of the old URL (elsevier.com/journals/softwarex/2352-7110/guide-for-authors) is from 26 October 2022.
2. **The 2022 guide says:** "Maximum word count: 3000 excluding: title, authors, affiliations, references, metadata tables and including: abstract, running text, captions, footnotes. Maximum six (6) figures." It also says "We will only consider software articles that are submitted using SoftwareX's template" and asks authors to "strictly adhere to the Template for Original Software Publications."
3. **The current template (copyright 2026) says:** "Our word limit is 4,000 words." and "Maximum 6 pages in total (excluding metadata, tables, figures, references) with a 4000-word limit (we ask that more priority is placed on the word limit versus the page count)." It adds "in exceptional circumstances, we can be flexible with the page numbers and word limit." The abstract is "Ca. 100 words".
4. **Reading.** The template is newer than the 2022 guide, and the guide defers to the template, so 4,000 words is most likely the current limit. A web search snippet of an Elsevier journal page still says 3,000; that page now returns 404, so its date is unknown. Whether the 4,000 words include the abstract and captions is not stated in the template; the 2022 guide counted both.
5. **v23 measured both ways** (LaTeX source; tables, listings and display math removed; inline math and citations count as one word each; script in the session scratchpad, wc/count.py):

   | Part | Running text | Captions |
   |---|---|---|
   | Abstract | 215 | |
   | 1 Motivation and significance | 874 | |
   | 2 Software description | 2092 | 190 |
   | 3 Illustrative examples | 1059 | 246 |
   | 4 Model validation | 970 | 172 |
   | 5 Impact | 494 | |
   | 6 Conclusions | 588 | |
   | Total | 6292 | 608 |

   Sections 1 to 6 without abstract and captions are 6077 words. The §H count (5919) used a stricter filter for inline code; the difference does not change the conclusion. Against the 2026 template (4,000) the paper is about 2,000 words over. Against the 2022 rule (3,000 including abstract and captions, total 6,900) it is about 3,900 over.
6. **Safe target.** A version of at most 4,000 words including abstract and captions satisfies the template under either reading. Confirming the exact rule needs the current guide, opened by a person in a logged-in browser, or a question to the editorial office.

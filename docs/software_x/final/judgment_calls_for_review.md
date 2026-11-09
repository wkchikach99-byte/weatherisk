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
   With 5 starts Python was worse than R by more than 0.03 at some cells. The sweep found that 50 plain LHS starts is the smallest tested budget that meets the criterion on 10 selection and 10 held-out seeds (worst margin +0.0287, close to 0.03).
   Plain LHS designs are not nested, so 40 starts failed where 20 passed. A nested design was being tested at the time of writing.
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
7. **Test counts.** The paper reports 238 tests (Level 1: 116, Level 2: 16, Level 3: 106). With R removed from the search path, the 237 tests that run by default pass in 67 s on an Apple M4 laptop; the remaining test is the slow end-to-end run. Recount if tests change before the tag.
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

   Without detrending both counts are at or next to the published ones. With STL, k_LEC is 14, far from 24. k_LEC depends on the optimizer through the local estimates, so 23 versus 24 is within what we would expect from a different optimizer (inference, not tested). The March 2026 STL run with the old R backend gave k_LEC 21 and k_EDC 142. Both jobs completed on 64 CPUs (measured, `sacct`): STL 36:02:48 wall time, no-detrend 25:58:03. Memory needs two numbers. `sacct` MaxRSS for the whole job is 95 GiB (STL) and 85 GiB (no-detrend). The largest single process, from `/usr/bin/time -v`, is 38.2 GB (STL) and 39.8 GB (no-detrend). An earlier version of this note called the per-process figure the job's max RSS, which was wrong. The paper's "approximately 4 hours" was removed on 2026-11-09 and Table 2 now gives the measured times. Step 6 runs one worker per cluster, LEC clusters first, so the largest LEC cluster (3247 and 6230 cells) ran alone for hours while the other cores idled (measured with `ps` inside the jobs on 2026-11-08; average CPU use 2214 % and 3016 % of 6400 %). The cluster counts in the final `pipeline_results.npz` equal the step-5 counts above.

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

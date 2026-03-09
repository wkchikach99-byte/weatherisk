# Optimization Learnings

## 2026-03-09 — Fix 1: Clustering memory retention

- Change: CMIP6 clustering now reuses the condensed distance vector for thresholding and linkage, lowers the default LEC chunk size, and drops `dm_*` / `hc_*` artifacts by default unless explicitly requested.
- Why: This reduces peak memory without changing the clustering math, threshold rule, linkage method, or labels.
- Validation: full test suite passed (`148 passed`).
- Benchmark case: synthetic full-pipeline benchmark, `36 years`, `12x12`, `4 workers`.
- Result: total time `5.400s`, peak process-tree RSS `0.583 GiB`.
- Interpretation: the memory fix did not materially slow the benchmarked pipeline, but it also did not make clustering faster. This is expected because the refactor targets retention and duplication, not arithmetic cost.
- Gotcha: benchmarking multiprocessing code from a heredoc or stdin entrypoint fails on macOS/Python 3.14 spawn mode. A real script entrypoint is required.
- Gotcha: monkeypatched tests that return synthetic clustering artifacts must mirror the new default behavior, otherwise they can report false regressions.

## 2026-03-09 — Fix 2: Vectorized pair assembly

- Change: replaced Python list-building in CMIP6 local MLE setup and shared density helpers with equivalent NumPy reshape, tile, and repeat operations.
- Why: this removes Python overhead while keeping the pair ordering, likelihood values, optimizer, and parameter bounds unchanged.
- Validation: full test suite passed (`148 passed`).
- Benchmark case: synthetic full-pipeline benchmark, `36 years`, `12x12`, `4 workers`.
- Result: total time `5.309s`, peak process-tree RSS `0.583 GiB`.
- Interpretation: the benchmark signal was modest. Total runtime improved slightly, but the local-MLE step itself stayed within noise on this small synthetic case. That suggests optimizer cost still dominates over pair-array setup at this scale.
- Gotcha: small multiprocessing benchmarks are noisy enough that a single run should be interpreted qualitatively, not as proof of a large speedup.

## 2026-03-09 — Fix 3: Condensed EDC path

- Change: when CMIP6 clustering artifacts are not retained, the EDC step now computes the condensed distance vector directly instead of building the full square matrix first.
- Why: linkage and thresholding only need the upper triangle, so the square matrix was an avoidable memory allocation in the common pipeline path.
- Validation: full test suite passed (`148 passed`) and an explicit test now checks condensed EDC output against `squareform(_edc_matrix_flat(...))`.
- Benchmark case: synthetic full-pipeline benchmark, `36 years`, `12x12`, `4 workers`.
- Result: total time `5.316s`, peak process-tree RSS `0.579 GiB`.
- Interpretation: the synthetic benchmark shows a small additional memory reduction and a slightly faster clustering step. The total runtime remained effectively flat because clustering is still a small share of wall-clock time on this benchmark size.
- Gotcha: a clear memory improvement may show up more strongly on full-scale grids than on small synthetic benchmarks, where process overhead and timing jitter dominate.

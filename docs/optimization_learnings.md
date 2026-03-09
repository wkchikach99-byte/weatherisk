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

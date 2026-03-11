# Figure 9 Parity Progress Log

Date: 2026-03-11

## Scope Note

This document records operational CMIP6 Figure 9 validation work only.
It is not the authority for strict constituent-function Python-to-R parity.

The authoritative strict acceptance standard is:

- `docs/python_r_parity_migration_plan.md`

Any statement in this file about "parity" means mini-pipeline or operational
Figure 9 agreement under the dedicated CMIP6 validation suite, not the higher
standard required by the migration plan.

## Goal

Bring the Python and Rust CMIP6 Figure 9 pipeline into high-accuracy operational parity with the historical R workflow used in the Extremes paper, then validate a clean full-data rerun.

## Current status

- The CMIP6 parity test module now passes fully.
- The local clean Figure 9 rerun was launched in an isolated output directory and kept separate from older outputs.
- A clean albedo launcher script has been prepared for the production full-data rerun.

## Verified Operational Validation State

Command used for validation:

```bash
/Users/wikchi001/Desktop/weatherisk/.venv/bin/pytest tests/test_cmip6_rparity.py -q
```

Result:

```text
30 passed in 26.63s
```

This means the dedicated CMIP6 operational validation suite is green for both
backends.

It does not by itself imply that all constituent functions satisfy the strict
raw-R parity contract in `docs/python_r_parity_migration_plan.md`.

## Root-cause fixes completed

### 1. Local MLE neighborhood traversal

- The CMIP6 local-estimation path was not traversing neighbors in the same order as the original R stencil logic.
- This mattered because the local likelihood surface is often flat or weakly identified, so accumulation order can change the selected optimum.
- The Python path now follows the original R row/column grid-index traversal exactly.

### 2. Fréchet precision normalization at the local-MLE boundary

- The upstream detrending, annual maxima, and Fréchet steps were already extremely close to R, but last-bit differences were still enough to push some local optimizations onto different equivalent minima.
- Normalizing Fréchet inputs to the effective historical R precision before local MLE closed the remaining mini-pipeline parity gap.

### 3. Test expectation alignment

- One older parity assertion still expected Python to beat R's local MLE at many cells.
- That is not the right requirement for strict constituent-function parity work.
- The test was updated so the requirement is now: Python/Rust must not be materially worse than the historical R result, while the full chain must reproduce the stored reference outputs.

## What now matches

- STL detrending on the CMIP6 mini reference.
- Annual maxima on the CMIP6 mini reference.
- GEV to Fréchet transformation closely enough for downstream parity.
- Raw local MLE estimates against the mini R reference.
- Smoothed local estimates against the mini R reference.
- Full mini EDC chain against the R reference.
- Full mini LEC chain against the R reference.
- Python and Rust backend agreement on the CMIP6 parity path.

## Full-run execution state

### Local clean run

- Output directory: `output/cmip6_fig9_clean_20260311_105434`
- Log file: `output/cmip6_fig9_clean_20260311_105434/run.log`

Reason for keeping it:

- It preserves a clean local artifact trail.
- It is useful for immediate debugging and inspection.

Reason it is not the preferred production path:

- The full Figure 9 rerun is a much better fit for albedo.
- The local process showed very low current CPU activity during Step 1a, so it is not the efficient route for the final production rerun.

### Recommended production run target

- Use albedo for the clean full-data Figure 9 rerun after strict parity claims
	are interpreted via `docs/python_r_parity_migration_plan.md` rather than this
	operational note.
- Keep the local run only as a secondary diagnostic artifact.

## Clean albedo launcher

Prepared script:

- `hpc/run_fig9_clean.slurm`

Recommended submit command:

```bash
sbatch --export=ALL,FIG9_RUN_ID=cmip6_fig9_clean_20260311,FIG9_OUTPUT_DIR=output/cmip6_fig9_clean_20260311 hpc/run_fig9_clean.slurm
```

Expected log files on albedo:

- `logs/fig9_clean_<jobid>.out`
- `logs/fig9_clean_<jobid>.err`

## Interpretation relative to the Extremes paper

The current state supports the following claim precisely:

- The operational CMIP6 Figure 9 path now reproduces the historical R
	reference chain at the dedicated CMIP6 validation-test level, including the
	previously failing local-estimation and LEC stages.

This document does not claim that strict constituent-function parity is
complete. That claim must be evaluated only through
`docs/python_r_parity_migration_plan.md`.

The remaining task is not algorithmic parity repair.
The remaining task is a clean full-data reproduction run and output inspection against the paper-level expectations.

## Next actions

1. Submit the clean albedo Figure 9 job.
2. Monitor the albedo log until the run completes.
3. Inspect produced `k_LEC`, `k_EDC`, saved arrays, and final figures against the paper and previous R outputs.
4. Record any remaining paper-level differences, if any, as output-validation findings rather than parity-code defects.
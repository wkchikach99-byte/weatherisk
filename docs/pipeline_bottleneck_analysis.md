# Pipeline Data Flow & Bottleneck Analysis

> **Context**: Reproducing Figure 9 from Contzen et al. (2025, *Extremes* 28:713–737)
> using AWI-ESM-1-1-LR historical precipitation on Albedo HPC (16 CPUs).
>
> **Current runtime**: ~4–6 hours.  **Target**: <30 minutes.

---

## 1. End-to-End Data Flow

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                    RAW INPUT                                    │
 │  165 NetCDF files × 12 months each                             │
 │  Shape after loading: (1872 months, 96 lat, 192 lon)           │
 │  = 1872 × 18,432 cells = 34.5 M float64 values (~276 MB)      │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 1a: DE-TREND                                    ~2s  ✅  │
 │                                                                 │
 │  Input:  (1872, 96, 192)  float64 3-D array                    │
 │  Method: Vectorised NumPy — monthly climatology + running mean  │
 │  Output: (1872, 96, 192)  de-trended precipitation              │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: implicit (NumPy BLAS)                             │
 │  Bottleneck? NO — already fast                                  │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 1b: ANNUAL MAXIMA                               <1s  ✅  │
 │                                                                 │
 │  Input:  (1872, 96, 192)                                        │
 │  Method: Group by year → max per group (12 months → 1 value)   │
 │  Output: (156 years, 96, 192)                                   │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: implicit (NumPy)                                  │
 │  Bottleneck? NO                                                 │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 2: GEV FIT → FRÉCHET TRANSFORM              ~885s  🔴   │
 │                                                                 │
 │  Input:  (156, 96, 192)  reshaped to (156, 18432)              │
 │  Method: For EACH of 18,432 cells:                              │
 │          • scipy.stats.genextreme.fit(column)  → (loc, sc, sh) │
 │          • scipy.stats.genextreme.cdf(column)  → probabilities │
 │          • Z = -1/log(P)                        → Fréchet       │
 │  Output: (156, 18432) Fréchet-transformed values                │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: ❌ SERIAL LOOP — one cell at a time              │
 │                                                                 │
 │  WHY IT'S SLOW:                                                 │
 │  • genextreme.fit() does numerical MLE internally (Nelder-Mead) │
 │  • Called 18,432 times × ~50 iterations each                    │
 │  • Each call has Python overhead: input validation, wrapping    │
 │  • ZERO parallelism — 15 of 16 CPUs sit idle                   │
 │                                                                 │
 │  Bottleneck? YES — 15 min on 1 core, could be ~1 min on 16     │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 3: LOCAL MLE (PAIRWISE COMPOSITE LIKELIHOOD)  ~3-5h  🔴🔴│
 │                                                                 │
 │  THE DOMINANT BOTTLENECK — ~80% of total pipeline time          │
 │                                                                 │
 │  Input:  (156, 18432) Fréchet data + (18432, 2) grid coords    │
 │                                                                 │
 │  For EACH of 18,432 cells:                                      │
 │    1. Find neighbours within ε=5 grid-point radius              │
 │       → typically ~78 neighbours per cell                       │
 │    2. Build pair arrays: zi, zj, x_lag, y_lag                   │
 │       → 78 neighbours × 156 years = ~12,168 pairs              │
 │    3. Run 3 multi-start L-BFGS-B optimisations:                │
 │       Each iteration evaluates:                                 │
 │       ┌──────────────────────────────────────────────────────┐  │
 │       │  pairwise_density_summand(zi, zj, x, y, df, α,a,b,γ)│  │
 │       │                                                      │  │
 │       │  For EACH of ~12,168 pairs:                          │  │
 │       │    • cov_fkt_2d(x, y, α, a, b, γ)  — 1 exp() call   │  │
 │       │    • (z2/z1)^(1/df), (z1/z2)^(1/df) — 2 pow() calls │  │
 │       │    • scipy.stats.t.pdf(m, df+1)     — 2 calls  🐌    │  │
 │       │    • scipy.stats.t.cdf(m, df+1)     — 2 calls  🐌    │  │
 │       │    • _dtdiff(m, df+1)               — 2 calls        │  │
 │       │    • ~30 arithmetic operations                       │  │
 │       │                                                      │  │
 │       │  Total per iteration: ~12K × 40 ops = ~500K FLOPs    │  │
 │       └──────────────────────────────────────────────────────┘  │
 │       × ~50-200 L-BFGS-B iterations per start                  │
 │       × 3 starts                                                │
 │       = ~75M–300M FLOPs PER CELL                                │
 │                                                                 │
 │  Total: 18,432 × ~200M = ~3.7 TRILLION floating-point ops      │
 │                                                                 │
 │  Output: (18432, 3) — estimated (a, b, γ) per cell             │
 │                                                                 │
 │  Data type: dense NumPy ndarrays (all in-memory)                │
 │  Parallelism: ✅ multiprocessing.Pool(16) with imap_unordered  │
 │                                                                 │
 │  WHY IT'S STILL SLOW DESPITE 16 WORKERS:                       │
 │                                                                 │
 │  🐌 scipy.stats.t.pdf() and t.cdf() are the #1 cost centre:   │
 │     • They call _validate_args → _argcheck → broadcasting      │
 │     • Then call betainc (for CDF) — a compiled but complex fn  │
 │     • For df=6 (fixed!), we could replace with a closed-form   │
 │       formula that's 5–10× faster                               │
 │                                                                 │
 │  🐌 Python multiprocessing IPC overhead:                        │
 │     • Each result is pickled → sent via pipe → unpickled        │
 │     • With chunksize=20, that's 18432/20 = 921 IPC round-trips │
 │                                                                 │
 │  🐌 maxiter=10000 is overkill:                                  │
 │     • L-BFGS-B typically converges in 50–200 iterations         │
 │     • The extra headroom wastes time on pathological cells      │
 │                                                                 │
 │  Bottleneck? YES — THE bottleneck. Everything else is noise.    │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 4: SPATIAL SMOOTHING                          ~10s  ✅   │
 │                                                                 │
 │  Input:  (18432, 3) local estimates                             │
 │  Method: Moving average within radius=2 grid points             │
 │          Angular wrapping for γ                                 │
 │  Output: (18432, 3) smoothed estimates                          │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: ❌ serial loop, but N is small enough             │
 │  Bottleneck? NO — fast enough as-is                             │
 └─────────────────────────┬───────────────────────────────────────┘
                           │
                           ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 5a: LEC DISSIMILARITY MATRIX                 ~5-10min 🟡 │
 │                                                                 │
 │  Input:  (18432, 3) smoothed estimates                          │
 │  Method: Rasterise normalised ellipses on half-circle grid      │
 │          Jaccard (1 - IoU) for each pair → 18432² matrix        │
 │  Output: (18432, 18432) symmetric dissimilarity matrix          │
 │          = 2.7 GB float64                                       │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: ✅ chunked vectorisation (256 rows at a time)    │
 │                                                                 │
 │  WHY IT'S MODERATELY SLOW:                                      │
 │  • Intermediate arrays mask_i, mask_j: (256, 18432, ~300)       │
 │    = ~1.3 GB per chunk — memory-bound, causes cache thrashing   │
 │  • 72 chunks × broadcast comparisons                            │
 │                                                                 │
 │  Bottleneck? MODERATE — could be improved with numba            │
 └────────────────────────┬────────────────────────────────────────┘
                          │
                          ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 5b: EDC DISSIMILARITY MATRIX                 ~5-10min 🟡 │
 │                                                                 │
 │  Input:  (156, 18432) Fréchet data                              │
 │  Method: Rank each cell's 156 values                            │
 │          Madogram: mean |rank_i - rank_j| for all pairs         │
 │          Convert to extremal coefficient                        │
 │  Output: (18432, 18432) symmetric dissimilarity matrix          │
 │                                                                 │
 │  Data type: dense NumPy ndarray                                 │
 │  Parallelism: ❌ SERIAL Python loop over 18,432 rows           │
 │                                                                 │
 │  WHY IT'S UNNECESSARILY SLOW:                                   │
 │  • Uses a manual loop: for i in range(n_cells-1): ...           │
 │  • scipy.spatial.distance.cdist('cityblock') does the same      │
 │    thing but in compiled C — 10–50× faster                      │
 │  • This exact optimisation already exists in clustering.py's    │
 │    c_extrcoeff_matrix() but IS NOT USED in the pipeline!        │
 │                                                                 │
 │  Bottleneck? YES — easy fix (use cdist)                         │
 └────────────────────────┬────────────────────────────────────────┘
                          │
                          ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 5c: HIERARCHICAL CLUSTERING                    ~30s  ✅  │
 │                                                                 │
 │  Input:  (18432, 18432) dissimilarity matrices (LEC + EDC)      │
 │  Method: scipy.cluster.hierarchy.linkage('average')             │
 │          + quantile threshold → k clusters                      │
 │  Output: labels_lec (k≈24), labels_edc (k≈104)                │
 │                                                                 │
 │  Data type: condensed distance vector (169M entries)            │
 │  Parallelism: ✅ scipy linkage is compiled C                   │
 │  Bottleneck? NO                                                 │
 └────────────────────────┬────────────────────────────────────────┘
                          │
                          ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 6: IN-CLUSTER RE-ESTIMATION                  ~30-60min 🔴│
 │                                                                 │
 │  Input:  Fréchet data + cluster labels                          │
 │  Method: For EACH cluster (k_LEC=24, k_EDC=104):               │
 │          • Extract cells belonging to cluster                   │
 │          • Build ALL pairwise (z_i, z_j) within cluster         │
 │          • Run global MLE: same pairwise_density_summand()      │
 │            with 3 multi-start L-BFGS-B                          │
 │  Output: (a, b, γ) per cluster                                 │
 │                                                                 │
 │  Data type: dense NumPy ndarrays                                │
 │  Parallelism: ❌ SERIAL loop over clusters                     │
 │                                                                 │
 │  WHY IT'S SLOW:                                                 │
 │  • 104 EDC clusters + 24 LEC clusters = 128 MLE optimisations  │
 │  • Each cluster may have ~177 cells (18432/104)                 │
 │  • Pairwise within cluster: C(177,2) × 156 = ~2.4M pairs       │
 │  • Same slow scipy.stats.t.pdf/cdf bottleneck as Step 3         │
 │  • All 128 clusters run sequentially — 15 CPUs idle             │
 │                                                                 │
 │  Bottleneck? YES — parallelise across clusters                  │
 └────────────────────────┬────────────────────────────────────────┘
                          │
                          ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │  STEP 7: PLOT GENERATION                             ~10s  ✅  │
 │                                                                 │
 │  Input:  labels + lat/lon grid                                  │
 │  Method: matplotlib + cartopy world maps                        │
 │  Output: PNG/PDF figures                                        │
 │                                                                 │
 │  Bottleneck? NO                                                 │
 └─────────────────────────────────────────────────────────────────┘
```

---

## 2. Bottleneck Summary

| Step | Current Time | % of Total | Root Cause | Fix Difficulty |
|------|-------------|------------|------------|----------------|
| **3. Local MLE** | **3–5 hours** | **~80%** | `scipy.stats.t.pdf/cdf` called billions of times through L-BFGS-B | Medium |
| **2. GEV fit** | **15 min** | **~5%** | Serial loop, no parallelism | Easy |
| **6. In-cluster** | **30–60 min** | **~10%** | Serial over clusters, same slow density | Easy |
| **5b. EDC matrix** | **5–10 min** | **~3%** | Python loop instead of `cdist` | Trivial |
| **5a. LEC matrix** | **5–10 min** | **~2%** | Memory-bound chunked broadcast | Medium |

**Total current**: ~4–6 hours.

---

## 3. The Inner Hot Loop — Why `pairwise_density_summand` Is So Slow

The absolute #1 cost is inside `pairwise_density_summand()`, called from Steps 3 and 6.
Here is a breakdown of what happens per function call (on ~12,000 pairs):

```
pairwise_density_summand(zi, zj, x, y, df=5, α=1, a, b, γ)
│
├── cov_fkt_2d(x, y, 1, a, b, γ)        →  exp(-sqrt(Q))        ~5 μs
│     Pure NumPy arithmetic — FAST
│
├── (z2/z1)**(1/5), (z1/z2)**(1/5)       →  np.power()           ~10 μs
│     Two element-wise power calls — FAST
│
├── scipy.stats.t.pdf(m1, 6)             →  SLOW                 ~200 μs  🐌
│     Internally: _validate_args → _logpdf → exp → ...
│     For df=6, this is just:
│       C₆ × (1 + x²/6)^(-7/2)
│     which is ONE line of NumPy!
│
├── scipy.stats.t.cdf(m1, 6)             →  SLOW                 ~500 μs  🐌🐌
│     Internally: _validate_args → betainc (special function)
│     For df=6 (integer), this can be computed with a polynomial!
│
├── scipy.stats.t.pdf(m2, 6)             →  another              ~200 μs  🐌
├── scipy.stats.t.cdf(m2, 6)             →  another              ~500 μs  🐌🐌
│
├── _dtdiff(m1, 6)                       →  gamma_fn import!     ~100 μs  🐌
│     from scipy.special import gamma as gamma_fn  ← EVERY CALL
│     The import is cached but still has overhead
│
├── _dtdiff(m2, 6)                       →  another              ~100 μs  🐌
│
├── ~30 arithmetic ops (+, ×, /)          →  np vectorised         ~20 μs
│
└── np.log + np.maximum                   →  np vectorised          ~5 μs

Total per call: ~1,640 μs ≈ 1.6 ms
```

**Key insight**: `t.pdf` + `t.cdf` account for **~85% of the time** in this function.
For fixed `df=6`, both have closed-form expressions that can be computed in ~50 μs
with pure NumPy. That's a **10–15× speedup** on the inner loop alone.

---

## 4. Would Polars Help?

### What Polars Is Good At

Polars excels at:
- **DataFrame operations**: filter, group-by, join, aggregate on tabular data
- **Lazy evaluation**: query planning, predicate pushdown, projection pruning
- **Multi-threaded execution**: parallel column operations on Arrow-backed data
- **String/categorical processing**: much faster than pandas

### Where Our Pipeline's Data Lives

| Step | Data Structure | Shape | Type |
|------|---------------|-------|------|
| Load | xarray Dataset → NumPy 3-D array | (1872, 96, 192) | Dense float64 |
| De-trend | NumPy 3-D array | (1872, 96, 192) | Dense float64 |
| Annual max | NumPy 3-D → 2-D array | (156, 18432) | Dense float64 |
| GEV fit | NumPy column → scipy.stats | (156,) per cell | Dense float64 |
| Local MLE | NumPy vectors → scipy.optimize | (12168,) per cell | Dense float64 |
| Dissimilarity | NumPy 2-D matrix | (18432, 18432) | Dense float64 |
| Clustering | scipy condensed vector | (169M,) | Dense float64 |

### Verdict: Polars Would NOT Help Here

**The data is never tabular.** Our pipeline operates on:
1. **3-D grids** (time × lat × lon) — Polars doesn't handle N-D arrays
2. **Dense numerical matrices** (dissimilarity, covariance) — this is NumPy/BLAS territory
3. **Per-cell numerical optimisation** (L-BFGS-B) — this is scipy.optimize territory
4. **Scientific special functions** (t-distribution PDF/CDF) — this is scipy.stats territory

Polars would only be useful if we were doing things like:
- Joining metadata tables
- Filtering/grouping CSV logs
- String processing on station names

**None of our bottlenecks involve tabular data processing.** Converting to Polars
would add complexity without addressing the actual performance problems.

### What WOULD Help (In Order of Impact)

| Technology | Where | Expected Speedup | Why |
|-----------|-------|-------------------|-----|
| **Numba `@jit`** | `pairwise_density_summand` | **10–20×** | JIT-compiles the inner loop to native code. No scipy overhead. Eliminates `t.pdf/cdf` call overhead entirely |
| **Direct t-dist formulas** | `pairwise_density_summand` | **5–10×** | Replace `scipy.stats.t.pdf/cdf(x, 6)` with 2 lines of NumPy. Zero Python dispatch overhead |
| **multiprocessing for GEV** | `_compute_frechet_global` | **10–15×** | Embarrassingly parallel — each cell is independent |
| **multiprocessing for clusters** | `_incluster_reestimate_cmip6` | **10–15×** | Each cluster's MLE is independent |
| **`scipy.spatial.distance.cdist`** | `_edc_matrix_flat` | **10–50×** | Replace Python loop with compiled C. Already exists in `clustering.py`! |
| **Reduce `maxiter`** | `_local_mle_one_cmip6` | **2–3×** | 10000 → 2000. L-BFGS-B converges in <200 iterations |
| **SLURM array job** | Step 3 | **N×** | Split 18,432 cells across N separate SLURM jobs. Each job handles 18432/N cells. Embarrassingly parallel across nodes |

---

## 5. The Ideal Fast Path

```
Current:  ~4–6 hours on 16 CPUs
                    │
                    ▼
Fix 1: Direct t-distribution formulas (no scipy.stats overhead)
        Step 3: 3–5h → 30–60 min       (5–10× on inner loop)
        Step 6: 30–60 min → 5–10 min
                    │
                    ▼
Fix 2: Parallelize GEV fitting (Pool across 16 CPUs)
        Step 2: 15 min → 1–2 min        (16× parallelism)
                    │
                    ▼
Fix 3: Parallelize in-cluster re-estimation
        Step 6: 5–10 min → 1 min         (16× parallelism)
                    │
                    ▼
Fix 4: Use cdist for EDC matrix
        Step 5b: 5–10 min → 10–30 sec   (compiled C)
                    │
                    ▼
Fix 5: Reduce maxiter 10000 → 2000
        Step 3: 30–60 min → 15–30 min   (2× less wasted iterations)
                    │
                    ▼
Fix 6: Numba @jit on pairwise_density_summand
        Step 3: 15–30 min → 5–10 min    (JIT native code)
                    │
                    ▼
                ~15–25 minutes total on 16 CPUs

Fix 7 (optional): SLURM array job across multiple nodes
        Step 3: 5–10 min → 1–2 min      (N nodes)
                    │
                    ▼
                ~5–10 minutes total
```

---

## 6. Recommended Implementation Order

### Phase 1 — Quick Wins (1–2 hours of coding)

1. **Replace `scipy.stats.t.pdf/cdf` with direct formulas** in `density.py`
2. **Reduce `maxiter` from 10000 to 2000** in `cmip6_pipeline.py`
3. **Parallelize GEV fitting** in `_compute_frechet_global()`
4. **Use `cdist`** in `_edc_matrix_flat()` (copy pattern from `clustering.py`)
5. **Parallelize in-cluster re-estimation** in `_incluster_reestimate_cmip6()`

**Expected result**: 4–6 hours → **20–40 minutes**

### Phase 2 — Numba (2–4 hours of coding)

6. **`@numba.jit(nopython=True)`** on `pairwise_density_summand` and `cov_fkt_2d`
7. **`@numba.jit(parallel=True)`** with `prange` for the outer cell loop

**Expected result**: 20–40 min → **5–15 minutes**

### Phase 3 — Multi-Node (if needed)

8. **SLURM array job**: split cells into chunks, each chunk = one SLURM task
9. **Merge results**: collect partial `.npy` files and concatenate

**Expected result**: 5–15 min → **2–5 minutes** (scales with node count)

---

## 7. Why NOT These Alternatives?

| Alternative | Why Not |
|-------------|---------|
| **Polars** | Data is never tabular. All bottlenecks are NumPy arrays + scipy optimisation. Polars can't help. |
| **Dask** | Our arrays are small enough to fit in memory (276 MB). Dask's overhead outweighs any benefit for 18K cells. Could help for multi-node Step 3, but SLURM array jobs are simpler. |
| **GPU (CuPy/JAX)** | The inner loop is an optimiser calling a density function. GPUs excel at bulk matrix ops but not at sequential L-BFGS-B iterations with Python callbacks. Would need a full rewrite to JAX's `jax.scipy.optimize`. |
| **pandas** | Same as Polars — data isn't tabular. |
| **R (original)** | The R code has the same bottleneck (it uses `optim(method="L-BFGS-B")` with the same density). R is typically 2–5× slower than NumPy for the same operations. |
| **C/Fortran extension** | Maximum performance, but high development cost. Numba gives 80% of the benefit with 10% of the effort. |

# Methodology Notes — CPC Pipeline Figures

These notes explain the step-by-step pipeline that produces the figures
in `docs/figures/`, assess their mathematical and climatological
correctness, and discuss interpretation caveats.  Intended as a
reference for the SoftwareX paper text.

---

## Pipeline Walkthrough: From Raw Data to Figures

### Step 1: Load CPC Data

- **Input**: NOAA CPC daily precipitation NetCDF files (2000–2019), 0.5° resolution
- **Sub-region**: lat [30°N, 65°N], lon [5°E, 55°E] — Europe + Middle East
- **Coarsening**: every 4th cell → ~2° grid
- Region matches the map extents in the figures. Coarsening factor is
  standard for computational tractability.

### Step 2: Block Maxima → GEV → Fréchet

- **Block maxima**: 365-day blocks → 20 annual maxima per cell
- **GEV fit**: Per-cell MLE via `scipy.stats.genextreme.fit()`
- **Fréchet transform**: Z = −1 / log F(x; μ, σ, ξ)  where F is the
  fitted GEV CDF
- **Clipping**: Fréchet values clipped to [0.05, n_blocks²]
  (Padoan et al. 2010 safeguard)
- **Land mask**: cells with >50% NaN are excluded
- **Coordinate normalisation**: geographic coords mapped to [−5, 5] × [−5, 5]
  — matches the R code convention exactly (verified by unit tests)

**Assessment**: Mathematically correct. The GEV→Fréchet marginal
transform is the standard approach (de Haan & Ferreira 2006,
Padoan et al. 2010). The clipping at n² prevents extreme outliers
from dominating the likelihood — a known best practice.

### Step 3: Local MLE of (a, b, γ)

- For each grid cell, neighbours within normalised radius ε=3.0 are selected
- Pairwise composite log-likelihood is maximised over (a, b, γ) using
  multi-start L-BFGS-B (3 starts, Latin hypercube sampling)
- The density is the Schlather (2002) bivariate max-stable density with
  non-stationary covariance via Σ⁻¹ parametrised by the ellipse (a, b, γ)

**Assessment**: Correct implementation of the Justus (2025) §3.2
methodology. Cross-validated against R output to 1e-10 tolerance
(28 tests passed).

### Step 4: Spatial Smoothing

- Moving-average over neighbours within radius 2.0 (normalised units)
- Angular parameter γ uses wrapped averaging in [−π/2, π/2]

**Assessment**: Correct. Matches R code. Verified by
`test_smoothing_a_and_b` and `test_smoothing_gamma` tests.

### Step 5: LEC & EDC Clustering

**LEC** (→ `lec_clusters_map.png`):
- Jaccard-like ellipse-overlap dissimilarity D₂: rasterises each cell's
  normalised ellipse on a half-circle grid, computes 1 − IoU between
  all pairs
- Hierarchical agglomerative clustering (Ward's linkage)
- k determined by 30%-quantile threshold of pairwise distances → k=26

**EDC** (→ `edc_clusters_map.png`):
- Rank-based madogram → extremal coefficients → dissimilarity D₁
- Same clustering + threshold method → k=18

Implementation note: the current Python code uses SciPy average linkage
via `weatherisk.clustering.clustering()`, not Ward linkage.

**Assessment**: Implementation verified against R reference data
(cluster counts match, merge heights match, assignments equivalent).
The 30%-quantile threshold is exactly the Justus (2025) method.

### Step 6: In-cluster Re-estimation

- For each cluster, global pairwise CL MLE is run on all cell pairs
  within the cluster to get refined (a, b, γ)

### Step 7: Risk Metrics

**Figures produced**: `lec_risk_es_map.png`, `edc_risk_es_map.png`

The **loss variable** for cluster A at time t is:

    L_t^A = max_{s in A} Z_t(s)

where Z_t(s) is the Fréchet-transformed block maximum at cell s, time t.

- VaR₉₅ = F_L⁻¹(0.95) — the 95th percentile of the spatial block maxima
- ES₉₅ = E[L | L > VaR₉₅] — the conditional mean above VaR

**Assessment**:
- Mathematically correct. ES is the standard coherent risk measure
  (Artzner et al. 1999, McNeil et al. 2005).
- The ES is on the **Fréchet scale** — see interpretation section below.
- Most clusters show ES₉₅ in the 300–400 range, with a few at 50–250.
  With 20 annual maxima and spatial max over many cells, the 95th
  percentile picks the 2nd-largest value; ES averages the top 1–2.
  Larger clusters have higher spatial maxima.

### Step 8: GDP-Weighted Risk

**Figures produced**: `lec_risk_gdp_map.png`, `edc_risk_gdp_map.png`,
`gdp_exposure_map.png`

    Risk(s) = ES₉₅^{k(s)} × GDP(s)

- Each cell s inherits the ES of its cluster k(s)
- GDP PPP per cell from Kummu et al. (2018), regridded by summing all
  fine (5 arc-min) cells within each coarsened pipeline cell boundary
- The product gives a USD-denominated risk index per cell

### Summary Panel

**Figure produced**: `lec_edc_summary_panel.png`

6-panel figure: LEC clusters, EDC clusters, parameter a, parameter b,
LEC ES₉₅, EDC ES₉₅.

## Current Code Mapping

For the implemented CPC maps pipeline, the corresponding internal
functions are:

- Step 1: `weatherisk.cpc_pipeline._load_subregion`
- Step 2: `weatherisk.cpc_pipeline._compute_frechet`
- Step 3: `weatherisk.cpc_pipeline._run_local_estimation`
- Step 4: `weatherisk.cpc_pipeline._smooth_estimates`
- Step 5: `weatherisk.cpc_pipeline._run_clustering`
- Step 6: `weatherisk.cpc_pipeline._incluster_reestimate`
- Step 7: `weatherisk.cpc_pipeline._cluster_risk`
- Step 8: optional GDP weighting inside `weatherisk.cpc_pipeline.run_cpc_pipeline`
- Plot generation: `weatherisk.cpc_pipeline.generate_maps`

---

## Key Interpretation Points

### Why ES is on the Fréchet Scale, Not in Physical Precipitation (mm/day)

In Step 2, we transform each cell's annual block maxima from physical
precipitation values (mm/day) into **unit Fréchet margins** via:

    Z = −1 / log F_GEV(x; μ, σ, ξ)

This is a **probability integral transform** — it maps each cell's data
onto a common heavy-tailed scale where the marginal distribution is
identical everywhere: P(Z ≤ z) = exp(−1/z).  After this transform, a
value of Z=20 means the same tail probability at every cell, regardless
of whether that cell receives 50 mm/day or 200 mm/day in physical terms.

**Why we do this:** The whole point of max-stable process modelling is to
separate **marginal behaviour** (how extreme is a single cell?) from
**dependence structure** (how do extremes co-occur across space?).  The
Fréchet transform removes the marginal differences so the pairwise
likelihood estimation (Step 3) captures only the spatial dependence —
the ellipse parameters (a, b, γ) describe *how* extremes spread across
space, not *how big* they are at any one location.

**What this means for ES:** When we compute ES₉₅ on the Fréchet data,
we're measuring the **severity of joint spatial extremes relative to
their own tail behaviour**, not in mm/day.  An ES₉₅ of 350 means: "when
the worst cell in this cluster exceeds its own 95th-percentile Fréchet
level, the average exceedance is 350 Fréchet units."  This captures
**dependence intensity** — clusters where many cells tend to be extreme
simultaneously will have higher spatial maxima and hence higher ES.

**Why not back-transform to mm/day?** Because each cell has a different
GEV.  A cluster's "loss" is L_t = max_{s in A} Z_t(s) — the max is
taken over cells with different marginals.  To express this in mm/day
you'd need to pick *which cell's* marginal to invert through, and the
answer would depend on which cell happened to be the worst in each year.
The Fréchet scale gives a **common currency** for comparing tail risk
across clusters that cover different climatic regimes.


### Why GDP × ES is a Hazard-Exposure Index, Not a Direct Loss Estimate

The formula is:

    Risk(s) = ES₉₅^{k(s)} × GDP(s)

For this to be a **direct monetary loss**, you'd need ES to represent
something like "fraction of GDP destroyed by an extreme event."  But ES₉₅
is a **dimensionless Fréchet-scale number** (typically 50–400 in our
maps).  It doesn't represent a damage fraction, a probability, or a
dollar amount.  It represents the **intensity of spatial tail dependence**
in the cluster.

So the product ES × GDP gives units of "Fréchet-units × USD" — a
**composite index** that ranks cells by "how intense are the joint
extremes in this region × how much economic value is at stake."  It
answers: **where do strong spatial dependence and high economic exposure
coincide?**

This is analogous to how seismic hazard maps multiply a ground-shaking
intensity index by population exposure — the result isn't "expected
deaths" but a **risk prioritisation score**.

**Why we plot it this way:** The value of this map is in **relative
comparison** across cells, not in the absolute numbers.  It correctly
identifies Western Europe (Germany, France, Benelux) as the highest-risk
zone because those cells have both:
- Moderate-to-high ES (precipitation extremes are spatially coherent —
  when it's extreme, it's extreme over large areas)
- Very high GDP (lots of economic assets are exposed)

Meanwhile, Central Asian cells with *higher* ES but negligible GDP
correctly show as low-risk.  And Mediterranean cells with high GDP but
*lower* ES (extremes are more localised) show as moderate-risk.


### What Would Make It a Direct Loss Estimate?

To get actual expected monetary losses, you'd need a **damage function**
D(x) that maps physical precipitation intensity x (mm/day) to a damage
fraction (0–1), then compute:

    Loss(s) = E[ D(X(s)) × GDP(s) ]

That requires: (1) back-transforming from Fréchet to physical units,
(2) a vulnerability/damage curve specific to precipitation (e.g., flood
damage functions from insurance or catastrophe models).  This is outside
the scope of the current methodology, which focuses on the **statistical
dependence structure** rather than impact modelling.

**Recommendation for the paper:** State clearly that the GDP-risk maps
provide a *spatially explicit hazard-exposure index* that identifies
regions where strong extremal dependence coincides with economic
concentration, rather than claiming they estimate actual losses.  This is
honest and still very useful for risk prioritisation.

---

## Figure-to-File Mapping

| Figure description                        | File in `docs/figures/`          |
|-------------------------------------------|----------------------------------|
| LEC spatial dependence clusters (k=26)    | `lec_clusters_map.png`           |
| EDC spatial dependence clusters (k=18)    | `edc_clusters_map.png`           |
| Estimated parameter a (dependence range)  | `param_a_map.png`                |
| Estimated parameter b (anisotropy)        | `param_b_map.png`                |
| Estimated parameter γ (rotation)          | `param_gamma_map.png`            |
| Tail-risk intensity ES₉₅, LEC clusters   | `lec_risk_es_map.png`            |
| Tail-risk intensity ES₉₅, EDC clusters   | `edc_risk_es_map.png`            |
| GDP-weighted risk, LEC clusters           | `lec_risk_gdp_map.png`           |
| GDP-weighted risk, EDC clusters           | `edc_risk_gdp_map.png`           |
| GDP exposure layer (Kummu et al. 2018)    | `gdp_exposure_map.png`           |
| 6-panel summary (clusters + params + ES)  | `lec_edc_summary_panel.png`      |
| Fig. 3 reproduction (stripes, smoothed)   | `fig3_smoothed.png`              |
| Fig. 3 full resolution                    | `fig3_full.png`                  |
| Fig. 3 sanity check (low-res)             | `fig3_sanity_check.png`          |

---

## Potential Concerns and Caveats

1. **ES₉₅ with only 20 observations**: With 20 annual block maxima per
   cluster, VaR₉₅ = the 19th order statistic, and ES₉₅ averages the
   top 1–2 values.  This is a small-sample extreme, making ES estimates
   **noisy** (high variance).  The paper should acknowledge this.

2. **Cluster-size effect on ES**: Larger clusters (more cells) produce
   higher spatial maxima simply because max over more cells is
   stochastically larger.  The ES₉₅ map is therefore partly a proxy for
   cluster size.  The GDP weighting partially mitigates this.

3. **Risk units are an index**: See interpretation section above.

4. **No spatial contiguity in clustering**: The LEC/EDC clustering uses
  average linkage on the dissimilarity matrix without spatial contiguity
   constraints.  Some clusters may contain geographically disconnected
   cells (visible in `lec_clusters_map.png`).  This is methodologically
   consistent with Justus (2025) but worth discussing.

5. **GEV fit quality**: No diagnostic plots or goodness-of-fit tests for
   the per-cell GEV fits.  With only 20 annual maxima, the GEV shape
   parameter estimates will be noisy.  A Q-Q plot or Anderson-Darling
   test per cell (or at least a sample) would strengthen validation.
   → See GEV diagnostics script (to be generated).

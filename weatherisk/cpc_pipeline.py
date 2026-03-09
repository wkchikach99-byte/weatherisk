"""Real-data pipeline: CPC climate data → LEC / EDC clusters → risk metrics.

Implements the full methodology from:
    Justus (2025), *Extremal dependence and local estimation clustering
    for non-stationary max-stable processes*, Extremes.

Pipeline steps:
    1. Load NOAA CPC daily NetCDF files for a year range
    2. Select geographic sub-region, coarsen grid
    3. Annual block maxima (365 d) → parametric GEV → unit Fréchet
    4. Local pairwise composite-likelihood MLE → (a, b, γ)
    5. Spatial moving-average smoothing  (angular wrap for γ)
    6. LEC / EDC clustering  (30 %-quantile threshold)
    7. In-cluster re-estimation of (a, b, γ)
    8. Risk metrics: VaR / ES per cluster  (added value)

All heavy computation is in :func:`run_cpc_pipeline`; plotting is
delegated to :mod:`weatherisk.map_plotting`.

Usage from CLI::

    weatherisk maps --data-dir data/netcdf --output-dir docs/figures

Usage from Python::

    from weatherisk.cpc_pipeline import PipelineConfig, run_cpc_pipeline
    cfg = PipelineConfig(data_dir="data/netcdf")
    result = run_cpc_pipeline(cfg)
"""

from __future__ import annotations

import os
import time
import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import minimize
from scipy.stats import qmc, rankdata

from weatherisk.extremes import block_maxima, fit_gev, to_frechet
from weatherisk.density import pairwise_density_summand
from weatherisk.clustering import (
    calc_distance_ellipses,
    clustering,
    cluster_number_threshold_method,
)
from weatherisk.risk import compute_var, compute_es


# ══════════════════════════════════════════════════════════════════
#  Configuration
# ══════════════════════════════════════════════════════════════════

@dataclass
class PipelineConfig:
    """All tuneable parameters for the CPC real-data pipeline.

    Defaults mirror the Justus (2025) Extremes paper, §4.2.
    """

    # I/O
    data_dir: str = "data/netcdf"
    output_dir: str = "docs/figures"
    year_start: int = 2000
    year_end: int = 2020          # exclusive

    # Climate variable
    variable: str = "precip"      # NetCDF variable name
    file_prefix: str = "precip"   # file naming: {prefix}.{year}.nc

    # Sub-region  (avoids prime-meridian wrap)
    lat_range: tuple[float, float] = (30.0, 65.0)
    lon_range: tuple[float, float] = (5.0, 55.0)

    # Resolution
    coarsen: int = 4              # every Nth cell → ~2° from 0.5°
    block_size: int = 365         # annual block maxima

    # Max-stable hyper-parameters  (§4.2)
    df: float = 5.0               # ν  degrees of freedom
    alpha: float = 1.0            # α  smoothness exponent

    # Local estimation  (§3.2) — normalised-coordinate units
    neighbor_radius: float = 3.0  # ε  neighbourhood for pairwise CL
    smoothing_radius: float = 2.0 # spatial smoothing
    mle_ensemble: int = 3         # multi-start L-BFGS-B

    # Clustering  (§3.3)
    quantile_threshold: float = 0.30  # 30 %-quantile of pairwise dists

    # Risk
    risk_level: float = 0.95      # p  for VaR / ES

    # GDP exposure  (optional)
    gdp_path: str | None = None   # path to Kummu et al. GDP NetCDF
    gdp_year: int = 2015          # year snapshot for exposure

    # Plotting
    dpi: int = 300

    @property
    def years(self) -> range:
        return range(self.year_start, self.year_end)

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """Map extent (lon_min, lon_max, lat_min, lat_max) with padding."""
        return (
            self.lon_range[0] - 2,
            self.lon_range[1] + 2,
            self.lat_range[0] - 2,
            self.lat_range[1] + 2,
        )

    @property
    def variable_label(self) -> str:
        """Human-readable variable name for titles."""
        _labels = {"tmax": "Max Temperature", "precip": "Precipitation",
                   "tmin": "Min Temperature"}
        return _labels.get(self.variable, self.variable)

    @property
    def unit(self) -> str:
        """Physical unit string."""
        _units = {"tmax": "°C", "precip": "mm/day", "tmin": "°C"}
        return _units.get(self.variable, "")


# ══════════════════════════════════════════════════════════════════
#  Step 1.  Load CPC data — sub-region & coarsen
# ══════════════════════════════════════════════════════════════════

def _load_subregion(cfg: PipelineConfig, *, verbose: bool = True):
    """Return *(daily, lats_1d, lons_1d)* for the coarsened sub-region."""
    import xarray as xr

    if verbose:
        print("=" * 60)
        print(f"  Step 1 : Load NOAA CPC {cfg.variable}")
        print("=" * 60)

    all_daily: list[np.ndarray] = []
    sub_lats = sub_lons = None

    for year in cfg.years:
        path = os.path.join(cfg.data_dir, f"{cfg.file_prefix}.{year}.nc")
        if not os.path.exists(path):
            if verbose:
                print(f"  WARNING: {path} not found — skipping")
            continue

        if verbose:
            print(f"  {year} … ", end="", flush=True)
        ds = xr.open_dataset(path)

        lats = ds.lat.values
        if lats[0] > lats[-1]:
            lat_sel = slice(float(cfg.lat_range[1]), float(cfg.lat_range[0]))
        else:
            lat_sel = slice(float(cfg.lat_range[0]), float(cfg.lat_range[1]))
        lon_sel = slice(float(cfg.lon_range[0]), float(cfg.lon_range[1]))

        da = ds[cfg.variable].sel(lat=lat_sel, lon=lon_sel)
        da = da.isel(lat=slice(None, None, cfg.coarsen),
                     lon=slice(None, None, cfg.coarsen))

        data = da.values.astype(np.float64)
        if sub_lats is None:
            sub_lats = da.lat.values
            sub_lons = da.lon.values

        all_daily.append(data)
        ds.close()
        if verbose:
            print(f"shape {data.shape}, "
                  f"valid {np.isfinite(data).mean():.0%}")

    daily = np.concatenate(all_daily, axis=0)
    if verbose:
        print(f"\n  Combined: {daily.shape}  "
              f"lat [{sub_lats.min():.1f}, {sub_lats.max():.1f}]  "
              f"lon [{sub_lons.min():.1f}, {sub_lons.max():.1f}]")
    return daily, sub_lats, sub_lons


# ══════════════════════════════════════════════════════════════════
#  Step 2.  Block maxima → GEV → Fréchet
# ══════════════════════════════════════════════════════════════════

def _compute_frechet(daily, lats, lons, cfg: PipelineConfig, *,
                     verbose: bool = True):
    """Return (frechet, bm, land_idx, coords, geo_coords, gev_par)."""
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 2 : Block maxima → GEV → Fréchet")
        print("=" * 60)

    n_days, n_lat, n_lon = daily.shape
    n_cells = n_lat * n_lon
    daily_flat = daily.reshape(n_days, n_cells)

    nan_frac = np.isnan(daily_flat).mean(axis=0)
    land_mask = nan_frac < 0.50
    land_idx = np.where(land_mask)[0]
    if verbose:
        print(f"  Grid cells : {n_cells}   land : {len(land_idx)}"
              f"  ({len(land_idx) / n_cells:.0%})")

    daily_land = daily_flat[:, land_idx].copy()
    for c in range(daily_land.shape[1]):
        col = daily_land[:, c]
        m = np.nanmean(col)
        col[np.isnan(col)] = m

    bm = block_maxima(daily_land, block_size=cfg.block_size)
    n_blocks = bm.shape[0]
    if verbose:
        print(f"  Block size : {cfg.block_size} d → {n_blocks} annual maxima")
        print(f"  Block-max range : [{bm.min():.1f}, {bm.max():.1f}] {cfg.unit}")

    n_land = len(land_idx)
    frechet = np.empty((n_blocks, n_land))
    gev_par = np.empty((n_land, 3))
    ok = np.ones(n_land, dtype=bool)

    n_fail = 0
    for c in range(n_land):
        try:
            loc, scale, shape = fit_gev(bm[:, c])
            gev_par[c] = [loc, scale, shape]
            fr = to_frechet(bm[:, c], loc, scale, shape)
            fr = np.clip(fr, 0.01, None)
            frechet[:, c] = fr
        except Exception:
            ok[c] = False
            n_fail += 1

    ok &= np.all(np.isfinite(frechet), axis=0)
    frechet = frechet[:, ok]
    land_idx = land_idx[ok]
    gev_par = gev_par[ok]
    bm = bm[:, ok]

    # Clip extreme Fréchet tails  (Padoan et al. 2010)
    frechet = np.clip(frechet, 0.05, float(n_blocks ** 2))

    # Normalise coordinates to [-5, 5] × [-5, 5]
    lat_grid = np.repeat(lats, len(lons))
    lon_grid = np.tile(lons, len(lats))
    raw_lat = lat_grid[land_idx]
    raw_lon = lon_grid[land_idx]

    lat_min, lat_max = raw_lat.min(), raw_lat.max()
    lon_min, lon_max = raw_lon.min(), raw_lon.max()
    norm_lat = -5.0 + 10.0 * (raw_lat - lat_min) / max(lat_max - lat_min, 1e-6)
    norm_lon = -5.0 + 10.0 * (raw_lon - lon_min) / max(lon_max - lon_min, 1e-6)

    coords = np.column_stack([norm_lat, norm_lon])
    geo_coords = np.column_stack([raw_lat, raw_lon])

    if verbose:
        print(f"  GEV failures : {n_fail}")
        print(f"  Valid cells : {len(land_idx)}")
        print(f"  Fréchet shape: {frechet.shape}  "
              f"range [{frechet.min():.2f}, {frechet.max():.2f}]")
        print(f"  Normalised coords: "
              f"lat [{coords[:, 0].min():.1f}, {coords[:, 0].max():.1f}]"
              f"  lon [{coords[:, 1].min():.1f}, {coords[:, 1].max():.1f}]")
    return frechet, bm, land_idx, coords, geo_coords, gev_par


# ══════════════════════════════════════════════════════════════════
#  Step 3.  Local pairwise CL MLE  (a, b, γ)
# ══════════════════════════════════════════════════════════════════

def _local_mle_one(frechet, cidx, coords, cfg: PipelineConfig):
    """Estimate (a, b, γ) at a single cell."""
    n_blocks = frechet.shape[0]

    dlat = coords[:, 0] - coords[cidx, 0]
    dlon = coords[:, 1] - coords[cidx, 1]
    dists = np.sqrt(dlat ** 2 + dlon ** 2)
    nb = np.where((dists > 0.01) & (dists <= cfg.neighbor_radius))[0]

    if len(nb) < 3:
        return np.array([1.0, 0.0, 0.0])

    z_c = frechet[:, cidx]
    zi_l, zj_l, xl_l, yl_l = [], [], [], []
    for j in nb:
        z_nb = frechet[:, j]
        zi_l.extend(z_nb)
        zj_l.extend(z_c)
        xl_l.extend([coords[j, 1] - coords[cidx, 1]] * n_blocks)
        yl_l.extend([coords[j, 0] - coords[cidx, 0]] * n_blocks)

    zi = np.asarray(zi_l)
    zj = np.asarray(zj_l)
    xl = np.asarray(xl_l)
    yl = np.asarray(yl_l)

    good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
    zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]
    if len(zi) < 5:
        return np.array([1.0, 0.0, 0.0])

    lo = np.array([0.01, 0.0, -np.pi / 2])
    hi = np.array([15.0, 15.0, np.pi / 2])

    from weatherisk.backend import neg_log_likelihood_sum as _nll_sum

    def neg_llh(p):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            v = _nll_sum(zi, zj, xl, yl, cfg.df, cfg.alpha, p[0], p[1], p[2])
        return v if np.isfinite(v) else 1e20

    sampler = qmc.LatinHypercube(d=3, seed=42 + cidx)
    starts = qmc.scale(sampler.random(n=max(cfg.mle_ensemble, 5)), lo, hi)

    best_v, best_p = np.inf, np.array([1.0, 0.0, 0.0])
    for s in range(cfg.mle_ensemble):
        try:
            r = minimize(
                neg_llh, starts[s], method="L-BFGS-B",
                bounds=list(zip(lo, hi)),
                options={"maxiter": 10000, "ftol": 1e-10},
            )
            if r.fun < best_v:
                best_v, best_p = r.fun, r.x.copy()
        except Exception:
            pass
    return best_p


def _run_local_estimation(frechet, coords, cfg: PipelineConfig, *,
                          verbose: bool = True):
    """Step 3: local MLE at every cell."""
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 3 : Local MLE  (a, b, γ)")
        print("=" * 60)

    n = frechet.shape[1]
    est = np.zeros((n, 3))
    t0 = time.time()
    for c in range(n):
        if verbose and c % max(1, n // 20) == 0:
            print(f"    {c + 1:4d}/{n}  ({time.time() - t0:.0f} s)")
        est[c] = _local_mle_one(frechet, c, coords, cfg)

    if verbose:
        print(f"  a ∈ [{est[:, 0].min():.3f}, {est[:, 0].max():.3f}]")
        print(f"  b ∈ [{est[:, 1].min():.3f}, {est[:, 1].max():.3f}]")
        print(f"  γ ∈ [{np.degrees(est[:, 2]).min():.1f}°,"
              f" {np.degrees(est[:, 2]).max():.1f}°]")
    return est


# ══════════════════════════════════════════════════════════════════
#  Step 4.  Spatial smoothing
# ══════════════════════════════════════════════════════════════════

def _smooth_estimates(est, coords, cfg: PipelineConfig, *,
                      verbose: bool = True):
    """Spatial moving-average smoothing; γ wrapped in [-π/2, π/2]."""
    if verbose:
        print("\n" + "=" * 60)
        print(f"  Step 4 : Spatial smoothing  (radius {cfg.smoothing_radius})")
        print("=" * 60)

    n = len(est)
    out = np.empty_like(est)
    for c in range(n):
        d = np.sqrt((coords[:, 0] - coords[c, 0]) ** 2 +
                    (coords[:, 1] - coords[c, 1]) ** 2)
        nb = np.where(d <= cfg.smoothing_radius)[0]

        nbe = est[nb]
        out[c, 0] = nbe[:, 0].mean()
        out[c, 1] = nbe[:, 1].mean()

        cg = est[c, 2]
        ang = nbe[:, 2].copy()
        ang = np.where(ang < cg - np.pi / 2, ang + np.pi, ang)
        ang = np.where(ang > cg + np.pi / 2, ang - np.pi, ang)
        mg = ang.mean()
        if mg < -np.pi / 2:
            mg += np.pi
        elif mg > np.pi / 2:
            mg -= np.pi
        out[c, 2] = mg

    if verbose:
        print(f"  Smoothed a ∈ [{out[:, 0].min():.3f}, {out[:, 0].max():.3f}]")
        print(f"  Smoothed b ∈ [{out[:, 1].min():.3f}, {out[:, 1].max():.3f}]")
    return out


# ══════════════════════════════════════════════════════════════════
#  Step 5.  LEC / EDC clustering
# ══════════════════════════════════════════════════════════════════

def _edc_matrix(frechet):
    """Rank-based madogram extremal-coefficient matrix D₁."""
    n_blocks, n_cells = frechet.shape
    ranks = np.column_stack(
        [rankdata(frechet[:, s]) for s in range(n_cells)]
    ).T  # (n_cells, n_blocks)

    ec = np.zeros((n_cells, n_cells))
    for i in range(n_cells - 1):
        diff = np.abs(ranks[i] - ranks[i + 1:])
        v = diff.mean(axis=1) / (2.0 * (n_blocks + 1))
        denom = 1.0 - 2.0 * v
        denom[denom <= 0] = 1e-12
        ec[i, i + 1:] = np.minimum(1.0, (1.0 + 2.0 * v) / denom - 1.0)
    return ec + ec.T


def _run_clustering(smoothed, frechet, cfg: PipelineConfig, *,
                    verbose: bool = True):
    """Step 5: LEC (D₂) and EDC (D₁) with quantile threshold."""
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 5 : LEC & EDC clustering")
        print("=" * 60)

    # --- LEC ---
    if verbose:
        print("  Computing LEC dissimilarity …")
    dm_lec = calc_distance_ellipses(smoothed, res=21)
    hc_lec = clustering(dm_lec)
    vec_lec = dm_lec[np.triu_indices_from(dm_lec, k=1)]
    thr_lec = np.quantile(vec_lec, cfg.quantile_threshold)
    k_lec = cluster_number_threshold_method(hc_lec, thr_lec)
    k_lec = max(2, k_lec)
    labels_lec = fcluster(hc_lec, t=k_lec, criterion="maxclust")
    if verbose:
        print(f"  LEC → k = {k_lec}  "
              f"(30 %-quantile threshold = {thr_lec:.3f})")

    # --- EDC ---
    if verbose:
        print("  Computing EDC dissimilarity …")
    dm_edc = _edc_matrix(frechet)
    hc_edc = clustering(dm_edc)
    vec_edc = dm_edc[np.triu_indices_from(dm_edc, k=1)]
    thr_edc = np.quantile(vec_edc, cfg.quantile_threshold)
    k_edc = cluster_number_threshold_method(hc_edc, thr_edc)
    k_edc = max(2, k_edc)
    labels_edc = fcluster(hc_edc, t=k_edc, criterion="maxclust")
    if verbose:
        print(f"  EDC → k = {k_edc}  "
              f"(30 %-quantile threshold = {thr_edc:.5f})")

    return dict(
        labels_lec=labels_lec, k_lec=k_lec,
        hc_lec=hc_lec, dm_lec=dm_lec,
        labels_edc=labels_edc, k_edc=k_edc,
        hc_edc=hc_edc, dm_edc=dm_edc,
    )


# ══════════════════════════════════════════════════════════════════
#  Step 6.  In-cluster re-estimation
# ══════════════════════════════════════════════════════════════════

def _incluster_reestimate(frechet, coords, labels, cfg, tag, *,
                          verbose: bool = True):
    """Re-estimate (a, b, γ) per cluster via global pairwise CL MLE."""
    from weatherisk.density import pairwise_density_optim

    unique_cl = sorted(np.unique(labels))
    if verbose:
        print(f"  {tag}: re-estimating {len(unique_cl)} clusters …")

    results: dict[int, np.ndarray] = {}
    for cl in unique_cl:
        mask = labels == cl
        n_cl = int(mask.sum())
        if n_cl < 3:
            results[cl] = np.array([1.0, 0.0, 0.0])
            continue

        z_cl = frechet[:, mask].T
        X_cl = coords[mask, 1]
        Y_cl = coords[mask, 0]

        try:
            est = pairwise_density_optim(
                z_cl, cfg.df, cfg.alpha, X_cl, Y_cl,
                upper_bounds=(15.0, 15.0),
                max_dist=4.0 * cfg.neighbor_radius,
                ensemble=3,
            )
        except Exception:
            est = np.array([1.0, 0.0, 0.0])

        results[cl] = est
        if verbose:
            print(f"    cl {cl:2d} ({n_cl:3d} cells)  "
                  f"a={est[0]:.3f}  b={est[1]:.3f}  "
                  f"γ={np.degrees(est[2]):.1f}°")
    return results


# ══════════════════════════════════════════════════════════════════
#  Step 7.  Risk metrics
# ══════════════════════════════════════════════════════════════════

def _cluster_risk(frechet, labels, p: float = 0.95):
    """VaR and ES on the Fréchet scale, per cluster.

    Loss: L_t = max_{s in A} Z_t(s)  (spatial block max).
    """
    out = []
    for cl in sorted(np.unique(labels)):
        mask = labels == cl
        cmax = frechet[:, mask].max(axis=1)
        out.append(dict(
            cluster=int(cl),
            n_cells=int(mask.sum()),
            var=compute_var(cmax, p),
            es=compute_es(cmax, p),
        ))
    return out


# ══════════════════════════════════════════════════════════════════
#  Public orchestrator
# ══════════════════════════════════════════════════════════════════

def run_cpc_pipeline(
    cfg: PipelineConfig | None = None,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Execute the full CPC real-data pipeline.

    Parameters
    ----------
    cfg : PipelineConfig, optional
        Pipeline configuration.  Uses defaults if *None*.
    verbose : bool
        Print progress banners.

    Returns
    -------
    dict
        Keys: ``lat``, ``lon``, ``frechet``, ``bm``, ``smoothed``,
        ``labels_lec``, ``labels_edc``, ``k_lec``, ``k_edc``,
        ``risk_lec``, ``risk_edc``, ``incl_lec``, ``incl_edc``,
        ``hc_lec``, ``hc_edc``, ``dm_lec``, ``dm_edc``,
        ``gev_params``, ``config``.
    """
    if cfg is None:
        cfg = PipelineConfig()

    t_start = time.time()
    if verbose:
        print("\n" + "=" * 60)
        print(f"  LEC / EDC Risk Pipeline  —  real CPC {cfg.variable} data")
        print("  Following Justus (2025, Extremes) methodology")
        print("=" * 60)

    # 1
    daily, lats, lons = _load_subregion(cfg, verbose=verbose)

    # 2
    frechet, bm, land_idx, coords, geo_coords, gev_p = \
        _compute_frechet(daily, lats, lons, cfg, verbose=verbose)
    lat_v = geo_coords[:, 0]
    lon_v = geo_coords[:, 1]
    # keep reference to full daily array for grid shape in return dict
    _daily_shape = daily.shape

    # 3
    est = _run_local_estimation(frechet, coords, cfg, verbose=verbose)

    # 4
    sm = _smooth_estimates(est, coords, cfg, verbose=verbose)

    # 5
    cl = _run_clustering(sm, frechet, cfg, verbose=verbose)

    # 6
    if verbose:
        print("\n" + "=" * 60)
        print("  Step 6 : In-cluster re-estimation")
        print("=" * 60)
    incl_lec = _incluster_reestimate(
        frechet, coords, cl["labels_lec"], cfg, "LEC", verbose=verbose)
    incl_edc = _incluster_reestimate(
        frechet, coords, cl["labels_edc"], cfg, "EDC", verbose=verbose)

    # 7
    if verbose:
        print("\n" + "=" * 60)
        print(f"  Step 7 : Risk metrics  "
              f"(VaR_{cfg.risk_level:.0%} / ES_{cfg.risk_level:.0%})")
        print("=" * 60)
    risk_lec = _cluster_risk(frechet, cl["labels_lec"], cfg.risk_level)
    risk_edc = _cluster_risk(frechet, cl["labels_edc"], cfg.risk_level)

    if verbose:
        print("\n  ── LEC ──")
        for r in risk_lec:
            print(f"    cl {r['cluster']:2d}  ({r['n_cells']:3d} cells)  "
                  f"VaR={r['var']:.2f}  ES={r['es']:.2f}")
        print("\n  ── EDC ──")
        for r in risk_edc:
            print(f"    cl {r['cluster']:2d}  ({r['n_cells']:3d} cells)  "
                  f"VaR={r['var']:.2f}  ES={r['es']:.2f}")

    # 8  GDP-weighted cell-level risk  (optional)
    gdp_per_cell = None
    risk_gdp_lec = None
    risk_gdp_edc = None

    if cfg.gdp_path is not None:
        if verbose:
            print("\n" + "=" * 60)
            print("  Step 8 : GDP-weighted cell-level risk")
            print("=" * 60)

        from weatherisk.gdp import gdp_for_land_cells as _gdp_for_land

        gdp_per_cell = _gdp_for_land(
            cfg.gdp_path, lats, lons, land_idx,
            year=cfg.gdp_year, verbose=verbose,
        )

        # Build ES lookup per cell: each cell inherits its cluster's ES
        es_map_lec = np.zeros(len(land_idx))
        for r in risk_lec:
            es_map_lec[cl["labels_lec"] == r["cluster"]] = r["es"]

        es_map_edc = np.zeros(len(land_idx))
        for r in risk_edc:
            es_map_edc[cl["labels_edc"] == r["cluster"]] = r["es"]

        # Risk(s) = ES_95^{k(s)} × GDP(s)
        risk_gdp_lec = es_map_lec * gdp_per_cell
        risk_gdp_edc = es_map_edc * gdp_per_cell

        if verbose:
            def _fmt(v):
                if v >= 1e12:
                    return f"${v / 1e12:.1f}T"
                if v >= 1e9:
                    return f"${v / 1e9:.1f}B"
                if v >= 1e6:
                    return f"${v / 1e6:.1f}M"
                return f"${v:,.0f}"

            print(f"\n  GDP-weighted risk  (ES₉₅ × GDP per cell):")
            print(f"    LEC  total = {_fmt(risk_gdp_lec.sum())}")
            print(f"    LEC  max cell = {_fmt(risk_gdp_lec.max())}")
            print(f"    EDC  total = {_fmt(risk_gdp_edc.sum())}")
            print(f"    EDC  max cell = {_fmt(risk_gdp_edc.max())}")

    elapsed = time.time() - t_start
    if verbose:
        print(f"\n  Pipeline completed in {elapsed:.0f} s  "
              f"({elapsed / 60:.1f} min)")

    return dict(
        lat=lat_v,
        lon=lon_v,
        frechet=frechet,
        bm=bm,
        smoothed=sm,
        estimates=est,
        labels_lec=cl["labels_lec"],
        labels_edc=cl["labels_edc"],
        k_lec=cl["k_lec"],
        k_edc=cl["k_edc"],
        risk_lec=risk_lec,
        risk_edc=risk_edc,
        incl_lec=incl_lec,
        incl_edc=incl_edc,
        hc_lec=cl["hc_lec"],
        hc_edc=cl["hc_edc"],
        dm_lec=cl["dm_lec"],
        dm_edc=cl["dm_edc"],
        gev_params=gev_p,
        config=cfg,
        # grid metadata — needed for filled-region (pcolormesh) maps
        lats_1d=lats,
        lons_1d=lons,
        n_lat=_daily_shape[1],
        n_lon=_daily_shape[2],
        land_idx=land_idx,
        # GDP-weighted risk  (None when gdp_path not set)
        gdp_per_cell=gdp_per_cell,
        risk_gdp_lec=risk_gdp_lec,
        risk_gdp_edc=risk_gdp_edc,
    )


# ══════════════════════════════════════════════════════════════════
#  Plot generation helper  (called by CLI)
# ══════════════════════════════════════════════════════════════════

def generate_maps(result: dict[str, Any], *, verbose: bool = True) -> list[str]:
    """Generate all standard Cartopy maps from pipeline results.

    Parameters
    ----------
    result : dict
        Output of :func:`run_cpc_pipeline`.
    verbose : bool
        Print confirmation for each saved file.

    Returns
    -------
    list of str
        Paths to all saved PDF files.
    """
    from weatherisk.map_plotting import (
        plot_cluster_map_geo,
        plot_field_map,
        plot_risk_map,
        plot_summary_panel,
    )

    cfg: PipelineConfig = result["config"]
    out = cfg.output_dir
    os.makedirs(out, exist_ok=True)
    ext = cfg.extent
    saved: list[str] = []

    if verbose:
        print("\n" + "=" * 60)
        print("  Generating Cartopy maps")
        print("=" * 60)

    def _save(path):
        saved.append(path)
        if verbose:
            print(f"  ✓ {os.path.basename(path)}")

    lat = result["lat"]
    lon = result["lon"]
    sm = result["smoothed"]

    # Grid metadata for filled-region (pcolormesh) rendering
    gk = dict(
        lats_1d=result.get("lats_1d"),
        lons_1d=result.get("lons_1d"),
        n_lat=result.get("n_lat"),
        n_lon=result.get("n_lon"),
        land_idx=result.get("land_idx"),
    )

    vl = cfg.variable_label  # e.g. "Precipitation"
    yr = f"{cfg.year_start}–{cfg.year_end - 1}"

    # 1. LEC clusters
    p = os.path.join(out, "lec_clusters_map.png")
    plot_cluster_map_geo(
        lat, lon, result["labels_lec"], result["k_lec"],
        extent=ext, **gk,
        title=(f"Spatial Dependence Clusters — {vl} Extremes\n"
               f"LEC (Ellipse-Overlap D₂), k={result['k_lec']}  ({yr})"),
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 2. EDC clusters
    p = os.path.join(out, "edc_clusters_map.png")
    plot_cluster_map_geo(
        lat, lon, result["labels_edc"], result["k_edc"],
        extent=ext, **gk,
        title=(f"Spatial Dependence Clusters — {vl} Extremes\n"
               f"EDC (Madogram D₁), k={result['k_edc']}  ({yr})"),
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 3. Parameter a
    p = os.path.join(out, "param_a_map.png")
    plot_field_map(
        lat, lon, sm[:, 0], extent=ext, **gk,
        title=f"Estimated Dependence Range (Semi-Minor Axis a)\n{vl} Extremes  ({yr})",
        label="a",
        cmap="viridis", save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 4. Parameter b
    p = os.path.join(out, "param_b_map.png")
    plot_field_map(
        lat, lon, sm[:, 1], extent=ext, **gk,
        title=f"Estimated Anisotropy (Semi-Major Extension b)\n{vl} Extremes  ({yr})",
        label="b",
        cmap="magma", save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 5. Parameter γ
    p = os.path.join(out, "param_gamma_map.png")
    plot_field_map(
        lat, lon, np.degrees(sm[:, 2]), extent=ext, **gk,
        title=f"Estimated Dependence Orientation (Rotation γ)\n{vl} Extremes  ({yr})",
        label="γ  (degrees)",
        cmap="twilight", vmin=-90, vmax=90,
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 6. ES₉₅ per LEC cluster
    p = os.path.join(out, "lec_risk_es_map.png")
    plot_risk_map(
        lat, lon, result["labels_lec"], result["risk_lec"],
        metric="es", extent=ext, **gk,
        title=(f"Tail-Risk Intensity (ES₉₅) — {vl} Hazard, LEC Clusters\n"
               f"Expected Shortfall of Fréchet Block Maxima per Cluster  ({yr})"),
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 7. ES₉₅ per EDC cluster
    p = os.path.join(out, "edc_risk_es_map.png")
    plot_risk_map(
        lat, lon, result["labels_edc"], result["risk_edc"],
        metric="es", extent=ext, **gk,
        title=(f"Tail-Risk Intensity (ES₉₅) — {vl} Hazard, EDC Clusters\n"
               f"Expected Shortfall of Fréchet Block Maxima per Cluster  ({yr})"),
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    # 8. GDP-weighted risk maps  (if available)
    if result.get("risk_gdp_lec") is not None:
        gy = result['config'].gdp_year

        p = os.path.join(out, "lec_risk_gdp_map.png")
        plot_field_map(
            lat, lon, result["risk_gdp_lec"], extent=ext, **gk,
            title=(f"Exposure-Weighted {vl} Risk — LEC Clusters\n"
                   f"ES₉₅ × GDP PPP per Cell  ({yr}, GDP {gy})"),
            label="Risk  (USD PPP)",
            cmap="YlOrRd", save_path=p, dpi=cfg.dpi,
        )
        _save(p)

        p = os.path.join(out, "edc_risk_gdp_map.png")
        plot_field_map(
            lat, lon, result["risk_gdp_edc"], extent=ext, **gk,
            title=(f"Exposure-Weighted {vl} Risk — EDC Clusters\n"
                   f"ES₉₅ × GDP PPP per Cell  ({yr}, GDP {gy})"),
            label="Risk  (USD PPP)",
            cmap="YlOrRd", save_path=p, dpi=cfg.dpi,
        )
        _save(p)

        p = os.path.join(out, "gdp_exposure_map.png")
        plot_field_map(
            lat, lon, result["gdp_per_cell"], extent=ext, **gk,
            title=(f"Economic Exposure — GDP PPP per Grid Cell\n"
                   f"Kummu et al. (2018), {gy} Snapshot (Constant 2011 USD)"),
            label="GDP PPP  (USD)",
            cmap="YlGn", save_path=p, dpi=cfg.dpi,
        )
        _save(p)

    # 9. Summary panel
    p = os.path.join(out, "lec_edc_summary_panel.png")
    plot_summary_panel(
        lat, lon,
        result["labels_lec"], result["labels_edc"],
        result["k_lec"], result["k_edc"],
        sm, result["risk_lec"], result["risk_edc"],
        extent=ext, **gk,
        suptitle=(f"Spatial Extreme-Dependence Clustering & Tail Risk\n"
                  f"CPC {vl} ({yr}) — LEC (k={result['k_lec']}) / "
                  f"EDC (k={result['k_edc']})"),
        save_path=p, dpi=cfg.dpi,
    )
    _save(p)

    return saved

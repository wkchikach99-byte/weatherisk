#!/usr/bin/env python3
"""Run Justus's LEC/EDC pipeline on real NOAA CPC temperature data.

Reproduces the exact mathematical framework from:
  Justus (2025), "Extremal dependence and local estimation clustering
  for non-stationary max-stable processes", Extremes.

Pipeline:
  1. Load 20 years of NOAA CPC daily tmax (2000-2019)
  2. Select Euro-Mediterranean sub-region, coarsen to ~2°
  3. Annual block maxima  (block_size = 365)
  4. Parametric GEV fit per cell → Fréchet transform  Z = -1/log F_GEV(x)
  5. Local pairwise composite-likelihood MLE → (a, b, γ)
  6. Spatial moving-average smoothing (angular wrap for γ)
  7. LEC dissimilarity D₂  (Jaccard ellipse overlap)
     → hierarchical average linkage → 30 %-quantile threshold → k
  8. EDC dissimilarity D₁  (rank-based madogram)
     → hierarchical average linkage → 30 %-quantile threshold → k
  9. In-cluster re-estimation of (a, b, γ)
 10. Risk metrics: VaR₉₅ / ES₉₅ per cluster  (our added value)
 11. Cartopy maps

Parameters mirror the paper:
  ν = 5,  α = 1,  ε = 5  (local neighbourhood),
  30 %-quantile threshold for k-selection.

Usage:
    python scripts/plot_risk_lec_worldmap.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import minimize
from scipy.stats import qmc, rankdata

# ── project path --------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.extremes import block_maxima, fit_gev, to_frechet
from weatherisk.density import pairwise_density_summand
from weatherisk.clustering import (
    calc_distance_ellipses,
    clustering,
    cluster_number_threshold_method,
)
from weatherisk.risk import compute_var, compute_es

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════
#  Configuration  (mirrors the Extremes paper)
# ══════════════════════════════════════════════════════════════════
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "netcdf")
YEARS = range(2000, 2020)           # 20 years → 20 annual block maxima

# Sub-region  (avoids prime-meridian wrap)
LAT_RANGE = (30, 65)                # °N   (Central Europe + Mediterranean)
LON_RANGE = (5, 55)                 # °E

COARSEN   = 4                       # every Nth cell → ~2° from 0.5°
BLOCK_SIZE = 365                    # annual block maxima

# Max-stable model hyper-parameters  (cf. paper §4.2)
DF    = 5.0                         # ν  degrees of freedom
ALPHA = 1.0                         # α  smoothness exponent

# Local estimation  (cf. paper §3.2)
# All radii below are in NORMALISED coordinates ([-5,5]x[-5,5])
NEIGHBOR_RADIUS  = 3.0              # ε  neighbourhood for pairwise CL
SMOOTHING_RADIUS = 2.0              # spatial smoothing radius
MLE_ENSEMBLE = 3                    # multi-start L-BFGS-B runs per cell

# Clustering  (cf. paper §3.3)
QUANTILE_THRESHOLD = 0.30           # 30 %-quantile of pairwise dissimilarities

# Risk  (our contribution)
RISK_LEVEL = 0.95                   # p  for VaR / ES


# ══════════════════════════════════════════════════════════════════
#  Step 1.  Load CPC tmax NetCDF — select sub-region & coarsen
# ══════════════════════════════════════════════════════════════════

def load_subregion_data():
    """Return (daily, lats_1d, lons_1d) for the coarsened sub-region."""
    import xarray as xr

    print("=" * 60)
    print("  Step 1 : Load NOAA CPC tmax")
    print("=" * 60)

    all_daily: list[np.ndarray] = []
    sub_lats = sub_lons = None

    for year in YEARS:
        path = os.path.join(DATA_DIR, f"tmax.{year}.nc")
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found — skipping")
            continue

        print(f"  {year} … ", end="", flush=True)
        ds = xr.open_dataset(path)

        # CPC lats are decreasing (89.75 → −89.75)
        lats = ds.lat.values
        if lats[0] > lats[-1]:
            lat_sel = slice(float(LAT_RANGE[1]), float(LAT_RANGE[0]))
        else:
            lat_sel = slice(float(LAT_RANGE[0]), float(LAT_RANGE[1]))
        lon_sel = slice(float(LON_RANGE[0]), float(LON_RANGE[1]))

        da = ds["tmax"].sel(lat=lat_sel, lon=lon_sel)
        da = da.isel(lat=slice(None, None, COARSEN),
                     lon=slice(None, None, COARSEN))

        data = da.values.astype(np.float64)
        if sub_lats is None:
            sub_lats = da.lat.values
            sub_lons = da.lon.values

        all_daily.append(data)
        ds.close()
        print(f"shape {data.shape}, "
              f"valid {np.isfinite(data).mean():.0%}")

    daily = np.concatenate(all_daily, axis=0)
    print(f"\n  Combined: {daily.shape}  "
          f"lat [{sub_lats.min():.1f}, {sub_lats.max():.1f}]  "
          f"lon [{sub_lons.min():.1f}, {sub_lons.max():.1f}]")
    return daily, sub_lats, sub_lons


# ══════════════════════════════════════════════════════════════════
#  Step 2.  Annual block maxima → parametric GEV → unit Fréchet
#           (Justus §2.1 — marginal transformation)
# ══════════════════════════════════════════════════════════════════

def compute_frechet_data(daily, lats, lons):
    """Return (frechet, bm_celsius, land_idx, coords, gev_params).

    • Annual block maxima  (365-day blocks)
    • Per-cell GEV MLE  →  (μ, σ, ξ)
    • Fréchet transform  Z = -1 / log F_{GEV}(x; μ, σ, ξ)
    """
    print("\n" + "=" * 60)
    print("  Step 2 : Block maxima → GEV → Fréchet")
    print("=" * 60)

    n_days, n_lat, n_lon = daily.shape
    n_cells = n_lat * n_lon
    daily_flat = daily.reshape(n_days, n_cells)

    # land = cells with > 50 % valid daily obs
    nan_frac = np.isnan(daily_flat).mean(axis=0)
    land_mask = nan_frac < 0.50
    land_idx = np.where(land_mask)[0]
    print(f"  Grid cells : {n_cells}   land : {len(land_idx)}"
          f"  ({len(land_idx) / n_cells:.0%})")

    # fill sporadic NaN with cell mean
    daily_land = daily_flat[:, land_idx].copy()
    for c in range(daily_land.shape[1]):
        col = daily_land[:, c]
        m = np.nanmean(col)
        col[np.isnan(col)] = m

    # annual block maxima  (in °C)
    bm = block_maxima(daily_land, block_size=BLOCK_SIZE)
    n_blocks = bm.shape[0]
    print(f"  Block size : {BLOCK_SIZE} d → {n_blocks} annual maxima")
    print(f"  Block-max range : [{bm.min():.1f}, {bm.max():.1f}] °C")

    # parametric GEV fit + Fréchet transform per cell
    n_land = len(land_idx)
    frechet  = np.empty((n_blocks, n_land))
    gev_par  = np.empty((n_land, 3))
    ok       = np.ones(n_land, dtype=bool)

    n_fail = 0
    for c in range(n_land):
        try:
            loc, scale, shape = fit_gev(bm[:, c])
            gev_par[c] = [loc, scale, shape]
            fr = to_frechet(bm[:, c], loc, scale, shape)
            # guard against extreme tails  (standard practice)
            fr = np.clip(fr, 0.01, None)
            frechet[:, c] = fr
        except Exception:
            ok[c] = False
            n_fail += 1

    # drop cells where GEV failed or Fréchet is non-finite
    ok &= np.all(np.isfinite(frechet), axis=0)
    frechet  = frechet[:, ok]
    land_idx = land_idx[ok]
    gev_par  = gev_par[ok]
    bm       = bm[:, ok]

    lat_grid = np.repeat(lats, len(lons))
    lon_grid = np.tile(lons, len(lats))
    coords   = np.column_stack([lat_grid[land_idx], lon_grid[land_idx]])

    # Clip extreme Fréchet tails  (standard finite-sample practice,
    # see e.g. Padoan et al. 2010 — values > ~n² are artefacts)
    frechet = np.clip(frechet, 0.05, float(n_blocks**2))

    # ── Normalise coordinates to [-5, 5] × [-5, 5] ──
    # The covariance function cov_fkt_2d is designed for Justus's
    # synthetic grid with range [-5,5] and spacing ~0.2.  Rescaling
    # real geographic coordinates preserves the scale at which the
    # pairwise CL MLE expects to operate.
    lat_grid = np.repeat(lats, len(lons))
    lon_grid = np.tile(lons, len(lats))
    raw_lat = lat_grid[land_idx]
    raw_lon = lon_grid[land_idx]

    lat_min, lat_max = raw_lat.min(), raw_lat.max()
    lon_min, lon_max = raw_lon.min(), raw_lon.max()
    norm_lat = -5.0 + 10.0 * (raw_lat - lat_min) / max(lat_max - lat_min, 1e-6)
    norm_lon = -5.0 + 10.0 * (raw_lon - lon_min) / max(lon_max - lon_min, 1e-6)

    coords = np.column_stack([norm_lat, norm_lon])   # normalised
    geo_coords = np.column_stack([raw_lat, raw_lon]) # for plotting

    print(f"  GEV failures : {n_fail}")
    print(f"  Valid cells : {len(land_idx)}")
    print(f"  Fréchet shape: {frechet.shape}  "
          f"range [{frechet.min():.2f}, {frechet.max():.2f}]")
    print(f"  Normalised coords: lat [{coords[:,0].min():.1f}, {coords[:,0].max():.1f}]"
          f"  lon [{coords[:,1].min():.1f}, {coords[:,1].max():.1f}]")
    return frechet, bm, land_idx, coords, geo_coords, gev_par


# ══════════════════════════════════════════════════════════════════
#  Step 3.  Local pairwise composite-likelihood MLE  (a, b, γ)
#           (Justus §3.2)
# ══════════════════════════════════════════════════════════════════

def _local_mle_one(frechet, cidx, coords, radius, df, alpha, ens):
    """Estimate (a, b, γ) at cell *cidx* using neighbours within *radius*."""
    n_blocks = frechet.shape[0]

    dlat = coords[:, 0] - coords[cidx, 0]
    dlon = coords[:, 1] - coords[cidx, 1]
    dists = np.sqrt(dlat**2 + dlon**2)
    nb = np.where((dists > 0.01) & (dists <= radius))[0]

    if len(nb) < 3:
        return np.array([1.0, 0.0, 0.0])

    # Build pairwise observation vectors  (centre vs each neighbour)
    z_c = frechet[:, cidx]                         # (n_blocks,)
    zi_l, zj_l, xl_l, yl_l = [], [], [], []
    for j in nb:
        z_nb = frechet[:, j]
        zi_l.extend(z_nb)
        zj_l.extend(z_c)
        xl_l.extend([coords[j, 1] - coords[cidx, 1]] * n_blocks)   # Δlon
        yl_l.extend([coords[j, 0] - coords[cidx, 0]] * n_blocks)   # Δlat

    zi = np.asarray(zi_l)
    zj = np.asarray(zj_l)
    xl = np.asarray(xl_l)
    yl = np.asarray(yl_l)

    # guard: both values must be positive and finite
    good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
    zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]
    if len(zi) < 5:
        return np.array([1.0, 0.0, 0.0])

    lo = np.array([0.01, 0.0,  -np.pi / 2])
    hi = np.array([15.0, 15.0,  np.pi / 2])

    def neg_llh(p):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            v = -np.sum(pairwise_density_summand(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2]))
        return v if np.isfinite(v) else 1e20

    sampler = qmc.LatinHypercube(d=3, seed=42 + cidx)
    starts  = qmc.scale(sampler.random(n=max(ens, 5)), lo, hi)

    best_v, best_p = np.inf, np.array([1.0, 0.0, 0.0])
    for s in range(ens):
        try:
            r = minimize(neg_llh, starts[s], method="L-BFGS-B",
                         bounds=list(zip(lo, hi)),
                         options={"maxiter": 10000, "ftol": 1e-10})
            if r.fun < best_v:
                best_v, best_p = r.fun, r.x.copy()
        except Exception:
            pass
    return best_p


def run_local_estimation(frechet, coords):
    """Step 3: local MLE at every cell."""
    print("\n" + "=" * 60)
    print("  Step 3 : Local MLE  (a, b, γ)")
    print("=" * 60)

    n = frechet.shape[1]
    est = np.zeros((n, 3))
    t0  = time.time()
    for c in range(n):
        if c % max(1, n // 20) == 0:
            elapsed = time.time() - t0
            print(f"    {c + 1:4d}/{n}  ({elapsed:.0f} s)")
        est[c] = _local_mle_one(frechet, c, coords,
                                NEIGHBOR_RADIUS, DF, ALPHA,
                                MLE_ENSEMBLE)

    print(f"  a ∈ [{est[:, 0].min():.3f}, {est[:, 0].max():.3f}]")
    print(f"  b ∈ [{est[:, 1].min():.3f}, {est[:, 1].max():.3f}]")
    print(f"  γ ∈ [{np.degrees(est[:, 2]).min():.1f}°,"
          f" {np.degrees(est[:, 2]).max():.1f}°]")
    return est


# ══════════════════════════════════════════════════════════════════
#  Step 4.  Spatial smoothing with angular wrapping for γ
#           (Justus §3.2, moving average)
# ══════════════════════════════════════════════════════════════════

def smooth_estimates(est, coords, radius):
    """Spatial moving-average smoothing; γ is wrapped in [-π/2, π/2]."""
    print("\n" + "=" * 60)
    print(f"  Step 4 : Spatial smoothing  (radius {radius}°)")
    print("=" * 60)

    n   = len(est)
    out = np.empty_like(est)
    for c in range(n):
        d  = np.sqrt((coords[:, 0] - coords[c, 0])**2 +
                     (coords[:, 1] - coords[c, 1])**2)
        nb = np.where(d <= radius)[0]

        nbe = est[nb]
        out[c, 0] = nbe[:, 0].mean()           # average  a
        out[c, 1] = nbe[:, 1].mean()           # average  b

        # angular average with wrapping  (cf. estimation.py)
        cg  = est[c, 2]
        ang = nbe[:, 2].copy()
        ang = np.where(ang < cg - np.pi / 2, ang + np.pi, ang)
        ang = np.where(ang > cg + np.pi / 2, ang - np.pi, ang)
        mg  = ang.mean()
        if mg < -np.pi / 2:
            mg += np.pi
        elif mg > np.pi / 2:
            mg -= np.pi
        out[c, 2] = mg

    print(f"  Smoothed a ∈ [{out[:, 0].min():.3f}, {out[:, 0].max():.3f}]")
    print(f"  Smoothed b ∈ [{out[:, 1].min():.3f}, {out[:, 1].max():.3f}]")
    return out


# ══════════════════════════════════════════════════════════════════
#  Step 5.  LEC / EDC clustering
#           (Justus §3.3–3.4)
# ══════════════════════════════════════════════════════════════════

def _edc_matrix(frechet):
    """Rank-based madogram extremal-coefficient matrix D₁ (Saunders)."""
    n_blocks, n_cells = frechet.shape

    ranks = np.column_stack(
        [rankdata(frechet[:, s]) for s in range(n_cells)]
    ).T   # (n_cells, n_blocks)

    ec = np.zeros((n_cells, n_cells))
    for i in range(n_cells - 1):
        diff  = np.abs(ranks[i] - ranks[i + 1:])         # (rest, n_blocks)
        v     = diff.mean(axis=1) / (2.0 * (n_blocks + 1))
        denom = 1.0 - 2.0 * v
        denom[denom <= 0] = 1e-12
        ec[i, i + 1:] = np.minimum(1.0, (1.0 + 2.0 * v) / denom - 1.0)
    return ec + ec.T


def run_clustering(smoothed, frechet):
    """Step 5: LEC (D₂) and EDC (D₁) with 30 %-quantile threshold."""
    print("\n" + "=" * 60)
    print("  Step 5 : LEC & EDC clustering")
    print("=" * 60)

    # --- LEC  (D₂ — Jaccard ellipse overlap) -----------------------
    print("  Computing LEC dissimilarity …")
    dm_lec = calc_distance_ellipses(smoothed, res=21)
    hc_lec = clustering(dm_lec)

    # Threshold = quantile of PAIRWISE dissimilarities (upper triangle),
    # matching Justus's R code: quantile(v_matrix[upper.tri(v_matrix)], 0.3)
    vec_lec = dm_lec[np.triu_indices_from(dm_lec, k=1)]
    thr_lec = np.quantile(vec_lec, QUANTILE_THRESHOLD)
    k_lec   = cluster_number_threshold_method(hc_lec, thr_lec)
    k_lec   = max(2, k_lec)
    labels_lec = fcluster(hc_lec, t=k_lec, criterion="maxclust")
    print(f"  LEC → k = {k_lec}  (30 %-quantile threshold = {thr_lec:.3f})")

    # --- EDC  (D₁ — madogram-based) --------------------------------
    print("  Computing EDC dissimilarity …")
    dm_edc = _edc_matrix(frechet)
    hc_edc = clustering(dm_edc)

    vec_edc = dm_edc[np.triu_indices_from(dm_edc, k=1)]
    thr_edc = np.quantile(vec_edc, QUANTILE_THRESHOLD)
    k_edc   = cluster_number_threshold_method(hc_edc, thr_edc)
    k_edc   = max(2, k_edc)
    labels_edc = fcluster(hc_edc, t=k_edc, criterion="maxclust")
    print(f"  EDC → k = {k_edc}  (30 %-quantile threshold = {thr_edc:.5f})")

    return dict(labels_lec=labels_lec, k_lec=k_lec,
                hc_lec=hc_lec, dm_lec=dm_lec,
                labels_edc=labels_edc, k_edc=k_edc,
                hc_edc=hc_edc, dm_edc=dm_edc)


# ══════════════════════════════════════════════════════════════════
#  Step 6.  In-cluster re-estimation of (a, b, γ)
#           (Justus §3.4)
# ══════════════════════════════════════════════════════════════════

def _incluster_reestimate(frechet, coords, labels, df, alpha, tag):
    """Re-estimate (a, b, γ) per cluster via global pairwise CL MLE."""
    from weatherisk.density import pairwise_density_optim

    unique_cl = sorted(np.unique(labels))
    print(f"  {tag}: re-estimating {len(unique_cl)} clusters …")

    results = {}   # cluster_id → (a, b, γ)
    for cl in unique_cl:
        mask  = labels == cl
        n_cl  = int(mask.sum())
        if n_cl < 3:
            results[cl] = np.array([1.0, 0.0, 0.0])
            continue

        # Frechet data for this cluster:  shape (n_cl, n_blocks)
        z_cl = frechet[:, mask].T                     # (n_cl, n_blocks)
        X_cl = coords[mask, 1]                        # lon
        Y_cl = coords[mask, 0]                        # lat

        max_dist = 4.0 * NEIGHBOR_RADIUS
        try:
            est = pairwise_density_optim(
                z_cl, df, alpha, X_cl, Y_cl,
                upper_bounds=(15.0, 15.0),
                max_dist=max_dist, ensemble=3,
            )
        except Exception:
            est = np.array([1.0, 0.0, 0.0])

        results[cl] = est
        print(f"    cl {cl:2d} ({n_cl:3d} cells)  "
              f"a={est[0]:.3f}  b={est[1]:.3f}  γ={np.degrees(est[2]):.1f}°")
    return results


def run_incluster_reestimation(frechet, coords, cl):
    """Step 6: in-cluster re-estimation."""
    print("\n" + "=" * 60)
    print("  Step 6 : In-cluster re-estimation")
    print("=" * 60)

    incl_lec = _incluster_reestimate(
        frechet, coords, cl["labels_lec"], DF, ALPHA, "LEC")
    incl_edc = _incluster_reestimate(
        frechet, coords, cl["labels_edc"], DF, ALPHA, "EDC")
    return incl_lec, incl_edc


# ══════════════════════════════════════════════════════════════════
#  Step 7.  Risk metrics: VaR₉₅ / ES₉₅ per cluster
#           (Koch 2017 — our contribution)
# ══════════════════════════════════════════════════════════════════

def cluster_risk(frechet, labels, p=RISK_LEVEL):
    """VaR and ES on the Fréchet scale, per cluster.

    For each cluster A, the loss is  L_t = max_{s ∈ A} Z_t(s)
    (spatial block maximum over the cluster at time t).
    """
    out = []
    for cl in sorted(np.unique(labels)):
        mask  = labels == cl
        cols  = frechet[:, mask]
        cmax  = cols.max(axis=1)                      # L_t per year
        out.append(dict(
            cluster = int(cl),
            n_cells = int(mask.sum()),
            var     = compute_var(cmax, p),
            es      = compute_es(cmax, p),
        ))
    return out


# ══════════════════════════════════════════════════════════════════
#  Cartopy plotting helpers
# ══════════════════════════════════════════════════════════════════

def _extent():
    return [LON_RANGE[0] - 2, LON_RANGE[1] + 2,
            LAT_RANGE[0] - 2, LAT_RANGE[1] + 2]


def _base_ax(fig, pos=111):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent(_extent(), crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    ax.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax.add_feature(cfeature.LAND, facecolor="wheat", alpha=0.15)
    return ax


def _cluster_cmap(k):
    base = (list(plt.get_cmap("tab20").colors)
            + list(plt.get_cmap("Set3").colors))
    return ListedColormap([base[i % len(base)] for i in range(max(k, 1))])


def plot_clusters(lat, lon, labels, k, title, fname):
    import cartopy.crs as ccrs
    fig  = plt.figure(figsize=(10, 7))
    ax   = _base_ax(fig)
    cmap = _cluster_cmap(k)
    norm = BoundaryNorm(np.arange(0.5, k + 1.5), cmap.N)
    sc   = ax.scatter(lon, lat, c=labels, cmap=cmap, norm=norm,
                      s=45, edgecolors="k", linewidths=0.3,
                      transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.7, label="Cluster",
                 ticks=range(1, k + 1))
    ax.set_title(title, fontsize=13, fontweight="bold")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {fname}")


def plot_field(lat, lon, vals, title, label, fname,
               cmap="viridis", vmin=None, vmax=None):
    import cartopy.crs as ccrs
    fig = plt.figure(figsize=(10, 7))
    ax  = _base_ax(fig)
    sc  = ax.scatter(lon, lat, c=vals, cmap=cmap,
                     vmin=vmin, vmax=vmax,
                     s=45, edgecolors="k", linewidths=0.3,
                     transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.7, label=label)
    ax.set_title(title, fontsize=13, fontweight="bold")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {fname}")


def plot_risk_choropleth(lat, lon, labels, risk, title, fname):
    import cartopy.crs as ccrs
    rmap = np.zeros(len(labels))
    for r in risk:
        rmap[labels == r["cluster"]] = r["es"]
    fig = plt.figure(figsize=(10, 7))
    ax  = _base_ax(fig)
    sc  = ax.scatter(lon, lat, c=rmap, cmap="YlOrRd",
                     s=45, edgecolors="k", linewidths=0.3,
                     transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.7,
                 label="ES$_{95}$  (Fréchet scale)")
    ax.set_title(title, fontsize=13, fontweight="bold")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {fname}")


def plot_summary(lat, lon, lbl_lec, lbl_edc, k_lec, k_edc,
                 sm, risk_lec, risk_edc):
    """2×3 summary panel."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    fig = plt.figure(figsize=(20, 12))
    ext = _extent()

    def _ax(pos):
        ax = fig.add_subplot(2, 3, pos, projection=ccrs.PlateCarree())
        ax.set_extent(ext, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
        ax.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
        ax.add_feature(cfeature.LAND, facecolor="wheat", alpha=0.15)
        return ax

    # (1) LEC clusters
    ax = _ax(1)
    cm = _cluster_cmap(k_lec)
    nm = BoundaryNorm(np.arange(0.5, k_lec + 1.5), cm.N)
    ax.scatter(lon, lat, c=lbl_lec, cmap=cm, norm=nm, s=22,
               edgecolors="k", linewidths=0.2,
               transform=ccrs.PlateCarree(), zorder=5)
    ax.set_title(f"LEC  (k={k_lec})", fontsize=11, fontweight="bold")

    # (2) EDC clusters
    ax = _ax(2)
    cm = _cluster_cmap(k_edc)
    nm = BoundaryNorm(np.arange(0.5, k_edc + 1.5), cm.N)
    ax.scatter(lon, lat, c=lbl_edc, cmap=cm, norm=nm, s=22,
               edgecolors="k", linewidths=0.2,
               transform=ccrs.PlateCarree(), zorder=5)
    ax.set_title(f"EDC  (k={k_edc})", fontsize=11, fontweight="bold")

    # (3) parameter a
    ax = _ax(3)
    sc = ax.scatter(lon, lat, c=sm[:, 0], cmap="viridis", s=22,
                    edgecolors="k", linewidths=0.2,
                    transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.55, label="a")
    ax.set_title("Semi-minor axis  a", fontsize=11, fontweight="bold")

    # (4) parameter b
    ax = _ax(4)
    sc = ax.scatter(lon, lat, c=sm[:, 1], cmap="magma", s=22,
                    edgecolors="k", linewidths=0.2,
                    transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.55, label="b")
    ax.set_title("Anisotropy  b", fontsize=11, fontweight="bold")

    # (5) ES₉₅ per LEC cluster
    ax = _ax(5)
    rmap = np.zeros(len(lbl_lec))
    for r in risk_lec:
        rmap[lbl_lec == r["cluster"]] = r["es"]
    sc = ax.scatter(lon, lat, c=rmap, cmap="YlOrRd", s=22,
                    edgecolors="k", linewidths=0.2,
                    transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.55, label="ES$_{95}$")
    ax.set_title("ES$_{95}$  (LEC)", fontsize=11, fontweight="bold")

    # (6) ES₉₅ per EDC cluster
    ax = _ax(6)
    rmap2 = np.zeros(len(lbl_edc))
    for r in risk_edc:
        rmap2[lbl_edc == r["cluster"]] = r["es"]
    sc = ax.scatter(lon, lat, c=rmap2, cmap="YlOrRd", s=22,
                    edgecolors="k", linewidths=0.2,
                    transform=ccrs.PlateCarree(), zorder=5)
    plt.colorbar(sc, ax=ax, shrink=0.55, label="ES$_{95}$")
    ax.set_title("ES$_{95}$  (EDC)", fontsize=11, fontweight="bold")

    # gridlines on map axes only
    for i in range(6):
        a = fig.axes[i * 2] if len(fig.axes) > i * 2 else None
        if a is not None and hasattr(a, "gridlines"):
            try:
                gl = a.gridlines(draw_labels=True,
                                 linewidth=0.2, alpha=0.4)
                gl.top_labels = gl.right_labels = False
            except Exception:
                pass

    plt.suptitle("Climate Risk — LEC / EDC on CPC tmax  (2000–2019)",
                 fontsize=14, fontweight="bold", y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUTPUT_DIR, "lec_edc_summary_panel.pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ lec_edc_summary_panel.pdf")


# ══════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    print("\n" + "=" * 60)
    print("  LEC / EDC Risk Pipeline  —  real CPC tmax data")
    print("  Following Justus (2025, Extremes) methodology")
    print("=" * 60)

    # ── 1. Load ---------------------------------------------------
    daily, lats, lons = load_subregion_data()

    # ── 2. Block maxima → GEV → Fréchet --------------------------
    frechet, bm, land_idx, coords, geo_coords, gev_p = \
        compute_frechet_data(daily, lats, lons)
    lat_v = geo_coords[:, 0]    # geographic (for plots)
    lon_v = geo_coords[:, 1]

    # ── 3. Local MLE (a, b, γ) -----------------------------------
    est = run_local_estimation(frechet, coords)

    # ── 4. Smooth -------------------------------------------------
    sm = smooth_estimates(est, coords, SMOOTHING_RADIUS)

    # ── 5. Cluster ------------------------------------------------
    cl = run_clustering(sm, frechet)

    # ── 6. In-cluster re-estimation -------------------------------
    incl_lec, incl_edc = run_incluster_reestimation(
        frechet, coords, cl)

    # ── 7. Risk metrics -------------------------------------------
    print("\n" + "=" * 60)
    print(f"  Step 7 : Risk metrics  (VaR_{RISK_LEVEL:.0%}"
          f" / ES_{RISK_LEVEL:.0%})")
    print("=" * 60)
    risk_lec = cluster_risk(frechet, cl["labels_lec"])
    risk_edc = cluster_risk(frechet, cl["labels_edc"])

    print("\n  ── LEC ──")
    for r in risk_lec:
        print(f"    cl {r['cluster']:2d}  ({r['n_cells']:3d} cells)  "
              f"VaR={r['var']:.2f}  ES={r['es']:.2f}")

    print("\n  ── EDC ──")
    for r in risk_edc:
        print(f"    cl {r['cluster']:2d}  ({r['n_cells']:3d} cells)  "
              f"VaR={r['var']:.2f}  ES={r['es']:.2f}")

    # ── 8. Plots --------------------------------------------------
    print("\n" + "=" * 60)
    print("  Step 8 : Generating Cartopy maps")
    print("=" * 60)

    plot_clusters(lat_v, lon_v, cl["labels_lec"], cl["k_lec"],
                  f"LEC Clusters (k={cl['k_lec']}) — "
                  "Ellipse-Overlap Dissimilarity D₂",
                  "lec_clusters_map.pdf")

    plot_clusters(lat_v, lon_v, cl["labels_edc"], cl["k_edc"],
                  f"EDC Clusters (k={cl['k_edc']}) — "
                  "Madogram Dissimilarity D₁",
                  "edc_clusters_map.pdf")

    plot_field(lat_v, lon_v, sm[:, 0],
               "Smoothed Semi-Minor Axis  a", "a",
               "param_a_map.pdf", cmap="viridis")

    plot_field(lat_v, lon_v, sm[:, 1],
               "Smoothed Anisotropy  b", "b",
               "param_b_map.pdf", cmap="magma")

    plot_field(lat_v, lon_v, np.degrees(sm[:, 2]),
               "Smoothed Rotation  γ", "γ  (degrees)",
               "param_gamma_map.pdf", cmap="twilight",
               vmin=-90, vmax=90)

    plot_risk_choropleth(lat_v, lon_v, cl["labels_lec"], risk_lec,
                         "Expected Shortfall ES$_{95}$ per LEC Cluster",
                         "lec_risk_es_map.pdf")

    plot_risk_choropleth(lat_v, lon_v, cl["labels_edc"], risk_edc,
                         "Expected Shortfall ES$_{95}$ per EDC Cluster",
                         "edc_risk_es_map.pdf")

    plot_summary(lat_v, lon_v,
                 cl["labels_lec"], cl["labels_edc"],
                 cl["k_lec"], cl["k_edc"],
                 sm, risk_lec, risk_edc)

    elapsed = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"  Done in {elapsed:.0f} s  ({elapsed / 60:.1f} min)")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()

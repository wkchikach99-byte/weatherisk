#!/usr/bin/env python3
"""GEV fit diagnostics for CPC pipeline block maxima.

Produces:
  1. Q-Q plots for a sample of grid cells (GEV theoretical vs empirical quantiles)
  2. Per-cell Anderson-Darling p-values with a summary histogram
  3. GEV shape parameter (ξ) map to check for spatial coherence

Usage:
    python scripts/gev_diagnostics.py
    python scripts/gev_diagnostics.py --n-sample 30 --output-dir docs/figures
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import genextreme, anderson, kstest


def load_block_maxima_and_gev(data_dir: str, year_start: int, year_end: int,
                               lat_range: tuple, lon_range: tuple,
                               coarsen: int, block_size: int):
    """Load CPC data, compute block maxima, fit GEV per cell.

    Returns bm, gev_params, lats, lons, land_idx.
    """
    import xarray as xr

    all_daily = []
    sub_lats = sub_lons = None

    for year in range(year_start, year_end):
        path = os.path.join(data_dir, f"precip.{year}.nc")
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found — skipping")
            continue

        ds = xr.open_dataset(path)
        lats = ds.lat.values
        if lats[0] > lats[-1]:
            lat_sel = slice(float(lat_range[1]), float(lat_range[0]))
        else:
            lat_sel = slice(float(lat_range[0]), float(lat_range[1]))
        lon_sel = slice(float(lon_range[0]), float(lon_range[1]))

        da = ds["precip"].sel(lat=lat_sel, lon=lon_sel)
        da = da.isel(lat=slice(None, None, coarsen),
                     lon=slice(None, None, coarsen))
        data = da.values.astype(np.float64)
        if sub_lats is None:
            sub_lats = da.lat.values
            sub_lons = da.lon.values
        all_daily.append(data)
        ds.close()

    daily = np.concatenate(all_daily, axis=0)
    n_days, n_lat, n_lon = daily.shape
    n_cells = n_lat * n_lon
    daily_flat = daily.reshape(n_days, n_cells)

    # Land mask
    nan_frac = np.isnan(daily_flat).mean(axis=0)
    land_mask = nan_frac < 0.50
    land_idx = np.where(land_mask)[0]
    daily_land = daily_flat[:, land_idx].copy()

    # Fill NaN with column mean
    for c in range(daily_land.shape[1]):
        col = daily_land[:, c]
        m = np.nanmean(col)
        col[np.isnan(col)] = m

    # Block maxima
    from weatherisk.extremes import block_maxima, fit_gev
    bm = block_maxima(daily_land, block_size=block_size)

    # Fit GEV per cell
    n_land = len(land_idx)
    gev_params = np.empty((n_land, 3))  # loc, scale, shape
    ok = np.ones(n_land, dtype=bool)

    for c in range(n_land):
        try:
            loc, scale, shape = fit_gev(bm[:, c])
            gev_params[c] = [loc, scale, shape]
        except Exception:
            ok[c] = False
            gev_params[c] = [np.nan, np.nan, np.nan]

    # Build geo coords
    lat_grid = np.repeat(sub_lats, len(sub_lons))
    lon_grid = np.tile(sub_lons, len(sub_lats))
    geo_lat = lat_grid[land_idx]
    geo_lon = lon_grid[land_idx]

    return bm, gev_params, ok, geo_lat, geo_lon, sub_lats, sub_lons, land_idx, n_lat, n_lon


def qq_plot_grid(bm, gev_params, sample_idx, geo_lat, geo_lon, save_path):
    """Create a grid of Q-Q plots for sampled cells."""
    n = len(sample_idx)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.2 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, cidx in enumerate(sample_idx):
        ax = axes[i]
        data = np.sort(bm[:, cidx])
        n_obs = len(data)
        loc, scale, shape = gev_params[cidx]
        c_scipy = -shape  # scipy sign convention

        # Theoretical quantiles
        probs = (np.arange(1, n_obs + 1) - 0.375) / (n_obs + 0.25)  # Blom plotting positions
        theo = genextreme.ppf(probs, c_scipy, loc=loc, scale=scale)

        ax.scatter(theo, data, s=18, alpha=0.7, edgecolors="k", linewidths=0.3, zorder=3)
        lims = [min(theo.min(), data.min()) * 0.9,
                max(theo.max(), data.max()) * 1.1]
        ax.plot(lims, lims, "r--", lw=1, alpha=0.7)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("GEV theoretical", fontsize=8)
        ax.set_ylabel("Observed", fontsize=8)
        ax.set_title(f"({geo_lat[cidx]:.1f}°N, {geo_lon[cidx]:.1f}°E)\n"
                      f"ξ={shape:.3f}", fontsize=8)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("GEV Q-Q Plots — Sample of CPC Precipitation Cells\n"
                 "(20 annual block maxima, mm/day)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    print(f"  Saved {save_path}")
    plt.close(fig)


def goodness_of_fit_tests(bm, gev_params, ok):
    """Run KS test (with estimated parameters) per cell.

    Returns p-values array (NaN for failed cells).
    """
    n_cells = bm.shape[1]
    pvals = np.full(n_cells, np.nan)

    for c in range(n_cells):
        if not ok[c]:
            continue
        loc, scale, shape = gev_params[c]
        c_scipy = -shape
        data = bm[:, c]
        try:
            stat, pval = kstest(data, "genextreme", args=(c_scipy, loc, scale))
            pvals[c] = pval
        except Exception:
            pass

    return pvals


def plot_pvalue_histogram(pvals, save_path):
    """Histogram of KS p-values across all cells."""
    valid = pvals[np.isfinite(pvals)]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(valid, bins=20, range=(0, 1), edgecolor="black", alpha=0.7,
            color="steelblue")
    ax.axhline(y=len(valid) / 20, color="red", ls="--", lw=1.5,
               label=f"Uniform expectation ({len(valid)/20:.0f})")
    ax.set_xlabel("KS test p-value")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"GEV Goodness-of-Fit: KS Test p-values\n"
                 f"({len(valid)} land cells, 20 annual block maxima)")
    ax.legend()

    n_reject_05 = (valid < 0.05).sum()
    n_reject_10 = (valid < 0.10).sum()
    ax.text(0.98, 0.95,
            f"Reject at α=0.05: {n_reject_05}/{len(valid)} ({100*n_reject_05/len(valid):.1f}%)\n"
            f"Reject at α=0.10: {n_reject_10}/{len(valid)} ({100*n_reject_10/len(valid):.1f}%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))

    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    print(f"  Saved {save_path}")
    plt.close(fig)


def plot_shape_map(gev_params, ok, geo_lat, geo_lon, lats_1d, lons_1d,
                   land_idx, n_lat, n_lon, save_path):
    """Map of GEV shape parameter ξ across the grid."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        has_cartopy = True
    except ImportError:
        has_cartopy = False

    shape_vals = gev_params[:, 2].copy()
    shape_vals[~ok] = np.nan

    if has_cartopy:
        from weatherisk.map_plotting import _to_grid
        grid, lat_e, lon_e = _to_grid(
            shape_vals, lats_1d, lons_1d, n_lat, n_lon, land_idx, fill=np.nan)

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent([lons_1d.min() - 2, lons_1d.max() + 2,
                       lats_1d.min() - 2, lats_1d.max() + 2],
                      crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=":")

        im = ax.pcolormesh(lon_e, lat_e, grid, cmap="RdBu_r",
                           vmin=-0.5, vmax=0.5,
                           transform=ccrs.PlateCarree())
        plt.colorbar(im, ax=ax, label="GEV shape ξ", shrink=0.7)
        ax.set_title("GEV Shape Parameter (ξ) — CPC Precipitation\n"
                      "Annual Block Maxima 2000–2019  (ξ>0: heavy tail, ξ<0: bounded)")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc = ax.scatter(geo_lon, geo_lat, c=shape_vals, cmap="RdBu_r",
                        vmin=-0.5, vmax=0.5, s=30)
        plt.colorbar(sc, ax=ax, label="GEV shape ξ")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("GEV Shape Parameter (ξ) — CPC Precipitation")

    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    print(f"  Saved {save_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="GEV fit diagnostics")
    parser.add_argument("--data-dir", default="data/netcdf")
    parser.add_argument("--output-dir", default="docs/figures")
    parser.add_argument("--n-sample", type=int, default=20,
                        help="Number of cells for Q-Q plots")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("  GEV Fit Diagnostics — CPC Precipitation")
    print("=" * 60)

    # Use same pipeline config as the main CPC pipeline
    bm, gev_params, ok, geo_lat, geo_lon, lats_1d, lons_1d, land_idx, n_lat, n_lon = \
        load_block_maxima_and_gev(
            data_dir=args.data_dir,
            year_start=2000, year_end=2020,
            lat_range=(30.0, 65.0), lon_range=(5.0, 55.0),
            coarsen=4, block_size=365,
        )

    n_land = bm.shape[1]
    n_ok = ok.sum()
    print(f"\n  Land cells: {n_land}  (GEV fit OK: {n_ok})")
    print(f"  Block maxima shape: {bm.shape}  (years × cells)")
    print(f"  GEV shape ξ range: [{gev_params[ok, 2].min():.3f}, "
          f"{gev_params[ok, 2].max():.3f}]")
    print(f"  GEV shape ξ mean:  {gev_params[ok, 2].mean():.3f}")
    print(f"  GEV shape ξ median: {np.median(gev_params[ok, 2]):.3f}")

    # 1. Q-Q plots for a random sample
    rng = np.random.default_rng(args.seed)
    ok_idx = np.where(ok)[0]
    n_sample = min(args.n_sample, len(ok_idx))
    sample = rng.choice(ok_idx, size=n_sample, replace=False)
    sample.sort()

    print(f"\n  Creating Q-Q plots for {n_sample} sampled cells...")
    qq_path = os.path.join(args.output_dir, "gev_qq_plots.png")
    qq_plot_grid(bm, gev_params, sample, geo_lat, geo_lon, qq_path)

    # 2. KS test for all cells
    print("  Running KS goodness-of-fit tests...")
    pvals = goodness_of_fit_tests(bm, gev_params, ok)
    valid_p = pvals[np.isfinite(pvals)]

    n_reject_05 = (valid_p < 0.05).sum()
    n_reject_10 = (valid_p < 0.10).sum()
    print(f"  KS test results ({len(valid_p)} cells):")
    print(f"    Reject at α=0.05: {n_reject_05} ({100*n_reject_05/len(valid_p):.1f}%)")
    print(f"    Reject at α=0.10: {n_reject_10} ({100*n_reject_10/len(valid_p):.1f}%)")
    print(f"    Median p-value: {np.median(valid_p):.3f}")

    hist_path = os.path.join(args.output_dir, "gev_ks_pvalues.png")
    plot_pvalue_histogram(pvals, hist_path)

    # 3. Shape parameter map
    print("  Creating GEV shape parameter map...")
    shape_path = os.path.join(args.output_dir, "gev_shape_map.png")
    plot_shape_map(gev_params, ok, geo_lat, geo_lon,
                   lats_1d, lons_1d, land_idx, n_lat, n_lon, shape_path)

    print("\n  Done.")


if __name__ == "__main__":
    main()

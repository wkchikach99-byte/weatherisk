"""Gridded GDP exposure layer for climate risk weighting.

Loads the Kummu et al. (2018) *Gridded global datasets for Gross
Domestic Product* at 5 arc-min resolution and regrids it to match the
coarsened CPC pipeline grid.

Reference
---------
Kummu M, Taka M, Guillaume JHA (2018).  Gridded global datasets for
Gross Domestic Product and Human Development Index over 1990–2015.
*Scientific Data* 5:180004.  doi:10.1038/sdata.2018.4

Data download: https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0

The GDP values are in constant 2011 international USD (PPP) per grid
cell.  Missing ocean / no-data cells are encoded as ``-9``.
"""

from __future__ import annotations

import os

import numpy as np


# ── public API ─────────────────────────────────────────────────

def load_gdp_grid(
    gdp_path: str,
    *,
    year: int = 2015,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load the Kummu et al. GDP PPP grid for a given year.

    Parameters
    ----------
    gdp_path : str
        Path to ``GDP_PPP_1990_2015_5arcmin_v2.nc``.
    year : int
        Year to extract (1990–2015).  Default 2015 (latest).

    Returns
    -------
    gdp : (n_lat, n_lon) float array
        GDP PPP per cell in constant 2011 USD.  Ocean/missing → 0.
    lats : (n_lat,) array
        Latitude centres (north → south).
    lons : (n_lon,) array
        Longitude centres (west → east).
    """
    import h5py

    with h5py.File(gdp_path, "r") as f:
        times = f["time"][:]
        idx = int(np.argmin(np.abs(times - year)))
        gdp = f["GDP_PPP"][idx].astype(np.float64)
        lats = f["latitude"][:].astype(np.float64)
        lons = f["longitude"][:].astype(np.float64)

    # Replace missing-value sentinel (-9) with 0
    gdp[gdp < 0] = 0.0
    return gdp, lats, lons


def regrid_gdp_to_pipeline(
    gdp: np.ndarray,
    gdp_lats: np.ndarray,
    gdp_lons: np.ndarray,
    pipeline_lats: np.ndarray,
    pipeline_lons: np.ndarray,
) -> np.ndarray:
    """Aggregate GDP onto the coarser CPC pipeline grid.

    For each pipeline cell, sums the GDP of all fine-resolution cells
    whose centres fall within the pipeline cell boundaries.  This gives
    the *total* GDP PPP within each coarse cell.

    Parameters
    ----------
    gdp : (n_lat_fine, n_lon_fine) array
        Fine-resolution GDP grid (5 arc-min).
    gdp_lats, gdp_lons : 1-D arrays
        Coordinates of the fine grid.
    pipeline_lats, pipeline_lons : 1-D arrays
        Coordinates of the coarsened CPC pipeline grid.

    Returns
    -------
    gdp_coarse : (n_pipeline_cells,) array
        Total GDP PPP per coarse pipeline cell (flattened, same order
        as ``np.repeat(lats, len(lons))`` / ``np.tile(lons, len(lats))``).
    """
    # Build coarse cell edges
    def _edges(c):
        c = np.asarray(c, dtype=float)
        d = np.diff(c)
        e = np.empty(len(c) + 1)
        e[1:-1] = c[:-1] + d / 2
        e[0] = c[0] - d[0] / 2
        e[-1] = c[-1] + d[-1] / 2
        return e

    lat_e = _edges(pipeline_lats)
    lon_e = _edges(pipeline_lons)
    n_clat = len(pipeline_lats)
    n_clon = len(pipeline_lons)

    coarse = np.zeros((n_clat, n_clon), dtype=np.float64)

    for i in range(n_clat):
        lat_lo = min(lat_e[i], lat_e[i + 1])
        lat_hi = max(lat_e[i], lat_e[i + 1])
        lat_mask = (gdp_lats >= lat_lo) & (gdp_lats < lat_hi)
        if not lat_mask.any():
            continue
        lat_sl = np.where(lat_mask)[0]
        i0, i1 = lat_sl[0], lat_sl[-1] + 1

        for j in range(n_clon):
            lon_lo = min(lon_e[j], lon_e[j + 1])
            lon_hi = max(lon_e[j], lon_e[j + 1])
            lon_mask = (gdp_lons >= lon_lo) & (gdp_lons < lon_hi)
            if not lon_mask.any():
                continue
            lon_sl = np.where(lon_mask)[0]
            j0, j1 = lon_sl[0], lon_sl[-1] + 1

            coarse[i, j] = gdp[i0:i1, j0:j1].sum()

    return coarse.ravel()


def gdp_for_land_cells(
    gdp_path: str,
    pipeline_lats: np.ndarray,
    pipeline_lons: np.ndarray,
    land_idx: np.ndarray,
    *,
    year: int = 2015,
    verbose: bool = True,
) -> np.ndarray:
    """One-call helper: load GDP, regrid, extract land cells.

    Parameters
    ----------
    gdp_path : str
        Path to Kummu et al. NetCDF file.
    pipeline_lats, pipeline_lons : 1-D arrays
        Coordinate vectors of the coarsened CPC grid.
    land_idx : 1-D int array
        Flat indices of valid land cells within the full grid.
    year : int
        GDP year (default 2015).
    verbose : bool
        Print progress.

    Returns
    -------
    gdp_land : (n_land,) array
        GDP PPP (constant 2011 USD) per valid pipeline cell.
    """
    if verbose:
        print(f"\n  Loading GDP PPP ({year}) from {os.path.basename(gdp_path)}")

    gdp_fine, g_lats, g_lons = load_gdp_grid(gdp_path, year=year)

    if verbose:
        print(f"  Fine grid: {gdp_fine.shape}  "
              f"({(gdp_fine > 0).sum():,} cells with GDP > 0)")
        print(f"  Regridding to pipeline grid "
              f"({len(pipeline_lats)}×{len(pipeline_lons)}) …")

    gdp_coarse = regrid_gdp_to_pipeline(
        gdp_fine, g_lats, g_lons, pipeline_lats, pipeline_lons,
    )
    gdp_land = gdp_coarse[land_idx]

    if verbose:
        nonzero = (gdp_land > 0).sum()
        total = gdp_land.sum()
        print(f"  Land cells with GDP > 0: {nonzero} / {len(gdp_land)}")
        print(f"  Total GDP in region: ${total:,.0f}")

    return gdp_land

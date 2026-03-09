"""CMIP6 data discovery, download, and loading for AWI-ESM-1-1-LR.

Handles automatic checking and downloading of AWI-ESM-1-1-LR monthly
precipitation data from ESGF, and loading/concatenation of multi-file
NetCDF datasets.

The data is used to reproduce Figure 9 from:
    Contzen et al. (2025), *Extremal dependence and local estimation
    clustering for non-stationary max-stable processes*, Extremes 28:713–737.

Model details (from paper §5):
    - AWI-ESM-1-1-LR, fully coupled global climate model
    - Atmospheric component: ECHAM6 on T63 grid (192×96, ~1.85°×1.85°)
    - Ocean component: FESOM (variable resolution)
    - Experiment: historical (1850–2005)
    - Variable: pr (monthly mean precipitation, kg m⁻² s⁻¹)
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np


# ── Default paths ──────────────────────────────────────────────

DEFAULT_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "data", "cmip6"
)

# Common locations on AWI Albedo / DKRZ where CMIP6 data is pre-staged
HPC_SEARCH_PATHS = [
    "/pool/data/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical",
    "/work/bk1099",
    "/work/ik1017",
    "/pool/data/CMIP6",
    "/scratch/{user}/cmip6",
]

# ESGF download info
ESGF_BASE_URL = (
    "https://esgf-data.dkrz.de/thredds/fileServer/cmip6/CMIP/"
    "AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212/"
)

# Expected file pattern
FILE_PATTERN = "pr_Amon_AWI-ESM-1-1-LR_historical_r1i1p1f1_gn_"

# Expected year chunks (CMIP6 standard splitting)
EXPECTED_CHUNKS = [
    "185001-185912",
    "186001-186912",
    "187001-187912",
    "188001-188912",
    "189001-189912",
    "190001-190912",
    "191001-191912",
    "192001-192912",
    "193001-193912",
    "194001-194912",
    "195001-195912",
    "196001-196912",
    "197001-197912",
    "198001-198912",
    "199001-199912",
    "200001-200512",
]


def _log(msg: str, verbose: bool = True) -> None:
    if verbose:
        print(msg, flush=True)


# ══════════════════════════════════════════════════════════════════
#  Data discovery
# ══════════════════════════════════════════════════════════════════

def find_cmip6_files(
    data_dir: str = DEFAULT_DATA_DIR,
    *,
    verbose: bool = True,
) -> list[str]:
    """Find AWI-ESM-1-1-LR monthly precipitation NetCDF files.

    Searches *data_dir* for files matching the expected naming convention.
    Returns sorted list of absolute paths.
    """
    data_path = Path(data_dir)
    if not data_path.exists():
        return []

    files = sorted(data_path.glob(f"{FILE_PATTERN}*.nc"))
    if files and verbose:
        _log(f"  Found {len(files)} AWI-ESM-1-1-LR precipitation file(s) "
             f"in {data_dir}")
    return [str(f) for f in files]


def find_on_hpc(*, verbose: bool = True) -> Optional[str]:
    """Search common HPC paths for pre-staged CMIP6 data.

    Returns the directory containing the files, or None.
    """
    user = os.environ.get("USER", "unknown")

    for pattern in HPC_SEARCH_PATHS:
        search_dir = pattern.format(user=user)
        if not os.path.isdir(search_dir):
            continue

        if verbose:
            _log(f"  Searching {search_dir} ...")

        # Direct check for the Amon/pr subdirectory
        for root, dirs, fnames in os.walk(search_dir):
            nc_files = [f for f in fnames if f.startswith(FILE_PATTERN)]
            if nc_files:
                _log(f"  ✅ Found {len(nc_files)} file(s) in {root}")
                return root

        # Also try find (faster for deep trees)
        try:
            result = subprocess.run(
                ["find", search_dir, "-name", f"{FILE_PATTERN}*.nc",
                 "-type", "f", "-maxdepth", 8],
                capture_output=True, text=True, timeout=30,
            )
            if result.stdout.strip():
                first = result.stdout.strip().split("\n")[0]
                found_dir = os.path.dirname(first)
                n = len(result.stdout.strip().split("\n"))
                _log(f"  ✅ Found {n} file(s) via find in {found_dir}")
                return found_dir
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

    return None


# ══════════════════════════════════════════════════════════════════
#  Download from ESGF
# ══════════════════════════════════════════════════════════════════

def download_from_esgf(
    data_dir: str = DEFAULT_DATA_DIR,
    *,
    verbose: bool = True,
) -> bool:
    """Download AWI-ESM-1-1-LR monthly precipitation from ESGF (DKRZ node).

    Uses wget or curl.  Supports resume (-c flag).

    Returns True if all files are present after download.
    """
    os.makedirs(data_dir, exist_ok=True)

    _log("\n📥 Downloading AWI-ESM-1-1-LR monthly precipitation from ESGF …",
         verbose)
    _log(f"   Source: DKRZ ESGF node", verbose)
    _log(f"   Target: {data_dir}\n", verbose)

    success = True
    for chunk in EXPECTED_CHUNKS:
        fname = f"{FILE_PATTERN}{chunk}.nc"
        target = os.path.join(data_dir, fname)

        if os.path.exists(target):
            _log(f"   ⏭️  {fname} (exists)", verbose)
            continue

        url = ESGF_BASE_URL + fname
        _log(f"   ⬇️  {fname} …", verbose)

        # Try wget first (more common on HPC), then curl
        downloaded = False
        for tool, cmd in [
            ("wget", ["wget", "-c", "-q", "--show-progress",
                       "-P", data_dir, url]),
            ("curl", ["curl", "-C", "-", "-L", "-o", target, url]),
        ]:
            try:
                subprocess.run(cmd, check=True, timeout=7200)
                downloaded = True
                break
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue
            except subprocess.TimeoutExpired:
                _log(f"   ⚠️  Timeout downloading {fname} with {tool}", verbose)

        if not downloaded:
            _log(f"   ❌ Failed to download {fname}", verbose)
            success = False

    return success


def symlink_hpc_data(
    source_dir: str,
    data_dir: str = DEFAULT_DATA_DIR,
    *,
    verbose: bool = True,
) -> None:
    """Create symlinks from HPC pre-staged data into our data directory."""
    os.makedirs(data_dir, exist_ok=True)

    for fname in os.listdir(source_dir):
        if fname.startswith(FILE_PATTERN) and fname.endswith(".nc"):
            src = os.path.join(source_dir, fname)
            dst = os.path.join(data_dir, fname)
            if not os.path.exists(dst):
                os.symlink(src, dst)
                _log(f"   🔗 {fname}", verbose)


# ══════════════════════════════════════════════════════════════════
#  Ensure data exists (main entry point)
# ══════════════════════════════════════════════════════════════════

def ensure_cmip6_data(
    data_dir: str = DEFAULT_DATA_DIR,
    *,
    allow_download: bool = True,
    verbose: bool = True,
) -> list[str]:
    """Check for data, discover on HPC, or download from ESGF.

    This is the main entry point.  Call this before any pipeline run.

    Parameters
    ----------
    data_dir : str
        Directory where data files should be (or will be downloaded to).
    allow_download : bool
        If True and data is not found locally or on HPC, attempt ESGF
        download.
    verbose : bool
        Print progress messages.

    Returns
    -------
    list of str
        Sorted paths to all found NetCDF files.

    Raises
    ------
    FileNotFoundError
        If no data could be found or downloaded.
    """
    _log("=" * 60, verbose)
    _log("  Checking for AWI-ESM-1-1-LR precipitation data …", verbose)
    _log("=" * 60, verbose)

    # 1. Check local directory
    files = find_cmip6_files(data_dir, verbose=verbose)
    if files:
        _log(f"  ✅ {len(files)} file(s) ready in {data_dir}", verbose)
        return files

    # 2. Search HPC paths
    _log("\n  Data not in local dir. Searching HPC paths …", verbose)
    hpc_dir = find_on_hpc(verbose=verbose)
    if hpc_dir:
        _log(f"\n  Creating symlinks from {hpc_dir} …", verbose)
        symlink_hpc_data(hpc_dir, data_dir, verbose=verbose)
        files = find_cmip6_files(data_dir, verbose=verbose)
        if files:
            return files

    # 3. Download from ESGF
    if allow_download:
        _log("\n  Not found on HPC. Downloading from ESGF …", verbose)
        download_from_esgf(data_dir, verbose=verbose)
        files = find_cmip6_files(data_dir, verbose=verbose)
        if files:
            return files

    # 4. Give up with helpful message
    msg = (
        f"AWI-ESM-1-1-LR precipitation data not found in {data_dir}.\n\n"
        "Manual download instructions:\n"
        "  1. Go to https://esgf-data.dkrz.de/search/cmip6-dkrz/\n"
        "  2. Search: source_id=AWI-ESM-1-1-LR, experiment_id=historical,\n"
        "     variable_id=pr, table_id=Amon, variant_label=r1i1p1f1\n"
        "  3. Download all files to: {data_dir}\n\n"
        "On Albedo HPC, the data may be at:\n"
        "  /pool/data/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/"
        "r1i1p1f1/Amon/pr/gn/\n"
    )
    raise FileNotFoundError(msg)


# ══════════════════════════════════════════════════════════════════
#  Load and concatenate
# ══════════════════════════════════════════════════════════════════

def load_monthly_precipitation(
    data_dir: str = DEFAULT_DATA_DIR,
    year_start: int = 1850,
    year_end: int = 2005,
    *,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load and concatenate monthly precipitation from AWI-ESM-1-1-LR.

    Parameters
    ----------
    data_dir : str
        Directory containing the NetCDF files.
    year_start, year_end : int
        Year range (inclusive on both ends).
    verbose : bool
        Print progress.

    Returns
    -------
    pr : ndarray, shape (n_months, n_lat, n_lon)
        Monthly mean precipitation rate (kg m⁻² s⁻¹).
    times : ndarray
        Time coordinate (datetime64).
    lats : ndarray, shape (n_lat,)
        Latitude values.
    lons : ndarray, shape (n_lon,)
        Longitude values.
    """
    import xarray as xr

    files = ensure_cmip6_data(data_dir, verbose=verbose)

    _log(f"\n  Loading {len(files)} file(s) …", verbose)
    ds = xr.open_mfdataset(files, combine="by_coords")

    # Select time range
    time_sel = slice(f"{year_start}-01-01", f"{year_end}-12-31")
    ds = ds.sel(time=time_sel)

    pr = ds["pr"].values.astype(np.float64)
    times = ds["time"].values
    lats = ds["lat"].values
    lons = ds["lon"].values

    ds.close()

    _log(f"  Loaded: shape={pr.shape}, "
         f"lat=[{lats.min():.1f}, {lats.max():.1f}], "
         f"lon=[{lons.min():.1f}, {lons.max():.1f}]", verbose)
    _log(f"  Time: {str(times[0])[:7]} → {str(times[-1])[:7]} "
         f"({len(times)} months)", verbose)

    # Verify expected grid size
    if pr.shape[1:] != (96, 192):
        _log(f"  ⚠️  Expected T63 grid (96×192), got {pr.shape[1:]}", verbose)

    return pr, times, lats, lons

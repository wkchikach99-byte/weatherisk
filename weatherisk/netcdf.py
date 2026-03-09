"""NetCDF climate data ingestion (CPC, ERA5, CMIP6).

Reads gridded climate data from NetCDF4 files into NumPy arrays,
handling longitude wrapping and missing data.
"""

from __future__ import annotations

import numpy as np


def load_climate_data(
    path: str,
    variable: str = "tmax",
) -> np.ndarray:
    """Load a climate variable from a NetCDF file.

    Parameters
    ----------
    path : str
        Path to a NetCDF4 file.
    variable : str
        Name of the data variable to extract.

    Returns
    -------
    ndarray, shape (time, lat, lon)
    """
    import xarray as xr

    ds = xr.open_dataset(path)
    da = ds[variable]
    data = da.values.astype(np.float64)
    ds.close()
    return data


def wrap_longitude(lons: np.ndarray) -> np.ndarray:
    """Convert longitudes from [0, 360) to [-180, 180)."""
    return np.where(lons > 180, lons - 360, lons)

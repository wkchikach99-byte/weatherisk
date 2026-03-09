"""Tests for weatherisk.netcdf — NetCDF data ingestion."""

import numpy as np
import pytest


class TestNetCDFReader:
    def test_load_synthetic_netcdf(self, tmp_path):
        """Create a tiny synthetic NetCDF and verify it loads correctly."""
        xr = pytest.importorskip("xarray")

        # Create synthetic CPC-like file
        lats = np.arange(-2, 3, 1.0)
        lons = np.arange(-2, 3, 1.0)
        time = np.arange(365)
        data = np.random.default_rng(0).random((365, 5, 5)).astype(np.float32)
        ds = xr.Dataset(
            {"tmax": (["time", "lat", "lon"], data)},
            coords={"time": time, "lat": lats, "lon": lons},
        )
        nc_path = tmp_path / "test.nc"
        ds.to_netcdf(nc_path)

        from weatherisk.netcdf import load_climate_data

        result = load_climate_data(str(nc_path), variable="tmax")
        assert result.shape == (365, 5, 5)

    def test_lon_wrapping(self):
        from weatherisk.netcdf import wrap_longitude

        lons = np.array([0, 90, 180, 270, 359])
        wrapped = wrap_longitude(lons)
        assert wrapped.min() >= -180
        assert wrapped.max() <= 180

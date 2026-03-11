"""Diagnose which cells fail the MLE NLL comparison for Python backend."""
import numpy as np
import pandas as pd
from pathlib import Path

CMIP6_MINI_REF = Path("tests/reference_data/cmip6_mini")

from weatherisk.cmip6_pipeline import CMIP6Config, _local_mle_one_cmip6
from weatherisk.backend import neg_log_likelihood_sum
import weatherisk.backend as _be

grid_df = pd.read_csv(CMIP6_MINI_REF / "grid_coordinates.csv")
frechet = pd.read_csv(CMIP6_MINI_REF / "frechet.csv").values
r_est = pd.read_csv(CMIP6_MINI_REF / "local_estimates_all.csv")
lec_params = pd.read_csv(CMIP6_MINI_REF / "lec_scalar_params.csv")
abstand = int(
    lec_params.loc[lec_params["parameter"] == "locest_abst", "value"].iloc[0]
)
grid_coords = np.column_stack([grid_df["Y"].values, grid_df["X"].values])

cfg = CMIP6Config(
    neighbor_radius=float(abstand),
    smoothing_radius=2.0,
    df=5.0,
    alpha=1.0,
)

n_grid = len(grid_df)
n_years = frechet.shape[0]
r_est_arr = np.column_stack(
    [r_est["a_est"].values, r_est["b_est"].values, r_est["g_est"].values]
)

# Test with Python backend
_be._USE_RUST = False
py_est = np.zeros((n_grid, 3))
for cidx in range(n_grid):
    py_est[cidx] = _local_mle_one_cmip6(frechet, cidx, grid_coords, cfg)

print("=== Python backend NLL comparison ===")
for cidx in range(n_grid):
    di = grid_coords[:, 0] - grid_coords[cidx, 0]
    dj = grid_coords[:, 1] - grid_coords[cidx, 1]
    dists = np.sqrt(di**2 + dj**2)
    nb = np.where((dists > 0.01) & (dists <= cfg.neighbor_radius))[0]
    if len(nb) < 3:
        continue
    z_c = frechet[:, cidx]
    zi = frechet[:, nb].T.reshape(-1)
    zj = np.tile(z_c, len(nb))
    xl = np.repeat(dj[nb], n_years).astype(np.float64)
    yl = np.repeat(di[nb], n_years).astype(np.float64)
    good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
    zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]
    nll_py = neg_log_likelihood_sum(
        zi, zj, xl, yl, 5.0, 1.0, *py_est[cidx]
    )
    nll_r = neg_log_likelihood_sum(
        zi, zj, xl, yl, 5.0, 1.0, *r_est_arr[cidx]
    )
    excess = nll_py - nll_r
    if excess > 0.01:
        print(
            f"  Cell {cidx}: NLL_py={nll_py:.6f}, NLL_r={nll_r:.6f}, "
            f"excess={excess:.6f}"
        )
        print(f"    py_est={py_est[cidx]}, r_est={r_est_arr[cidx]}")

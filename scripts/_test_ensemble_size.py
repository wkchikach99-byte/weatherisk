"""Test whether more ensemble starts fix the 4 failing cells."""
import numpy as np
import pandas as pd
from pathlib import Path

CMIP6_MINI_REF = Path("tests/reference_data/cmip6_mini")

from weatherisk.backend import neg_log_likelihood_sum, optimize_local_mle
import weatherisk.backend as _be

grid_df = pd.read_csv(CMIP6_MINI_REF / "grid_coordinates.csv")
frechet = pd.read_csv(CMIP6_MINI_REF / "frechet.csv").values
r_est = pd.read_csv(CMIP6_MINI_REF / "local_estimates_all.csv")
lec_params = pd.read_csv(CMIP6_MINI_REF / "lec_scalar_params.csv")
abstand = int(
    lec_params.loc[lec_params["parameter"] == "locest_abst", "value"].iloc[0]
)
grid_coords = np.column_stack([grid_df["Y"].values, grid_df["X"].values])

n_grid = len(grid_df)
n_years = frechet.shape[0]
r_est_arr = np.column_stack(
    [r_est["a_est"].values, r_est["b_est"].values, r_est["g_est"].values]
)

_be._USE_RUST = False

failing_cells = [8, 9, 13, 23]

for cidx in failing_cells:
    di = grid_coords[:, 0] - grid_coords[cidx, 0]
    dj = grid_coords[:, 1] - grid_coords[cidx, 1]
    dists = np.sqrt(di**2 + dj**2)
    nb = np.where((dists > 0.01) & (dists <= float(abstand)))[0]
    z_c = frechet[:, cidx]
    zi = frechet[:, nb].T.reshape(-1)
    zj = np.tile(z_c, len(nb))
    xl = np.repeat(dj[nb], n_years).astype(np.float64)
    yl = np.repeat(di[nb], n_years).astype(np.float64)
    good = (zi > 0) & (zj > 0) & np.isfinite(zi) & np.isfinite(zj)
    zi, zj, xl, yl = zi[good], zj[good], xl[good], yl[good]

    nll_r = neg_log_likelihood_sum(
        zi, zj, xl, yl, 5.0, 1.0, *r_est_arr[cidx]
    )

    for ens in [5, 10, 20, 50]:
        est = optimize_local_mle(zi, zj, xl, yl, 5.0, 1.0, ensemble=ens, seed=42)
        nll_py = neg_log_likelihood_sum(zi, zj, xl, yl, 5.0, 1.0, *est)
        excess = nll_py - nll_r
        status = "PASS" if excess <= 0.1 else "FAIL"
        print(f"  Cell {cidx} ens={ens:3d}: NLL_py={nll_py:.6f}, NLL_r={nll_r:.6f}, excess={excess:+.6f} [{status}]")
    print()

"""R-parity test for: pairwise_density_optim_local (local MLE optimizer)

R function: pairwise_density_optim_local in r_code/functions.R
Python equivalent: weatherisk.backend.optimize_local_mle

Criterion (likelihood, not parameters):
    At every selected grid point, the negative pairwise log-likelihood of
    Python's optimum must not exceed that of R's optimum by more than
    NLL_TOL = 0.03.  Both NLLs are evaluated by the same deterministic
    function on the same frozen pair arrays, so the test isolates the
    quality of the optimum.

Why not parameters:
    R draws starts with maximinLHS and runs optim(L-BFGS-B); Python draws
    starts with scipy LatinHypercube(seed=42) and runs SciPy's L-BFGS-B.
    When a sits at its lower bound the likelihood is flat in (b, g), so the
    two optimizers can stop at different parameters with equal likelihood.
    Measured on this fixture: at grid point 1 Python's optimum is 0.92 NLL
    units better than R's while its parameters differ by 0.72.

Start budget:
    The test uses the package default (CMIP6Config().mle_ensemble), not R's
    locest_ensemble, so it checks the budget a user actually gets.

Fixtures (frozen, from tests/generate_r_reference.R):
    tests/reference_data/local_optimizer_selected_outputs.csv
        grid_point, bounds and R's selected (a_est, b_est, g_est).
    tests/reference_data/local_optimizer_selected_inputs.bin
        int32 n_rows, then float64 arrays grid_point, zi, zj, xl, yl.
"""

import numpy as np
import pandas as pd
from tests.parity.conftest import (
    NLL_TOL,
    REF,
    package_default_ensemble,
    skip_if_no_ref,
)


def _load_selected_inputs():
    # Binary avoids 1-ULP disagreements between R's and pandas' CSV parsers.
    with open(REF / "local_optimizer_selected_inputs.bin", "rb") as f:
        n_rows = int(np.frombuffer(f.read(4), dtype=np.int32)[0])
        arrays = [np.frombuffer(f.read(n_rows * 8), dtype=np.float64) for _ in range(5)]
    return arrays


class TestParityPairwiseDensityOptimLocal:

    def setup_method(self):
        skip_if_no_ref()

    def test_selected_outputs_reach_r_likelihood(self):
        """Python's optimum is at least as likely as R's at every selected point."""
        from weatherisk.backend import neg_log_likelihood_sum, optimize_local_mle

        ref_outputs = pd.read_csv(REF / "local_optimizer_selected_outputs.csv")
        grid_points, zi_all, zj_all, xl_all, yl_all = _load_selected_inputs()
        ensemble = package_default_ensemble()

        excess = {}
        for _, row in ref_outputs.iterrows():
            gp = int(row["grid_point"])
            mask = grid_points == gp
            zi, zj, xl, yl = zi_all[mask], zj_all[mask], xl_all[mask], yl_all[mask]
            assert zi.size > 0, f"no frozen pair inputs for grid point {gp}"

            py_est = optimize_local_mle(
                zi, zj, xl, yl,
                df=5.0,
                alpha=1.0,
                ensemble=ensemble,
                seed=42,
                max_boundary_retries=5,
                lower_bounds=(float(row["lower_a"]), float(row["lower_b"])),
                upper_bounds=(float(row["upper_a"]), float(row["upper_b"])),
            )
            r_est = (row["a_est"], row["b_est"], row["g_est"])
            nll_py = neg_log_likelihood_sum(zi, zj, xl, yl, 5.0, 1.0, *py_est)
            nll_r = neg_log_likelihood_sum(zi, zj, xl, yl, 5.0, 1.0, *r_est)
            assert np.isfinite(nll_py) and np.isfinite(nll_r)
            excess[gp] = nll_py - nll_r

        worse = {gp: d for gp, d in excess.items() if d > NLL_TOL}
        assert not worse, (
            f"optimize_local_mle (ensemble={ensemble}) is worse than R by more "
            f"than {NLL_TOL} NLL units at grid points {worse}"
        )

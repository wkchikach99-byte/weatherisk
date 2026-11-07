"""Sweep the local-MLE start budget against the frozen R references.

For every (start design, ensemble, seed) the script runs
``weatherisk.backend.optimize_local_mle`` at each cell of two R fixtures
and reports, per cell, ``NLL(Python optimum) - NLL(R optimum)`` on the
same pair arrays. A configuration passes when no cell, at any seed,
exceeds ``--tol``.

Fixtures:
    mini      the 25 cells of ``tests/reference_data/cmip6_mini``; the pair
              arrays are captured from ``cmip6_pipeline._local_mle_one_cmip6``
              so the optimizer sees exactly the inputs the pipeline builds.
    selected  the grid points of ``local_optimizer_selected_*`` with their
              own bounds, as in ``tests/parity``.

The start design is switched by replacing ``backend._local_mle_starts``
for the duration of a run, so the script measures the shipped optimizer.
Designs ``sobol`` and ``uniform`` are nested (the first n starts do not
depend on the total drawn), so best-of-n cannot get worse as n grows;
``lhs`` and ``random-cd`` are not nested.

    PATH=/usr/bin:/bin .venv/bin/python scripts/optimizer_start_sweep.py \\
        --ensembles 5 10 20 --designs lhs sobol --seeds 42 0 1 --csv /tmp/sweep.csv
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
from scipy.stats import qmc

from weatherisk import backend
from weatherisk.backend import neg_log_likelihood_sum
from weatherisk.cmip6_pipeline import CMIP6Config, _local_mle_one_cmip6

REF = Path(__file__).resolve().parents[1] / "tests" / "reference_data"


def _lhs(n, seed):
    return qmc.LatinHypercube(d=3, seed=seed).random(n)


def _random_cd(n, seed):
    return qmc.LatinHypercube(d=3, seed=seed, optimization="random-cd").random(n)


def _sobol(n, seed):
    m = int(np.ceil(np.log2(max(n, 1))))
    return qmc.Sobol(d=3, scramble=True, seed=seed).random_base2(m)[:n]


def _uniform(n, seed):
    return np.random.default_rng(seed).random((n, 3))


DESIGNS = {"lhs": _lhs, "random-cd": _random_cd, "sobol": _sobol, "uniform": _uniform}


def load_mini() -> list[dict]:
    mini = REF / "cmip6_mini"
    grid = pd.read_csv(mini / "grid_coordinates.csv")
    frechet = pd.read_csv(mini / "frechet.csv").values
    r_est = pd.read_csv(mini / "local_estimates_all.csv")[["a_est", "b_est", "g_est"]].values
    scalars = pd.read_csv(mini / "lec_scalar_params.csv").set_index("parameter")["value"]
    cfg = CMIP6Config(
        neighbor_radius=float(scalars["locest_abst"]),
        smoothing_radius=float(scalars["smoothing_dist"]),
        df=5.0,
        alpha=1.0,
    )
    coords = np.column_stack([grid["Y"].values, grid["X"].values])

    captured: list[dict] = []

    def capture(zi, zj, xl, yl, df, alpha, **kw):
        captured.append({
            "name": f"mini:{len(captured)}",
            "args": (zi, zj, xl, yl, df, alpha),
            "kw": {"max_boundary_retries": kw["max_boundary_retries"]},
        })
        return np.zeros(3)

    with mock.patch.object(backend, "optimize_local_mle", capture):
        for c in range(len(grid)):
            _local_mle_one_cmip6(frechet, c, coords, cfg)
    assert len(captured) == len(grid), "every mini cell reaches the optimizer"
    for cell, r in zip(captured, r_est):
        cell["r_est"] = r
    return captured


def load_selected() -> list[dict]:
    out = pd.read_csv(REF / "local_optimizer_selected_outputs.csv")
    with open(REF / "local_optimizer_selected_inputs.bin", "rb") as f:
        n_rows = int(np.frombuffer(f.read(4), dtype=np.int32)[0])
        gp, zi, zj, xl, yl = (np.frombuffer(f.read(n_rows * 8), dtype=np.float64) for _ in range(5))
    cells = []
    for _, row in out.iterrows():
        m = gp == int(row["grid_point"])
        cells.append({
            "name": f"selected:{int(row['grid_point'])}",
            "args": (zi[m], zj[m], xl[m], yl[m], 5.0, 1.0),
            "kw": {
                "max_boundary_retries": 5,
                "lower_bounds": (float(row["lower_a"]), float(row["lower_b"])),
                "upper_bounds": (float(row["upper_a"]), float(row["upper_b"])),
            },
            "r_est": row[["a_est", "b_est", "g_est"]].to_numpy(float),
        })
    return cells


def run_config(cells, design, ensemble, seed):
    excess = np.empty(len(cells))
    t0 = time.perf_counter()
    with mock.patch.object(backend, "_local_mle_starts", DESIGNS[design]):
        for c, cell in enumerate(cells):
            p = backend.optimize_local_mle(*cell["args"], ensemble=ensemble, seed=seed, **cell["kw"])
            excess[c] = (
                neg_log_likelihood_sum(*cell["args"], *p)
                - neg_log_likelihood_sum(*cell["args"], *cell["r_est"])
            )
    return excess, time.perf_counter() - t0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ensembles", type=int, nargs="+", default=[5, 8, 10, 15, 20])
    ap.add_argument("--designs", nargs="+", choices=sorted(DESIGNS), default=["lhs", "random-cd"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[42, 0, 1])
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--csv", type=Path, help="write per-cell excess here")
    args = ap.parse_args()

    cells = load_mini() + load_selected()
    n_mini = sum(c["name"].startswith("mini:") for c in cells)
    rows = []
    summary = []
    for design in args.designs:
        for ens in args.ensembles:
            per_seed = [(seed, *run_config(cells, design, ens, seed)) for seed in args.seeds]
            rows += [
                (design, ens, seed, cells[c]["name"], e)
                for seed, excess, _ in per_seed
                for c, e in enumerate(excess)
            ]
            ex = np.array([e for _, e, _ in per_seed])
            secs = np.array([s for _, _, s in per_seed])
            worst_seed, worst_cell = np.unravel_index(ex.argmax(), ex.shape)
            s = {
                "design": design,
                "ensemble": ens,
                "worst_excess": ex.max(),
                "worst_at": f"seed {args.seeds[worst_seed]} {cells[worst_cell]['name']}",
                "worst_selected": ex[:, n_mini:].max(),
                "worse_cells_max_over_seeds": int((ex > args.tol).sum(axis=1).max()),
                "worse_cells_union": int((ex > args.tol).any(axis=0).sum()),
                "better_cells_min_over_seeds": int((ex < -args.tol).sum(axis=1).min()),
                "seeds_failing": int((ex > args.tol).any(axis=1).sum()),
                "sec_per_seed_median": float(np.median(secs)),
                "pass": bool(ex.max() <= args.tol),
            }
            summary.append(s)
            print(
                f"{design:9s} ens={ens:2d}  worst={s['worst_excess']:+.4f} ({s['worst_at']})  "
                f"selected worst={s['worst_selected']:+.4f}  "
                f"worse>tol max/seed={s['worse_cells_max_over_seeds']} union={s['worse_cells_union']}  "
                f"better min/seed={s['better_cells_min_over_seeds']}  "
                f"seeds failing={s['seeds_failing']}/{len(args.seeds)}  "
                f"{s['sec_per_seed_median']:.2f}s/{len(cells)} cells  "
                f"{'PASS' if s['pass'] else 'FAIL'}",
                flush=True,
            )

    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["design", "ensemble", "seed", "cell", "nll_excess"])
            w.writerows(rows)

    passing = [s for s in summary if s["pass"]]
    if passing:
        best = min(passing, key=lambda s: s["sec_per_seed_median"])
        print(f"cheapest passing: {best['design']} ensemble={best['ensemble']}")
    else:
        print("no configuration passes")


if __name__ == "__main__":
    main()

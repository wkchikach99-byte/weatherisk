"""End-to-end pipeline orchestration (replaces run.sh / SLURM workflow).

Provides run_pipeline() for stationary validation and
run_nonstationary_pipeline() for reproducing the full Contzen et al.
(2025) simulation study (Figures 3, 7, 8).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np

from weatherisk.grid import Grid
from weatherisk.covariance import cov_fkt_2d
from weatherisk.simulation import sim_expt_2d, sim_expt_2d_nonstat
from weatherisk.clustering import (
    calc_distance_ellipses,
    c_extrcoeff_matrix,
    clustering,
    cluster_number_threshold_method,
    quantile_threshold,
)
from weatherisk.density import pairwise_density_optim_local, run_local_mle_parallel
from weatherisk.estimation import (
    smooth_local_estimates,
    calc_estimates_in_clusters,
)
from weatherisk.parameters import ParameterPreset, get_preset
from scipy.cluster.hierarchy import fcluster


def run_pipeline(
    resolution: int = 10,
    n_sim: int = 10,
    df: float = 5.0,
    alpha: float = 1.0,
    seed: int = 42,
    output_dir: str | None = None,
    a: float = 2.0,
    b: float = 1.0,
    g: float = 0.0,
) -> dict[str, Any]:
    """Run the full validation pipeline on a synthetic grid.

    Parameters
    ----------
    resolution : int
        Grid resolution (points per axis).
    n_sim : int
        Number of simulations.
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent.
    seed : int
        Random seed.
    output_dir : str, optional
        Directory to save intermediate results.
    a, b, g : float
        Ellipse parameters for stationary simulation.

    Returns
    -------
    dict
        Contains 'clusters' (1-D array of length resolution^2),
        'sim_data', 'grid', 'linkage', and other results.
    """
    rng = np.random.default_rng(seed)
    grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=resolution)

    # Step 1: Simulate
    sim_data = sim_expt_2d(grid, df, alpha, a, b, g, n_sim=n_sim, rng=rng)

    # Step 2: Create simple local estimates from known parameters
    # For a small validation run, use ground-truth with small noise
    n_grid = grid.n_grid
    estimates = np.column_stack([
        np.full(n_grid, a) + rng.normal(0, 0.1, n_grid),
        np.full(n_grid, b) + rng.normal(0, 0.1, n_grid),
        np.full(n_grid, g) + rng.normal(0, 0.05, n_grid),
    ])
    estimates[:, 0] = np.maximum(estimates[:, 0], 0.01)
    estimates[:, 1] = np.maximum(estimates[:, 1], 0.0)

    # Step 3: Compute dissimilarity and cluster
    dm = calc_distance_ellipses(estimates, res=21)
    hc = clustering(dm)

    # Choose k based on grid size
    k = max(2, min(resolution, 10))
    labels = fcluster(hc, t=k, criterion="maxclust")
    # Convert 1-based to 0-based
    clusters = labels - 1

    result = {
        "clusters": clusters,
        "sim_data": sim_data,
        "grid": grid,
        "estimates": estimates,
        "linkage": hc,
        "dissimilarity": dm,
    }

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        np.save(os.path.join(output_dir, "clusters.npy"), clusters)
        np.save(os.path.join(output_dir, "estimates.npy"), estimates)

    return result


# ── Non-stationary pipeline (Contzen et al. 2025) ──────────────────


def run_nonstationary_pipeline(
    preset: str | ParameterPreset = "paper_stripes",
    resolution: int | None = None,
    n_sim: int | None = None,
    seed: int = 42,
    threshold_quantile: float = 0.30,
    output_dir: str | None = None,
    verbose: bool = True,
    n_workers: int = 1,
) -> dict[str, Any]:
    """Run the full non-stationary simulation + clustering pipeline.

    Reproduces the analysis from Contzen et al. (2025, Extremes 28:713–737):
    simulate non-stationary Huser-Genton process → local MLE at each grid
    cell → smooth → LEC & EDC clustering → in-cluster re-estimation.

    Parameters
    ----------
    preset : str or ParameterPreset
        Preset name (e.g. 'paper_stripes') or a ParameterPreset instance.
    resolution : int, optional
        Override preset resolution.
    n_sim : int, optional
        Override preset n_sim.
    seed : int
        Random seed.
    threshold_quantile : float
        Quantile (0–1) of the dissimilarity matrix upper triangle used
        as the clustering cut-off threshold (paper uses 0.30).
    output_dir : str, optional
        Save intermediate numpy arrays to this directory.
    verbose : bool
        Print progress messages.
    n_workers : int
        Number of parallel workers for local MLE (default 1 = serial).

    Returns
    -------
    dict
        Keys: grid, preset, a_matrix, b_matrix, g_matrix, sim_data,
        local_estimates, smoothed, dm_lec, hc_lec, labels_lec,
        k_lec, inclusters_lec, dm_edc, hc_edc, labels_edc,
        k_edc, inclusters_edc.
    """
    if isinstance(preset, str):
        p = get_preset(preset)
    else:
        p = preset

    res = resolution if resolution is not None else p.resolution
    ns = n_sim if n_sim is not None else p.n_sim

    import time as _time

    def _log(msg: str) -> None:
        if verbose:
            print(msg, flush=True)

    t_total = _time.perf_counter()

    _log(f"\n{'='*60}")
    _log(f"  Preset: {p.name}  (resolution={res}, n_sim={ns})")
    _log(f"{'='*60}")

    rng = np.random.default_rng(seed)
    grid = Grid(x_range=(-5, 5), y_range=(-5, 5), resolution=res)

    # Build true parameter fields
    a_matrix = p.a_func(grid.X, grid.Y)
    b_matrix = p.b_func(grid.X, grid.Y)
    g_matrix = p.g_func(grid.X, grid.Y)

    # ── Step 1: Simulate non-stationary max-stable process ──
    _log("  [1/6] Simulating non-stationary process...")
    t0 = _time.perf_counter()
    sim_data = sim_expt_2d_nonstat(
        grid, p.df, p.alpha, a_matrix, b_matrix, g_matrix,
        n_sim=ns, rng=rng,
    )
    _log(f"        {sim_data.shape}  ({_time.perf_counter() - t0:.1f}s)")

    # ── Step 2: Local MLE at each grid cell ──
    _log("  [2/6] Local MLE at each grid cell...")
    t0 = _time.perf_counter()
    local_estimates = run_local_mle_parallel(
        sim_data, grid, p.df, p.alpha,
        abstand=p.locest_abst,
        ensemble=p.locest_ensemble,
        n_workers=n_workers,
        verbose=verbose,
    )
    _log(f"        ({_time.perf_counter() - t0:.1f}s)")

    # ── Step 3: Smooth local estimates ──
    _log("  [3/6] Smoothing local estimates...")
    t0 = _time.perf_counter()
    smoothed = smooth_local_estimates(local_estimates, p.smoothing_dist, grid)
    _log(f"        ({_time.perf_counter() - t0:.1f}s)")

    # ── Step 4: LEC clustering (ellipse dissimilarity D2) ──
    _log("  [4/6] LEC clustering (ellipse dissimilarity)...")
    t0 = _time.perf_counter()
    dm_lec = calc_distance_ellipses(smoothed, res=21)
    hc_lec = clustering(dm_lec)
    thr_lec = quantile_threshold(dm_lec, threshold_quantile)
    k_lec = cluster_number_threshold_method(hc_lec, thr_lec)
    k_lec = max(2, k_lec)
    labels_lec = fcluster(hc_lec, t=k_lec, criterion="maxclust")
    _log(f"        k={k_lec}  ({_time.perf_counter() - t0:.1f}s)")

    # ── Step 5: EDC clustering (Saunders / madogram D1) ──
    _log("  [5/6] EDC clustering (madogram-based)...")
    t0 = _time.perf_counter()
    dm_edc = c_extrcoeff_matrix(sim_data, madogram=True)
    hc_edc = clustering(dm_edc)
    thr_edc = quantile_threshold(dm_edc, threshold_quantile)
    k_edc = cluster_number_threshold_method(hc_edc, thr_edc)
    k_edc = max(2, k_edc)
    labels_edc = fcluster(hc_edc, t=k_edc, criterion="maxclust")
    _log(f"        k={k_edc}  ({_time.perf_counter() - t0:.1f}s)")

    # ── Step 6: In-cluster re-estimation ──
    _log("  [6/6] In-cluster re-estimation...")
    t0 = _time.perf_counter()
    inclusters_lec = calc_estimates_in_clusters(
        sim_data, labels_lec, p.df, p.alpha, grid,
    )
    inclusters_edc = calc_estimates_in_clusters(
        sim_data, labels_edc, p.df, p.alpha, grid,
    )
    _log(f"        ({_time.perf_counter() - t0:.1f}s)")

    _log(f"  Total: {_time.perf_counter() - t_total:.1f}s")

    result = {
        "grid": grid,
        "preset": p,
        "a_matrix": a_matrix,
        "b_matrix": b_matrix,
        "g_matrix": g_matrix,
        "sim_data": sim_data,
        "local_estimates": local_estimates,
        "smoothed": smoothed,
        "dm_lec": dm_lec,
        "hc_lec": hc_lec,
        "k_lec": k_lec,
        "labels_lec": labels_lec,
        "inclusters_lec": inclusters_lec,
        "dm_edc": dm_edc,
        "hc_edc": hc_edc,
        "k_edc": k_edc,
        "labels_edc": labels_edc,
        "inclusters_edc": inclusters_edc,
    }

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        for key in ("labels_lec", "labels_edc", "local_estimates",
                     "smoothed", "inclusters_lec", "inclusters_edc"):
            np.save(os.path.join(output_dir, f"{key}.npy"), result[key])
        _log(f"  Results saved to {output_dir}")

    return result

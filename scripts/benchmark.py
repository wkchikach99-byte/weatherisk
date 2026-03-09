#!/usr/bin/env python3
"""Benchmark the vectorized pipeline components."""
import time
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from weatherisk.grid import Grid
from weatherisk.simulation import sim_expt_2d_nonstat
from weatherisk.clustering import c_extrcoeff_matrix, calc_distance_ellipses
from weatherisk.density import run_local_mle_parallel

RES = 21
N_SIM = 10
MLE_RES = 5  # smaller grid for MLE timing

print(f"=== Benchmark at {RES}x{RES}, {N_SIM} sims ===\n")

grid = Grid(resolution=RES)
a_mat = np.full_like(grid.X, 2.0)
b_mat = (grid.X + 5) / 2
g_mat = np.zeros_like(grid.X)

# 1. Simulation (cov build + Poisson spectral)
t0 = time.time()
sim = sim_expt_2d_nonstat(grid, 5.0, 1.0, a_mat, b_mat, g_mat,
                          n_sim=N_SIM, rng=np.random.default_rng(42))
t1 = time.time()
print(f"1. Simulation {RES}x{RES}x{N_SIM}: {t1-t0:.2f}s")

# 2. EDC matrix (madogram)
t2 = time.time()
edc = c_extrcoeff_matrix(sim, madogram=True)
t3 = time.time()
print(f"2. EDC matrix ({grid.n_grid} cells): {t3-t2:.3f}s")

# 3. LEC matrix (ellipse dissimilarity)
rng = np.random.default_rng(1)
estimates = rng.random((grid.n_grid, 3)) * [3, 2, np.pi] + [0.5, 0.1, -np.pi/2]
t4 = time.time()
lec = calc_distance_ellipses(estimates, res=11)
t5 = time.time()
print(f"3. LEC matrix ({grid.n_grid} cells, res=11): {t5-t4:.3f}s")

# 4. Local MLE — serial (smaller grid to finish quickly)
grid_sm = Grid(resolution=MLE_RES)
a_sm = np.full_like(grid_sm.X, 2.0)
b_sm = (grid_sm.X + 5) / 2
g_sm = np.zeros_like(grid_sm.X)
sim_sm = sim_expt_2d_nonstat(grid_sm, 5.0, 1.0, a_sm, b_sm, g_sm,
                             n_sim=N_SIM, rng=np.random.default_rng(42))

t6 = time.time()
est_serial = run_local_mle_parallel(sim_sm, grid_sm, 5.0, 1.0,
                                     abstand=3, ensemble=1,
                                     n_workers=1, verbose=False)
t7 = time.time()
print(f"4. Local MLE serial ({grid_sm.n_grid} cells, ens=1): {t7-t6:.2f}s")

# 5. Local MLE — parallel (8 workers)
t8 = time.time()
est_par = run_local_mle_parallel(sim_sm, grid_sm, 5.0, 1.0,
                                  abstand=3, ensemble=1,
                                  n_workers=8, verbose=False)
t9 = time.time()
print(f"5. Local MLE 8 workers ({grid_sm.n_grid} cells, ens=1): {t9-t8:.2f}s")

# Verify parallel == serial
diff = np.max(np.abs(est_serial - est_par))
print(f"\n   Max diff serial vs parallel: {diff:.2e}")

total_serial = (t1-t0) + (t3-t2) + (t5-t4) + (t7-t6)
total_par = (t1-t0) + (t3-t2) + (t5-t4) + (t9-t8)
print(f"\n   Total (serial MLE): {total_serial:.1f}s")
print(f"   Total (8-worker MLE): {total_par:.1f}s")
print(f"   Speedup on MLE step: {(t7-t6)/(t9-t8):.1f}x")

"""Inspect step2 checkpoint for original-scale precipitation data."""
import numpy as np

d = np.load("output/cmip6_fig9_local/checkpoints/step2.npz")
am = d["annual_max"]  # (156, 96, 192) in kg/m2/s
print(f"annual_max shape: {am.shape}")
print(f"annual_max range: [{am.min():.6f}, {am.max():.6f}] (kg/m2/s)")

am_mm = am * 86400  # convert to mm/day
print(f"In mm/day: [{am_mm.min():.2f}, {am_mm.max():.2f}]")
print(f"Mean annual max precip: {am_mm.mean():.2f} mm/day")
print(f"Median annual max: {np.median(am_mm):.2f} mm/day")
print(f"Size: {am.nbytes / 1e6:.1f} MB")

# Also check: can we quickly re-fit GEV for a single cell?
from scipy.stats import genextreme
cell_data = am_mm[:, 48, 96]  # equatorial cell
shape, loc, scale = genextreme.fit(cell_data)
print(f"\nExample GEV fit (equator cell 48,96):")
print(f"  xi={-shape:.4f}, mu={loc:.2f} mm/day, sigma={scale:.2f} mm/day")

# 100-year return level in mm/day
T = 100
rl100 = genextreme.isf(1/T, shape, loc=loc, scale=scale)
print(f"  100-year return level: {rl100:.1f} mm/day")

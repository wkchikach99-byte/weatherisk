"""Find ALL positions in R's dist_matrix that are corrupted by the linear-indexing bug."""
import numpy as np
import pandas as pd
from pathlib import Path

REF = Path("tests/reference_data/cmip6_mini")
r_smoothed = pd.read_csv(REF / "local_estimates_smoothed.csv")
smoothed = np.column_stack([
    r_smoothed["a_sm"].values,
    r_smoothed["b_sm"].values,
    r_smoothed["g_sm"].values,
])

from weatherisk.covariance import cov_fkt_2d

res = 21
xs = np.repeat(np.linspace(-1, 1, res), res)
ys = np.tile(np.linspace(-1, 1, res), res)
mask = ((xs**2 + ys**2) <= res**2) & ((ys > 0) | (xs > 0))
xs_m, ys_m = xs[mask], ys[mask]

n = len(smoothed)

# Simulate R's exact algorithm to find which matrix positions are written to
# by the linear-indexing bug
dist_r = np.zeros((n, n))
dist_correct = np.zeros((n, n))
r_bug_targets = set()

for i in range(n-1):
    for j in range(i+1, n):
        mx = max(smoothed[i,0]+smoothed[i,1], smoothed[j,0]+smoothed[j,1])
        cov_i = cov_fkt_2d(xs_m, ys_m, 1.0, smoothed[i,0]/mx, smoothed[i,1]/mx, smoothed[i,2])
        cov_j = cov_fkt_2d(xs_m, ys_m, 1.0, smoothed[j,0]/mx, smoothed[j,1]/mx, smoothed[j,2])
        ell1 = (cov_i > np.exp(-1)).astype(float)
        ell2 = (cov_j > np.exp(-1)).astype(float)

        if max(ell1.sum(), ell2.sum()) == 0:
            # R bug: dist_matrix[j] = 1  (j is 1-based in R)
            j_r = j + 1  # convert to R's 1-based j
            lin_idx = j_r - 1  # back to 0-based linear index
            row = lin_idx % n
            col = lin_idx // n
            dist_r[row, col] = 1.0  # overwrites whatever was there
            r_bug_targets.add((row, col))
            # Correct value should be at (i, j)
            dist_correct[i, j] = 1.0
        else:
            inter = (ell1 * ell2).sum() + 0.5
            union = (ell1 + ell2 - ell1 * ell2).sum() + 0.5
            d = 1.0 - inter / union
            dist_r[i, j] = d
            dist_correct[i, j] = d

# Apply symmetrization
final_r = 100.0 * (dist_r + dist_r.T)
final_correct = 100.0 * (dist_correct + dist_correct.T)

# Verify against reference
r_ref = pd.read_csv(REF / "ellipse_dissimilarity_matrix.csv").values
print(f"R simulation matches reference: {np.allclose(final_r, r_ref)}")

# Find all positions where R's bug affects the result
affected = np.abs(final_r - final_correct) > 0.01
print(f"\nTotal affected positions: {affected.sum()}")
print(f"Bug target (row,col) positions (0-based):")
for pos in sorted(r_bug_targets):
    print(f"  {pos} — R writes 1.0 here via linear index")

print(f"\nAll affected (i,j) positions in final matrix:")
aff_i, aff_j = np.where(affected)
for k in range(len(aff_i)):
    print(f"  ({aff_i[k]:2d},{aff_j[k]:2d}): R={final_r[aff_i[k],aff_j[k]]:.2f}  correct={final_correct[aff_i[k],aff_j[k]]:.2f}")

"""Size of the tied annual-maxima groups without detrending (run from the repo root).

Usage: python hpc/r_reference/tie_size_check.py DATA_DIR
"""

import sys

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, "hpc/r_reference")
from edc_preprocessing_sweep import load_pr  # noqa: E402

flat, _ = load_pr(sys.argv[1], 1850, 2014)
mx = flat[: 12 * 156].reshape(156, 12, -1).max(axis=1)
tied = [c for c in range(mx.shape[1]) if np.unique(mx[:, c]).size < 156]
sizes = [156 - np.unique(mx[:, c]).size for c in tied]
print("tied cells", len(tied), "years lost to ties per cell: min", min(sizes),
      "median", int(np.median(sizes)), "max", max(sizes))
print("cells whose minimum maximum is zero", sum(np.min(mx[:, c]) == 0 for c in tied))
r = {m: rankdata(mx[:, tied], method=m, axis=0) for m in ("min", "average", "max")}
print("max |min-average| rank", np.abs(r["min"] - r["average"]).max(),
      " max |min-max|", np.abs(r["min"] - r["max"]).max())

"""Do average and min-equivalent ranks give the same EDC labels on a Fig. 9 run?

Recomputes step 5 EDC from the step-2 checkpoint and compares with step5.npz.
Usage: python hpc/r_reference/edc_labels_rule_check.py CHECKPOINT_DIR
"""
import sys
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import rankdata

ck = sys.argv[1]
fr = np.load(f"{ck}/step2.npz")["frechet"]
s5 = np.load(f"{ck}/step5.npz")
n = fr.shape[0]
print("frechet", fr.shape, "cells with ties", int(sum(np.unique(fr[:, c]).size < n for c in range(fr.shape[1]))))
out = {}
for name, ranks in (("average", rankdata(fr, axis=0)), ("min-equivalent", rankdata(-fr, method="max", axis=0))):
    cond = pdist(ranks.T, metric="cityblock") / (n * 2.0 * (n + 1))
    thr = float(np.quantile(cond, 0.30))
    hc = linkage(cond, method="average")
    k = max(2, int(np.sum(hc[:, 2] > thr)))
    out[name] = fcluster(hc, t=k, criterion="maxclust")
    print(name, "k", k, "q30", f"{thr:.8f}", flush=True)
    del cond, hc
print("average labels == step5 checkpoint:", np.array_equal(out["average"], s5["labels_edc"]), "k_edc ckpt", int(s5["k_edc"]))
a, b = out["average"], out["min-equivalent"]
# same partition up to label names?
pairs = set(zip(a.tolist(), b.tolist()))
print("same partition (up to relabelling):", len(pairs) == len(set(a)) == len(set(b)))
print("identical labels:", np.array_equal(a, b), "cells with different label:", int(np.sum(a != b)))

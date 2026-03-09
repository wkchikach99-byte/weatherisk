"""Risk metrics: VaR and ES computation per cluster.

Value at Risk (VaR) and Expected Shortfall (ES) are standard tail-risk
measures applied to spatial extreme-value data grouped by cluster.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def compute_var(data: np.ndarray, p: float = 0.95) -> float:
    """Compute Value at Risk at level p (the p-quantile)."""
    return float(np.quantile(data, p))


def compute_es(data: np.ndarray, p: float = 0.95) -> float:
    """Compute Expected Shortfall at level p.

    ES_p = E[L | L > VaR_p].
    """
    var = compute_var(data, p)
    tail = data[data >= var]
    if len(tail) == 0:
        return var
    return float(np.mean(tail))


def compute_cluster_risk(
    data: np.ndarray,
    clusters: np.ndarray,
    p: float = 0.95,
) -> list[dict[str, Any]]:
    """Compute VaR and ES per cluster.

    Parameters
    ----------
    data : ndarray, shape (n_realisations, n_cells)
        Simulated or observed values.
    clusters : 1-D array of length n_cells
        Cluster label for each grid cell.
    p : float
        Probability level.

    Returns
    -------
    list[dict]
        One dict per cluster with 'cluster', 'var', 'es' keys.
    """
    result: list[dict[str, Any]] = []
    for label in sorted(np.unique(clusters)):
        mask = clusters == label
        cluster_max = data[:, mask].max(axis=1)
        var = compute_var(cluster_max, p)
        es = compute_es(cluster_max, p)
        result.append({"cluster": int(label), "var": var, "es": es})
    return result

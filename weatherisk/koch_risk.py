"""Koch (2016) spatial risk measures.

Computes the normalised spatially aggregated loss L_N, Value-at-Risk,
and Expected Shortfall per cluster following Koch (2016, Section 3).

Supports two cost functions E(z):
  - "indicator":  E(z) = 1{z > u_p}         (area fraction in exceedance)
  - "excess":     E(z) = max(z - u_p, 0)     (severity-weighted loss)

All operations are vectorised NumPy; the full computation for 163
clusters × 156 years runs in ~20 ms.

Usage
-----
    from weatherisk.koch_risk import compute_koch_risk
    results = compute_koch_risk("risk/inputs/pipeline_results.npz")
    results_excess = compute_koch_risk(
        "risk/inputs/pipeline_results.npz", cost="excess")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class KochRiskResult:
    """Container for Koch risk results for a single clustering."""

    clustering_name: str
    cost: str  # "indicator" or "excess"
    p: float
    q: float
    u_p: float
    years: np.ndarray
    clusters: list[dict[str, Any]] = field(default_factory=list)


def _cosine_weights(
    lats: np.ndarray,
    valid_idx: np.ndarray,
    n_lon: int,
) -> np.ndarray:
    """Area-proportional cosine-latitude weights for each valid cell."""
    lat_indices = (valid_idx // n_lon).astype(int)
    return np.cos(np.deg2rad(lats[lat_indices]))


def _frechet_threshold(p: float) -> float:
    """Fréchet quantile: u_p = 1 / (-log p)."""
    return 1.0 / (-np.log(p))


def compute_cluster_loss(
    frechet: np.ndarray,
    labels: np.ndarray,
    cos_weights: np.ndarray,
    u_p: float,
    cost: str = "indicator",
) -> dict[int, np.ndarray]:
    """Compute Koch's L_N per cluster per year.

    Parameters
    ----------
    frechet : (n_years, n_cells)
    labels : (n_cells,) cluster assignments
    cos_weights : (n_cells,) cosine-latitude weights
    u_p : Fréchet threshold
    cost : "indicator" for E(z) = 1{z > u_p},
           "excess"    for E(z) = max(z - u_p, 0)

    Returns
    -------
    dict mapping cluster_id -> L_N array of shape (n_years,)
    """
    if cost == "indicator":
        E = (frechet > u_p).astype(np.float64)
    elif cost == "excess":
        E = np.maximum(frechet - u_p, 0.0)
    else:
        raise ValueError(f"Unknown cost function: {cost!r}")

    cluster_loss = {}
    for cl in np.unique(labels):
        mask = labels == cl
        w = cos_weights[mask]
        w_normed = w / w.sum()
        cluster_loss[int(cl)] = E[:, mask] @ w_normed
    return cluster_loss


def compute_var(L_N: np.ndarray, q: float) -> float:
    """Empirical VaR_q: the ceil(n*q)-th order statistic."""
    n = len(L_N)
    k = int(np.ceil(n * q)) - 1
    return float(np.sort(L_N)[k])


def compute_es(L_N: np.ndarray, q: float) -> float:
    """Empirical ES_q: mean of values at or above VaR_q."""
    n = len(L_N)
    k = int(np.ceil(n * q)) - 1
    return float(np.sort(L_N)[k:].mean())


def compute_koch_risk(
    data_path: str | Path = "risk/inputs/pipeline_results.npz",
    p: float = 0.95,
    q: float = 0.95,
    clusterings: list[str] | None = None,
    cost: str = "indicator",
) -> dict[str, KochRiskResult]:
    """Compute Koch risk measures for all clusters.

    Parameters
    ----------
    data_path : path to pipeline_results.npz
    p : Fréchet threshold probability (defines u_p)
    q : confidence level for VaR/ES
    clusterings : which clusterings to use; default ["LEC", "EDC"]
    cost : "indicator" or "excess"

    Returns
    -------
    dict mapping clustering name -> KochRiskResult
    """
    if clusterings is None:
        clusterings = ["LEC", "EDC"]

    d = np.load(data_path)
    frechet = d["frechet"]
    lats = d["lats"]
    valid_idx = d["valid_idx"]
    years = d["years"]
    n_lon = int(d["lons"].shape[0])

    cos_weights = _cosine_weights(lats, valid_idx, n_lon)
    u_p = _frechet_threshold(p)

    label_map = {
        "LEC": d["labels_lec"],
        "EDC": d["labels_edc"],
    }

    results = {}
    for cname in clusterings:
        labels = label_map[cname]
        cluster_loss = compute_cluster_loss(
            frechet, labels, cos_weights, u_p, cost=cost,
        )

        result = KochRiskResult(
            clustering_name=cname,
            cost=cost,
            p=p,
            q=q,
            u_p=u_p,
            years=years,
        )
        for cl_id, L_N in sorted(cluster_loss.items()):
            n_cells = int((labels == cl_id).sum())
            result.clusters.append({
                "cluster": cl_id,
                "n_cells": n_cells,
                "L_N": L_N,
                "L_N_mean": float(L_N.mean()),
                "L_N_max": float(L_N.max()),
                "VaR": compute_var(L_N, q),
                "ES": compute_es(L_N, q),
            })
        results[cname] = result

    return results


# ─── Temporal trend analysis ─────────────────────────────────────────

def _mann_kendall(x: np.ndarray) -> tuple[float, float]:
    """Mann-Kendall trend test.

    Returns (S, p_value) where S is the test statistic and p_value
    is the two-sided significance level (normal approximation, valid
    for n >= 10).
    """
    n = len(x)
    s = 0
    for k in range(n - 1):
        s += int(np.sign(x[k + 1:] - x[k]).sum())
    # Variance under H0: no ties correction (continuous data)
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0
    # Two-sided p-value from standard normal
    from scipy.stats import norm
    p_value = 2.0 * norm.sf(abs(z))
    return float(s), float(p_value)


def compute_trend(
    years: np.ndarray,
    L_N: np.ndarray,
) -> dict[str, float]:
    """OLS slope + Mann-Kendall significance for an L_N time series.

    Returns
    -------
    dict with keys: slope, intercept, slope_per_decade, mk_S, mk_pvalue
    """
    n = len(years)
    t = years - years[0]  # centre at 0 for numerical stability
    slope = float(np.polyfit(t, L_N, 1)[0])
    intercept = float(np.polyfit(t, L_N, 1)[1])
    mk_s, mk_p = _mann_kendall(L_N)
    return {
        "slope": slope,
        "intercept": intercept,
        "slope_per_decade": slope * 10.0,
        "mk_S": mk_s,
        "mk_pvalue": mk_p,
    }


# ─── Cross-cluster joint exceedance analysis ────────────────────────

@dataclass
class CrossClusterResult:
    """Container for cross-cluster co-exceedance analysis."""

    clustering_name: str
    cost: str
    q: float
    cluster_ids: list[int]
    # (K, K) matrices — cluster_ids[i] is row/col i
    joint_prob: np.ndarray       # empirical P(both exceed)
    indep_prob: np.ndarray       # P_A × P_B  (independence baseline)
    coexceedance_ratio: np.ndarray  # joint / indep  (1.0 = independent)
    n_joint_years: np.ndarray    # count of years both exceed


def compute_cross_cluster_dependence(
    result: KochRiskResult,
) -> CrossClusterResult:
    """Compute pairwise co-exceedance statistics for all clusters.

    For each pair (A, B), defines "bad year" as L_N > VaR_q and computes:
      - joint exceedance probability  P(A bad ∩ B bad)
      - independence baseline         P(A bad) × P(B bad)
      - co-exceedance ratio           joint / indep

    A ratio >> 1 indicates systemic risk: bad years co-occur more than
    chance predicts.  A ratio ≈ 1 means the clusters' extreme years are
    independent.

    Parameters
    ----------
    result : KochRiskResult from compute_koch_risk for one clustering

    Returns
    -------
    CrossClusterResult with (K, K) matrices indexed by cluster order
    """
    clusters = result.clusters
    K = len(clusters)
    n_years = len(result.years)
    q = result.q

    # Boolean exceedance indicators: (K, n_years)
    exceed = np.zeros((K, n_years), dtype=bool)
    marginal_prob = np.zeros(K)
    cluster_ids = []

    for i, cl in enumerate(clusters):
        cluster_ids.append(cl["cluster"])
        exceed[i] = cl["L_N"] >= cl["VaR"]
        marginal_prob[i] = exceed[i].mean()

    # Pairwise joint exceedance
    joint_prob = np.zeros((K, K))
    indep_prob = np.zeros((K, K))
    n_joint = np.zeros((K, K), dtype=int)

    for i in range(K):
        for j in range(K):
            both = exceed[i] & exceed[j]
            n_joint[i, j] = int(both.sum())
            joint_prob[i, j] = both.mean()
            indep_prob[i, j] = marginal_prob[i] * marginal_prob[j]

    # Co-exceedance ratio (guard against division by zero)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(indep_prob > 0, joint_prob / indep_prob, 0.0)

    return CrossClusterResult(
        clustering_name=result.clustering_name,
        cost=result.cost,
        q=q,
        cluster_ids=cluster_ids,
        joint_prob=joint_prob,
        indep_prob=indep_prob,
        coexceedance_ratio=ratio,
        n_joint_years=n_joint,
    )

"""Risk-map loading, quantile banding, connected-component clustering, and region statistics.

Port of risk_pipeline_moreclusters_graph_fix.py with bug fixes from PRD section 11.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, label
from scipy.spatial import cKDTree


def load_and_grid(csv_path: str) -> dict[str, Any]:
    """Load a risk-map CSV and reshape to 2-D grids.

    Returns a dict with keys: df, lats, lons, ES, VaR, EXP,
    land_mask, lon_grid, lat_grid.
    """
    df = pd.read_csv(csv_path)
    req = {"lat", "lon", "VaR_95", "ES_95"}
    if not req.issubset(df.columns):
        raise ValueError(f"CSV must have {req}. Found: {df.columns.tolist()}")

    if "exposure" not in df.columns:
        df["exposure"] = 1.0

    df = df.sort_values(["lat", "lon"]).reset_index(drop=True)
    lats = np.sort(df["lat"].unique())
    lons = np.sort(df["lon"].unique())
    n_rows, n_cols = len(lats), len(lons)

    def to_grid(col):
        return df[col].to_numpy().reshape((n_rows, n_cols))

    ES = to_grid("ES_95")
    VaR = to_grid("VaR_95")
    EXP = to_grid("exposure")

    land_mask = np.where(np.nan_to_num(EXP, nan=0.0) > 0, 1.0, np.nan)
    ES = ES * land_mask
    VaR = VaR * land_mask

    lon_grid, lat_grid = np.meshgrid(lons, lats)
    return {
        "df": df,
        "lats": lats,
        "lons": lons,
        "ES": ES,
        "VaR": VaR,
        "EXP": EXP,
        "land_mask": land_mask,
        "lon_grid": lon_grid,
        "lat_grid": lat_grid,
    }


def smooth_field(
    a: np.ndarray,
    sigma: float,
    land_mask: np.ndarray,
) -> np.ndarray:
    """Gaussian smooth a 2-D field, respecting land mask."""
    if sigma and sigma > 0:
        # Use nan-aware smoothing to avoid leaking zeros into valid cells
        valid = np.isfinite(a) & np.isfinite(land_mask)
        a_filled = np.where(valid, a, 0.0)
        weights = valid.astype(float)
        a_s = gaussian_filter(a_filled, sigma=sigma)
        w_s = gaussian_filter(weights, sigma=sigma)
        w_s = np.where(w_s > 0, w_s, 1.0)
        a_s = a_s / w_s
        a_s[~valid] = np.nan
        return a_s
    return a.copy()


def quantile_bands(
    a: np.ndarray,
    q: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Assign grid cells to quantile-based bands.

    Parameters
    ----------
    a : 2-D array
        Field values (NaN = missing).
    q : int
        Number of bands.

    Returns
    -------
    bands : int array same shape as a
        Band index (0..q-1) or -1 for excluded cells.
    edges : 1-D array
        Quantile bin edges.
    """
    valid = np.isfinite(a) & (a > 0)
    A = a[valid]

    if A.size == 0:
        bands = np.full(a.shape, -1, dtype=int)
        edges = np.array([-np.inf, np.inf], dtype=float)
        return bands, edges

    q_eff = min(q, max(2, np.unique(A).size))
    edges = np.quantile(A, np.linspace(0, 1, q_eff + 1))
    edges[0] = -np.inf
    edges[-1] = np.inf

    bands = np.digitize(a, edges[1:-1], right=False)
    bands[~valid] = -1
    return bands, edges


def connected_patches(
    profile_code: np.ndarray,
    min_cells: int = 1,
) -> np.ndarray:
    """Label contiguous regions sharing the same band code.

    Parameters
    ----------
    profile_code : 2-D int array
        Band assignment per cell (-1 = excluded).
    min_cells : int
        Minimum patch size.  Smaller patches get label -2.

    Returns
    -------
    cluster_id : 2-D int array
        Cluster labels (>=0) or -1 (excluded) or -2 (tiny).
    """
    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    cluster_id = -np.ones_like(profile_code, dtype=np.int32)
    next_id = 0

    for p in np.unique(profile_code):
        if p < 0:
            continue
        mask = (profile_code == p).astype(np.int8)
        if mask.sum() == 0:
            continue
        L, num = label(mask, structure=structure)
        for k in range(1, num + 1):
            sel = L == k
            n = int(sel.sum())
            if n < min_cells:
                cluster_id[sel] = -2
            else:
                cluster_id[sel] = next_id
                next_id += 1

    return cluster_id


def merge_tiny_regions(
    cluster_id: np.ndarray,
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
) -> np.ndarray:
    """Merge patches labelled -2 into their nearest large region."""
    tiny_mask = cluster_id == -2
    if not tiny_mask.any():
        return cluster_id.copy()

    big_ids = np.unique(cluster_id[cluster_id >= 0])
    out = cluster_id.copy()

    if big_ids.size == 0:
        structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        L, num = label(tiny_mask.astype(np.int8), structure=structure)
        cid = 0
        for k in range(1, num + 1):
            out[L == k] = cid
            cid += 1
        return out

    big_centroids = []
    for cid_val in big_ids:
        sel = out == cid_val
        cx = lon_grid[sel].mean()
        cy = lat_grid[sel].mean()
        big_centroids.append([cx, cy])
    tree = cKDTree(np.array(big_centroids))

    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    L, num = label(tiny_mask.astype(np.int8), structure=structure)
    for k in range(1, num + 1):
        sel = L == k
        cx = lon_grid[sel].mean()
        cy = lat_grid[sel].mean()
        _, idx = tree.query([[cx, cy]], k=1)
        nearest_cid = big_ids[int(idx[0])]
        out[sel] = nearest_cid

    return out


def remap_ids_to_sequential(cluster_id: np.ndarray) -> tuple[np.ndarray, int]:
    """Remap cluster IDs to 0..K-1."""
    out = np.full(cluster_id.shape, -1, dtype=np.int32)
    valid = cluster_id >= 0
    uniq = np.unique(cluster_id[valid])
    for i, cid in enumerate(uniq):
        out[cluster_id == cid] = i
    return out, int(len(uniq))


def compute_cluster_stats(
    out_df: pd.DataFrame,
    grid_df: pd.DataFrame,
) -> pd.DataFrame:
    """Per-cluster ES/VaR stats including exposure-weighted means."""
    if "exposure" not in grid_df.columns:
        grid_df = grid_df.copy()
        grid_df["exposure"] = 1.0

    g = out_df.merge(
        grid_df[["lat", "lon", "ES_95", "VaR_95", "exposure"]],
        on=["lat", "lon"],
        how="left",
    )
    g = g[g["cluster"] >= 0].copy()

    gb = g.groupby("cluster", as_index=False, sort=True)

    basic = gb.agg(
        n_cells=("cluster", "size"),
        exposure_sum=("exposure", "sum"),
        ES_mean=("ES_95", "mean"),
        ES_median=("ES_95", "median"),
        ES_sum=("ES_95", "sum"),
        VaR_mean=("VaR_95", "mean"),
        VaR_median=("VaR_95", "median"),
        VaR_sum=("VaR_95", "sum"),
    )

    g["_wES"] = np.nan_to_num(g["ES_95"].values, nan=0.0) * np.clip(
        np.nan_to_num(g["exposure"].values, nan=0.0), 0, None
    )
    g["_wVaR"] = np.nan_to_num(g["VaR_95"].values, nan=0.0) * np.clip(
        np.nan_to_num(g["exposure"].values, nan=0.0), 0, None
    )

    wsum = gb["exposure"].sum().rename(columns={"exposure": "_wsum"})
    wesum = gb["_wES"].sum().rename(columns={"_wES": "_wesum"})
    wvsum = gb["_wVaR"].sum().rename(columns={"_wVaR": "_wvsum"})

    tmp = wsum.merge(wesum, on="cluster").merge(wvsum, on="cluster")
    tmp["ES_popw_mean"] = np.where(
        tmp["_wsum"] > 0, tmp["_wesum"] / tmp["_wsum"], np.nan
    )
    tmp["VaR_popw_mean"] = np.where(
        tmp["_wsum"] > 0, tmp["_wvsum"] / tmp["_wsum"], np.nan
    )
    tmp = tmp[["cluster", "ES_popw_mean", "VaR_popw_mean"]]

    stats = basic.merge(tmp, on="cluster", how="left")
    return stats

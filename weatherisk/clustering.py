"""Dissimilarity matrices, hierarchical clustering, and k-selection.

Implements the ellipse-shape Jaccard-like dissimilarity from
Justus's paper, the madogram-based Saunders method, and
hierarchical agglomerative clustering with threshold-based k.
"""

from __future__ import annotations

import numpy as np

from weatherisk.covariance import cov_fkt_2d


def calc_distance_ellipses(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Compute Jaccard-like ellipse overlap dissimilarity matrix.

    For each pair of local estimates, rasterise normalised ellipses on a
    half-circle grid and compute 1 - IoU (intersection over union).

    This version pre-computes the unscaled quadratic-form values for every
    grid point, so the inner pair loop only does threshold comparisons and
    boolean set operations — ~20×50× faster than the scalar version.

    Parameters
    ----------
    estimates : ndarray, shape (n, 3)
        Local estimates (a, b, g) per grid point.
    res : int
        Resolution of the rasterisation grid.
    chunk_size : int, optional
        Row-chunk size for the pairwise comparison loop. Smaller chunks use
        less peak memory at the cost of additional loop overhead.

    Returns
    -------
    ndarray, shape (n, n)
        Symmetric dissimilarity matrix (0 on diagonal), scaled 0-100.
    """
    n = estimates.shape[0]

    # Build evaluation points (half-circle for symmetry)
    xs = np.repeat(np.linspace(-1, 1, res), res)
    ys = np.tile(np.linspace(-1, 1, res), res)
    mask = ((xs ** 2 + ys ** 2) <= res ** 2) & ((ys > 0) | (xs > 0))
    xs = xs[mask]
    ys = ys[mask]
    n_pts = len(xs)

    # Pre-compute sqrt(Q_i(xs, ys)) for each grid point i
    # cov_fkt_2d(xs, ys, alpha=1, a, b, g) = exp(-sqrt(Q))
    # When we scale (a, b) by 1/mx, Q is multiplied by mx², so
    # sqrt(Q_scaled) = mx * sqrt(Q_original).
    # Threshold: exp(-mx * sq) > exp(-1) ⟺ sq < 1/mx
    sq = np.empty((n, n_pts))
    for i in range(n):
        a_i, b_i, g_i = estimates[i]
        sg, cg = np.sin(g_i), np.cos(g_i)
        ap = a_i + b_i
        if a_i == 0 and ap == 0:
            sq[i] = np.inf
            continue
        qf = (
            xs * xs * (sg * sg / (a_i * a_i) + cg * cg / (ap * ap))
            + 2 * xs * ys * sg * cg * (-1.0 / (a_i * a_i) + 1.0 / (ap * ap))
            + ys * ys * (cg * cg / (a_i * a_i) + sg * sg / (ap * ap))
        )
        sq[i] = np.sqrt(np.maximum(qf, 0.0))

    ab = estimates[:, 0] + estimates[:, 1]  # semi-major for each point

    # Compute dissimilarity matrix in chunks to limit memory
    dist_matrix = np.zeros((n, n))
    chunk = min(chunk_size or 256, n)

    for i_start in range(0, n, chunk):
        i_end = min(i_start + chunk, n)
        ch = i_end - i_start

        # mx[ci, j] = max(ab[i_start+ci], ab[j])
        mx = np.maximum(ab[i_start:i_end, None], ab[None, :])  # (ch, n)
        mx = np.maximum(mx, 1e-300)  # avoid div by zero
        thr = 1.0 / mx  # (ch, n)

        # mask_i[ci, j, p] = sq[i_start+ci, p] < thr[ci, j]
        sq_chunk = sq[i_start:i_end]  # (ch, n_pts)
        mask_i = sq_chunk[:, None, :] < thr[:, :, None]  # (ch, n, n_pts)
        mask_j = sq[None, :, :] < thr[:, :, None]  # (ch, n, n_pts)

        inter = (mask_i & mask_j).sum(axis=2).astype(np.float64) + 0.5
        union = (mask_i | mask_j).sum(axis=2).astype(np.float64) + 0.5

        d = np.where(union == 0.5, 1.0, 1.0 - inter / union)
        dist_matrix[i_start:i_end, :] = d

    np.fill_diagonal(dist_matrix, 0.0)
    return 100.0 * dist_matrix


def calc_distance_ellipses_condensed(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Condensed upper-triangle ellipse overlap dissimilarity vector.

    Same algorithm as :func:`calc_distance_ellipses` but outputs the
    condensed form directly (length ``n*(n-1)/2``), avoiding the
    full ``(n, n)`` matrix allocation.  This reduces peak memory
    from ``O(n²)`` float64 to ``O(n²/2)`` float64.

    Parameters
    ----------
    estimates : ndarray, shape (n, 3)
        Local estimates (a, b, g) per grid point.
    res : int
        Resolution of the rasterisation grid.
    chunk_size : int, optional
        Row-chunk size for the pairwise comparison loop.

    Returns
    -------
    ndarray, shape (n*(n-1)//2,)
        Condensed dissimilarity vector (scipy squareform order), 0–100.
    """
    n = estimates.shape[0]
    n_pairs = n * (n - 1) // 2

    # Build evaluation points (half-circle for symmetry)
    xs = np.repeat(np.linspace(-1, 1, res), res)
    ys = np.tile(np.linspace(-1, 1, res), res)
    mask = ((xs ** 2 + ys ** 2) <= res ** 2) & ((ys > 0) | (xs > 0))
    xs = xs[mask]
    ys = ys[mask]
    n_pts = len(xs)

    # Pre-compute sqrt(Q_i) for each grid point
    sq = np.empty((n, n_pts))
    for i in range(n):
        a_i, b_i, g_i = estimates[i]
        sg, cg = np.sin(g_i), np.cos(g_i)
        ap = a_i + b_i
        if a_i == 0 and ap == 0:
            sq[i] = np.inf
            continue
        qf = (
            xs * xs * (sg * sg / (a_i * a_i) + cg * cg / (ap * ap))
            + 2 * xs * ys * sg * cg * (-1.0 / (a_i * a_i) + 1.0 / (ap * ap))
            + ys * ys * (cg * cg / (a_i * a_i) + sg * sg / (ap * ap))
        )
        sq[i] = np.sqrt(np.maximum(qf, 0.0))

    ab = estimates[:, 0] + estimates[:, 1]

    condensed = np.empty(n_pairs)
    chunk = min(chunk_size or 256, n)

    for i_start in range(0, n, chunk):
        i_end = min(i_start + chunk, n)
        ch = i_end - i_start

        mx = np.maximum(ab[i_start:i_end, None], ab[None, :])
        mx = np.maximum(mx, 1e-300)
        thr = 1.0 / mx

        sq_chunk = sq[i_start:i_end]
        mask_i = sq_chunk[:, None, :] < thr[:, :, None]
        mask_j = sq[None, :, :] < thr[:, :, None]

        inter = (mask_i & mask_j).sum(axis=2).astype(np.float64) + 0.5
        union = (mask_i | mask_j).sum(axis=2).astype(np.float64) + 0.5

        d = np.where(union == 0.5, 1.0, 1.0 - inter / union)

        # Store only upper triangle into condensed vector
        for ci in range(ch):
            i = i_start + ci
            if i >= n - 1:
                break
            k_start = i * n - i * (i + 1) // 2
            n_j = n - i - 1
            condensed[k_start : k_start + n_j] = d[ci, i + 1 :]

    return 100.0 * condensed


def c_extrcoeff_matrix(
    sim_data: np.ndarray,
    madogram: bool = False,
) -> np.ndarray:
    """Compute rank-based extremal coefficient dissimilarity matrix.

    Saunders method using madogram of ranks.  This version uses
    ``scipy.stats.rankdata`` and ``scipy.spatial.distance.cdist`` for
    ~50–100× speedup over the scalar loop implementation.

    Parameters
    ----------
    sim_data : ndarray, shape (nrow, ncol, n_sim)
        Simulated data.
    madogram : bool
        If True, return raw madogram values instead of EC-1.

    Returns
    -------
    ndarray, shape (n_grid, n_grid)
        Symmetric matrix.
    """
    from scipy.spatial.distance import cdist
    from scipy.stats import rankdata

    nrow, ncol, n_sim = sim_data.shape
    n_grid = nrow * ncol
    # Flatten in Fortran (column-major) order to match R's as.vector()
    flat = sim_data.reshape(n_grid, n_sim, order='F')

    # Compute ranks per row — equivalent to the original
    #   rank[s, k] = sum(flat[s, k] <= flat[s, :])
    # For continuous data (no ties), any ranking method gives the same
    # absolute differences, so we use the fast argsort-based rankdata.
    rank_matrix = np.empty((n_grid, n_sim))
    for s in range(n_grid):
        rank_matrix[s] = rankdata(flat[s], method='max')

    # Pairwise mean absolute rank difference via cdist (C-optimised)
    # cdist 'cityblock' = sum|u_k - v_k|;  divide by n_sim for mean
    mado = cdist(rank_matrix, rank_matrix, metric='cityblock')
    mado /= n_sim * 2.0 * (n_sim + 1)

    if madogram:
        return mado

    # Convert to extremal coefficient minus 1
    with np.errstate(divide='ignore', invalid='ignore'):
        ec_matrix = np.minimum(1.0, (1.0 + 2.0 * mado) / (1.0 - 2.0 * mado) - 1.0)
    np.fill_diagonal(ec_matrix, 0.0)
    return ec_matrix


def clustering(
    dist_matrix: np.ndarray,
    method: str = "average",
) -> np.ndarray:
    """Hierarchical agglomerative clustering.

    Parameters
    ----------
    dist_matrix : ndarray, shape (n, n) or (n * (n - 1) / 2,)
        Symmetric dissimilarity matrix or its condensed upper triangle.
    method : str
        Linkage method (default: 'average').

    Returns
    -------
    ndarray
        Linkage matrix (scipy format).
    """
    from scipy.spatial.distance import squareform

    from weatherisk.backend import clustering_via_r

    if dist_matrix.ndim == 1:
        # Convert condensed to square for R
        dist_matrix = squareform(dist_matrix)

    return clustering_via_r(dist_matrix, method=method)


def cluster_number_threshold_method(
    hc: np.ndarray,
    threshold: float,
) -> int:
    """Determine number of clusters by cutting the dendrogram.

    Counts the number of merge heights exceeding the threshold.

    Parameters
    ----------
    hc : ndarray
        Linkage matrix.
    threshold : float
        Height threshold.

    Returns
    -------
    int
        Number of clusters (>= 1).
    """
    heights = hc[:, 2]
    k = int(np.sum(heights > threshold))
    return max(k, 1)


def quantile_threshold(
    dist_matrix: np.ndarray,
    quantile: float = 0.30,
) -> float:
    """Compute a dissimilarity threshold from a percentile of the upper triangle.

    This matches the approach in Contzen et al. (2025, Extremes 28:713–737)
    Section 4, where the 30-th percentile of pairwise dissimilarities is
    used as cut-off to determine the number of clusters.

    Parameters
    ----------
    dist_matrix : ndarray, shape (n, n) or (n * (n - 1) / 2,)
        Symmetric dissimilarity matrix or its condensed upper triangle.
    quantile : float
        Quantile in [0, 1] (default 0.30 = 30th percentile).

    Returns
    -------
    float
        The threshold value.
    """
    if dist_matrix.ndim == 1:
        return float(np.percentile(dist_matrix, quantile * 100))

    n = dist_matrix.shape[0]
    upper_tri = dist_matrix[np.triu_indices(n, k=1)]
    return float(np.percentile(upper_tri, quantile * 100))

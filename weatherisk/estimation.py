"""Local estimate smoothing, cropping, and in-cluster re-estimation.

Spatial moving-average smoothing of local parameter estimates, with
special angular wrapping for the rotation parameter gamma.
"""

from __future__ import annotations

import numpy as np

from weatherisk.grid import Grid


def smooth_local_estimates(
    estimates: np.ndarray,
    smoothing_dist: int,
    grid: Grid,
) -> np.ndarray:
    """Spatially smooth local estimates by moving average.

    Parameters
    ----------
    estimates : ndarray, shape (n_grid, 3)
        Local estimates (a, b, g) per grid point.
    smoothing_dist : int
        Radius in grid cells.  0 = no smoothing.
    grid : Grid
        Spatial grid.

    Returns
    -------
    ndarray, shape (n_grid, 3)
        Smoothed estimates.
    """
    if smoothing_dist == 0:
        return estimates.copy()

    n_grid = grid.n_grid
    result = np.empty_like(estimates)

    # Build neighbour offsets within the smoothing radius
    offsets = []
    for dx in range(-smoothing_dist, smoothing_dist + 1):
        for dy in range(-smoothing_dist, smoothing_dist + 1):
            if np.sqrt(dx * dx + dy * dy) <= smoothing_dist:
                offsets.append((dx, dy))

    for idx in range(n_grid):
        i, j = grid.number_grid(idx)
        # Collect neighbour indices
        neighbours = []
        for dx, dy in offsets:
            ni, nj = i + dy, j + dx
            if 0 <= ni < grid.nrow and 0 <= nj < grid.ncol:
                neighbours.append(grid.grid_number(ni, nj))

        nb_est = estimates[neighbours]

        # Average a and b normally
        result[idx, 0] = nb_est[:, 0].mean()
        result[idx, 1] = nb_est[:, 1].mean()

        # Average gamma with angular wrapping
        centre_g = estimates[idx, 2]
        angles = nb_est[:, 2].copy()
        # Centre angles around the current point's angle
        angles = np.where(
            angles < centre_g - np.pi / 2, angles + np.pi, angles
        )
        angles = np.where(
            angles > centre_g + np.pi / 2, angles - np.pi, angles
        )
        mean_g = angles.mean()
        # Re-centre to [-pi/2, pi/2]
        if mean_g < -np.pi / 2:
            mean_g += np.pi
        elif mean_g > np.pi / 2:
            mean_g -= np.pi
        result[idx, 2] = mean_g

    return result


def crop_matrix(
    values: np.ndarray,
    margin: int,
    grid: Grid,
) -> np.ndarray:
    """Crop a grid-shaped parameter vector by removing edge margins.

    Mirrors R's crop_matrix: reshapes the flat vector into the grid
    matrix (column-major / Fortran order), slices off *margin* rows
    and columns from each edge, then flattens back (column-major).

    Parameters
    ----------
    values : 1-D array of length n_grid
        Flat parameter values in column-major order.
    margin : int
        Number of rows/columns to remove from each edge.
    grid : Grid
        Spatial grid (provides nrow, ncol).

    Returns
    -------
    1-D array
        Cropped values, flattened in column-major order.
    """
    mat = values.reshape((grid.nrow, grid.ncol), order='F')
    cropped = mat[margin:grid.nrow - margin, margin:grid.ncol - margin]
    return cropped.ravel(order='F')


def crop_local_estimates(
    estimates: np.ndarray,
    margin: int,
    grid: Grid,
) -> np.ndarray:
    """Crop local estimates by removing edge margins from each column.

    Mirrors R's crop_local_estimates: applies crop_matrix to each of
    the three parameter columns (a, b, g).

    Parameters
    ----------
    estimates : ndarray, shape (n_grid, 3)
        Local estimates (a, b, g) per grid point.
    margin : int
        Number of rows/columns to remove from each edge.
    grid : Grid
        Spatial grid.

    Returns
    -------
    ndarray, shape (n_cropped, 3)
        Cropped estimates.
    """
    if margin == 0:
        return estimates.copy()
    a_crop = crop_matrix(estimates[:, 0], margin, grid)
    b_crop = crop_matrix(estimates[:, 1], margin, grid)
    g_crop = crop_matrix(estimates[:, 2], margin, grid)
    return np.column_stack([a_crop, b_crop, g_crop])


def calc_estimates_in_clusters(
    sim_data: np.ndarray,
    clusters: np.ndarray,
    df: float,
    alpha: float,
    grid: Grid,
    cluster_ids: list[int] | None = None,
    upper_bounds: tuple[float, float] = (15.0, 15.0),
) -> np.ndarray:
    """Re-estimate (a, b, g) within each cluster.

    Parameters
    ----------
    sim_data : ndarray, shape (nrow, ncol, n_sim)
        Simulation data.
    clusters : 1-D array of length n_grid
        Cluster label per grid point.
    df : float
        Degrees of freedom.
    alpha : float
        Smoothness exponent.
    grid : Grid
        Spatial grid.
    cluster_ids : list[int], optional
        Which clusters to process (default: all).
    upper_bounds : tuple
        Upper bounds for (a, b) in optimisation.

    Returns
    -------
    ndarray, shape (max_cluster+1, 5)
        Columns: a, b, g, n_cells, average_llh.
    """
    from weatherisk.density import pairwise_density_optim

    if cluster_ids is None:
        cluster_ids = list(range(1, int(clusters.max()) + 1))

    max_cl = int(clusters.max())
    results = np.full((max_cl + 1, 5), -np.inf)
    X_flat = grid.X_flat
    Y_flat = grid.Y_flat
    max_dist = 4.0 * ((grid.x_ax.max() - grid.x_ax.min()) / (grid.resolution - 1))

    for cl in cluster_ids:
        which_cl = np.where(clusters == cl)[0]
        if len(which_cl) < 5:
            continue

        # Extract simulation data for cluster cells
        z = np.empty((len(which_cl), sim_data.shape[2]))
        for k, idx in enumerate(which_cl):
            i, j = grid.number_grid(idx)
            z[k, :] = sim_data[i, j, :]

        est = pairwise_density_optim(
            z, df, alpha, X_flat[which_cl], Y_flat[which_cl],
            upper_bounds=upper_bounds, max_dist=max_dist, ensemble=3,
        )
        results[cl, 0:3] = est
        results[cl, 3] = len(which_cl)
        # Compute average log-likelihood (column 4), matching R's calc_estimates_in_clusters
        results[cl, 4] = llh_in_cluster(
            z, df, alpha, X_flat[which_cl], Y_flat[which_cl],
            est, max_dist=max_dist, average=True,
        )

    return results


def llh_in_cluster(
    sim_data: np.ndarray,
    df: float,
    alpha: float,
    X: np.ndarray,
    Y: np.ndarray,
    locest: np.ndarray,
    max_dist: float = 0.0,
    average: bool = False,
) -> float:
    """Log-likelihood of a cluster given its estimated parameters."""
    from weatherisk.backend import neg_log_likelihood_sum as _nll_sum

    n_grid = len(X)
    if sim_data.ndim == 3:
        z = np.empty((n_grid, sim_data.shape[2]))
        for k in range(n_grid):
            z[k, :] = sim_data.ravel()[k::n_grid][:sim_data.shape[2]]
    else:
        z = sim_data

    n_sim = z.shape[1]
    ilist, jlist = np.triu_indices(n_grid, k=1)

    Xlist = np.repeat(X[ilist] - X[jlist], n_sim)
    Ylist = np.repeat(Y[ilist] - Y[jlist], n_sim)

    zilist = []
    zjlist = []
    for k in range(len(ilist)):
        zilist.extend(z[ilist[k], :])
        zjlist.extend(z[jlist[k], :])
    zilist = np.array(zilist)
    zjlist = np.array(zjlist)

    if max_dist > 0:
        sel = Xlist ** 2 + Ylist ** 2 <= max_dist ** 2
        zilist = zilist[sel]
        zjlist = zjlist[sel]
        Xlist = Xlist[sel]
        Ylist = Ylist[sel]

    if len(zilist) == 0:
        return 0.0

    # neg_log_likelihood_sum returns -sum, so negate for log-likelihood
    lh = -_nll_sum(
        zilist, zjlist, Xlist, Ylist, df, alpha,
        locest[0], locest[1], locest[2],
    )
    if average:
        lh /= len(zilist)
    return lh

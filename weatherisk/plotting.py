"""Visualisation: heatmaps, cluster maps, dendrograms, choropleths, bar charts.

All plot functions return a matplotlib Figure object and optionally
display it (show=True).  This allows tests to verify figure creation
without opening a window.
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for tests
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from weatherisk.grid import Grid


def _make_cmap(n: int) -> ListedColormap:
    """Create a qualitative colormap with *n* distinct colours."""
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    return ListedColormap([base[i % len(base)] for i in range(max(n, 1))])


def plot_map(
    data: np.ndarray,
    grid: Grid,
    show: bool = True,
    title: str = "",
    cmap: str = "RdBu_r",
) -> plt.Figure:
    """Plot a gridded heatmap.

    Parameters
    ----------
    data : 2-D array, shape (nrow, ncol)
        Values to plot.
    grid : Grid
        Spatial grid for axis coordinates.
    show : bool
        Whether to display the plot.
    title : str
        Plot title.
    cmap : str
        Matplotlib colourmap name.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(grid.x_ax, grid.y_ax, data, shading="auto", cmap=cmap)
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if show:
        plt.show()
    return fig


def plot_cluster_map(
    clusters: np.ndarray,
    grid: Grid,
    show: bool = True,
    title: str = "Cluster map",
) -> plt.Figure:
    """Plot a cluster map with distinct colours.

    Parameters
    ----------
    clusters : 1-D array of length n_grid
        Cluster labels per grid point.
    grid : Grid
        Spatial grid.
    show : bool
        Whether to display the plot.
    title : str
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    cluster_grid = clusters.reshape(grid.nrow, grid.ncol, order='F').astype(float)
    k = int(clusters.max()) + 1
    cmap = _make_cmap(k)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(
        grid.x_ax, grid.y_ax, cluster_grid, shading="auto", cmap=cmap
    )
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if show:
        plt.show()
    return fig


def plot_dendrogram(
    linkage_matrix: np.ndarray,
    show: bool = True,
    title: str = "Dendrogram",
) -> plt.Figure:
    """Plot a dendrogram from a linkage matrix."""
    from scipy.cluster.hierarchy import dendrogram

    fig, ax = plt.subplots(figsize=(10, 5))
    dendrogram(linkage_matrix, ax=ax)
    ax.set_title(title)
    ax.set_ylabel("Distance")
    if show:
        plt.show()
    return fig


def plot_risk_choropleth(
    cluster_id: np.ndarray,
    metric_values: dict[int, float],
    lons: np.ndarray,
    lats: np.ndarray,
    show: bool = True,
    title: str = "Risk choropleth",
) -> plt.Figure:
    """Plot a choropleth map filled by a per-cluster metric."""
    filled = np.full(cluster_id.shape, np.nan, dtype=float)
    for cid, val in metric_values.items():
        filled[cluster_id == cid] = val

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh(lons, lats, filled, shading="auto", cmap="viridis")
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    if show:
        plt.show()
    return fig


def plot_bar_chart(
    labels: list[str],
    values: list[float],
    show: bool = True,
    title: str = "Bar chart",
    xlabel: str = "",
) -> plt.Figure:
    """Horizontal bar chart."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(labels, values)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    if show:
        plt.show()
    return fig


def plot_cluster_comparison(
    grid: Grid,
    labels_edc: np.ndarray,
    labels_lec: np.ndarray,
    inclusters_edc: np.ndarray,
    inclusters_lec: np.ndarray,
    true_field: np.ndarray,
    param_index: int = 1,
    param_name: str = "b",
    vmin: float | None = None,
    vmax: float | None = None,
    label: str = "",
    suptitle: str = "",
    show: bool = False,
    filename: str | None = None,
) -> plt.Figure:
    """Three-panel cluster comparison figure (cf. Fig. 3 Contzen et al. 2025).

    Panel (a): EDC clusters coloured by in-cluster estimated parameter.
    Panel (b): LEC clusters coloured by in-cluster estimated parameter.
    Panel (c): True parameter field for reference.

    Parameters
    ----------
    grid : Grid
        Spatial grid.
    labels_edc, labels_lec : 1-D arrays of length n_grid
        Cluster labels (1-based, from scipy.cluster.hierarchy.fcluster).
    inclusters_edc, inclusters_lec : ndarray, shape (max_cluster+1, 5)
        In-cluster parameter estimates (columns: a, b, g, n_cells, avg_llh).
    true_field : ndarray, shape (nrow, ncol)
        True parameter values.
    param_index : int
        Column index in inclusters arrays (0=a, 1=b, 2=gamma).
    param_name : str
        Parameter name for axis labels (e.g. 'b').
    vmin, vmax : float, optional
        Colour-scale limits.  Inferred from true_field if None.
    label : str
        Colourbar label.
    suptitle : str
        Figure super-title.
    show : bool
        Whether to call plt.show().
    filename : str, optional
        If given, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if vmin is None:
        vmin = float(np.nanmin(true_field))
    if vmax is None:
        vmax = float(np.nanmax(true_field))
    if not label:
        label = f"${param_name}$"

    cmap = plt.cm.Blues_r

    def _cluster_to_grid(labels, inclusters, pidx):
        mapped = np.full(len(labels), np.nan)
        for i, cl in enumerate(labels):
            if cl < inclusters.shape[0] and np.isfinite(inclusters[cl, pidx]):
                mapped[i] = inclusters[cl, pidx]
        return mapped.reshape(grid.nrow, grid.ncol, order='F')

    edc_map = _cluster_to_grid(labels_edc, inclusters_edc, param_index)
    lec_map = _cluster_to_grid(labels_lec, inclusters_lec, param_index)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # (a) EDC
    im0 = axes[0].pcolormesh(
        grid.x_ax, grid.y_ax, edc_map,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    cl_edc_2d = labels_edc.reshape(grid.nrow, grid.ncol, order='F').astype(float)
    axes[0].contour(
        grid.x_ax, grid.y_ax, cl_edc_2d,
        colors="black", linewidths=0.8, alpha=0.6,
    )
    axes[0].set_title(f"(a) EDC clusters — est. {param_name}", fontsize=10)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    # (b) LEC
    im1 = axes[1].pcolormesh(
        grid.x_ax, grid.y_ax, lec_map,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    cl_lec_2d = labels_lec.reshape(grid.nrow, grid.ncol, order='F').astype(float)
    axes[1].contour(
        grid.x_ax, grid.y_ax, cl_lec_2d,
        colors="black", linewidths=0.8, alpha=0.6,
    )
    axes[1].set_title(f"(b) LEC clusters — est. {param_name}", fontsize=10)
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    # (c) True field
    im2 = axes[2].pcolormesh(
        grid.x_ax, grid.y_ax, true_field,
        shading="auto", cmap=cmap, vmin=vmin, vmax=vmax,
    )
    axes[2].set_title(f"(c) True {param_name}", fontsize=10)
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("y")
    axes[2].set_aspect("equal")

    # Shared colourbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    plt.colorbar(im2, cax=cbar_ax, label=label)

    if suptitle:
        fig.suptitle(suptitle, fontsize=13, y=1.02)

    if filename:
        fig.savefig(filename, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig
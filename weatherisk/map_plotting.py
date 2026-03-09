"""Geographic map visualisation with Cartopy — filled-region style.

Provides Cartopy-based plotting functions for cluster maps, parameter
fields, risk choropleths, and multi-panel summary figures.  Uses
``pcolormesh`` to fill each grid cell with a solid colour (matching the
style of Justus 2025, Extremes, Figure 9) rather than scatter dots.

All plot functions return a :class:`matplotlib.figure.Figure` and
optionally save to PDF.  Designed for use with the CPC real-data
pipeline output.

Requires ``cartopy`` (>= 0.21) at runtime.

Example
-------
>>> from weatherisk.map_plotting import plot_cluster_map_geo
>>> fig = plot_cluster_map_geo(
...     lat, lon, labels, k=15,
...     lats_1d=lats_1d, lons_1d=lons_1d,
...     n_lat=n_lat, n_lon=n_lon, land_idx=land_idx,
...     title="LEC Clusters", save_path="lec.pdf",
... )
"""

from __future__ import annotations

import os
from typing import Sequence

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm


# ── colour helpers ─────────────────────────────────────────────

def _qualitative_cmap(k: int) -> ListedColormap:
    """Return a qualitative ListedColormap with *k* distinct colours."""
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    return ListedColormap([base[i % len(base)] for i in range(max(k, 1))])


# ── grid reconstruction ───────────────────────────────────────

def _to_grid(
    values: np.ndarray,
    lats_1d: np.ndarray,
    lons_1d: np.ndarray,
    n_lat: int,
    n_lon: int,
    land_idx: np.ndarray,
    fill: float = np.nan,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct a 2-D grid from land-only 1-D values.

    Parameters
    ----------
    values : (n_land,) array
        Value per valid land cell.
    lats_1d, lons_1d : 1-D arrays
        Latitude / longitude coordinate vectors of the coarsened grid.
    n_lat, n_lon : int
        Shape of the full grid.
    land_idx : (n_land,) int array
        Flat indices of the valid cells within the (n_lat × n_lon) grid.
    fill : float
        Value for ocean / missing cells (default NaN → transparent).

    Returns
    -------
    grid : (n_lat, n_lon) array
        2-D grid with *values* placed at *land_idx*.
    lat_edges : (n_lat + 1,) array
        Cell boundary latitudes for ``pcolormesh``.
    lon_edges : (n_lon + 1,) array
        Cell boundary longitudes for ``pcolormesh``.
    """
    grid = np.full(n_lat * n_lon, fill, dtype=float)
    grid[land_idx] = values
    grid = grid.reshape(n_lat, n_lon)

    def _edges(centres):
        c = np.asarray(centres, dtype=float)
        d = np.diff(c)
        edges = np.empty(len(c) + 1)
        edges[1:-1] = c[:-1] + d / 2
        edges[0] = c[0] - d[0] / 2
        edges[-1] = c[-1] + d[-1] / 2
        return edges

    return grid, _edges(lats_1d), _edges(lons_1d)


# ── base-axis constructor ─────────────────────────────────────

def _base_ax(
    fig: plt.Figure,
    extent: tuple[float, float, float, float],
    pos: int | tuple[int, int, int] = 111,
):
    """Add a GeoAxes with coastlines, borders and ocean shading."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent(list(extent), crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
    ax.add_feature(cfeature.LAND, facecolor="wheat", alpha=0.10)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
    return ax


def _add_gridlines(ax) -> None:
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False


# ── helper: pcolormesh on a GeoAxes ───────────────────────────

def _filled_plot(
    ax,
    grid: np.ndarray,
    lat_edges: np.ndarray,
    lon_edges: np.ndarray,
    *,
    cmap,
    norm=None,
    vmin=None,
    vmax=None,
    alpha: float = 1.0,
):
    """Draw *grid* with ``pcolormesh`` and return the QuadMesh."""
    import cartopy.crs as ccrs

    masked = np.ma.masked_invalid(grid)
    return ax.pcolormesh(
        lon_edges, lat_edges, masked,
        cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
        transform=ccrs.PlateCarree(), alpha=alpha,
        shading="flat",
    )


# ── public plotting functions ─────────────────────────────────

def plot_cluster_map_geo(
    lat: np.ndarray,
    lon: np.ndarray,
    labels: np.ndarray,
    k: int,
    *,
    lats_1d: np.ndarray | None = None,
    lons_1d: np.ndarray | None = None,
    n_lat: int | None = None,
    n_lon: int | None = None,
    land_idx: np.ndarray | None = None,
    extent: tuple[float, float, float, float] = (3, 57, 28, 67),
    title: str = "Clusters",
    save_path: str | None = None,
    dpi: int = 300,
) -> plt.Figure:
    """Filled-region cluster map on a Cartopy projection.

    Uses ``pcolormesh`` when grid metadata is provided, falling back to
    scatter dots otherwise.

    Parameters
    ----------
    lat, lon : 1-D arrays
        Geographic coordinates of each cell.
    labels : 1-D int array
        Cluster label per cell (1-based).
    k : int
        Number of clusters.
    lats_1d, lons_1d : optional 1-D arrays
        Coordinate vectors of the regular grid.
    n_lat, n_lon : optional int
        Grid dimensions.
    land_idx : optional 1-D int array
        Flat indices of valid cells.
    extent : (lon_min, lon_max, lat_min, lat_max)
        Map extent in degrees.
    title : str
        Plot title.
    save_path : str, optional
        If given, save figure as PDF/PNG.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(10, 7))
    ax = _base_ax(fig, extent)
    cmap = _qualitative_cmap(k)
    norm = BoundaryNorm(np.arange(0.5, k + 1.5), cmap.N)

    _has_grid = all(x is not None for x in (lats_1d, lons_1d, n_lat, n_lon, land_idx))
    if _has_grid:
        grid, lat_e, lon_e = _to_grid(
            labels.astype(float), lats_1d, lons_1d, n_lat, n_lon, land_idx)
        mesh = _filled_plot(ax, grid, lat_e, lon_e, cmap=cmap, norm=norm)
    else:
        mesh = ax.scatter(
            lon, lat, c=labels, cmap=cmap, norm=norm,
            s=45, edgecolors="k", linewidths=0.3,
            transform=ccrs.PlateCarree(), zorder=5,
        )

    plt.colorbar(mesh, ax=ax, shrink=0.7, label="Cluster",
                 ticks=range(1, k + 1))
    ax.set_title(title, fontsize=11, fontweight="bold")
    _add_gridlines(ax)

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_field_map(
    lat: np.ndarray,
    lon: np.ndarray,
    values: np.ndarray,
    *,
    lats_1d: np.ndarray | None = None,
    lons_1d: np.ndarray | None = None,
    n_lat: int | None = None,
    n_lon: int | None = None,
    land_idx: np.ndarray | None = None,
    extent: tuple[float, float, float, float] = (3, 57, 28, 67),
    title: str = "",
    label: str = "",
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    save_path: str | None = None,
    dpi: int = 300,
) -> plt.Figure:
    """Filled continuous-field map on a Cartopy projection.

    Parameters
    ----------
    lat, lon : 1-D arrays
    values : 1-D array
        Scalar value per cell.
    lats_1d, lons_1d, n_lat, n_lon, land_idx :
        Grid metadata for pcolormesh (optional — scatter fallback).
    extent, title, label, cmap, vmin, vmax, save_path, dpi :
        Display parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(10, 7))
    ax = _base_ax(fig, extent)

    _has_grid = all(x is not None for x in (lats_1d, lons_1d, n_lat, n_lon, land_idx))
    if _has_grid:
        grid, lat_e, lon_e = _to_grid(
            values, lats_1d, lons_1d, n_lat, n_lon, land_idx)
        mesh = _filled_plot(ax, grid, lat_e, lon_e, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        mesh = ax.scatter(
            lon, lat, c=values, cmap=cmap,
            vmin=vmin, vmax=vmax,
            s=45, edgecolors="k", linewidths=0.3,
            transform=ccrs.PlateCarree(), zorder=5,
        )

    plt.colorbar(mesh, ax=ax, shrink=0.7, label=label)
    ax.set_title(title, fontsize=11, fontweight="bold")
    _add_gridlines(ax)

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_risk_map(
    lat: np.ndarray,
    lon: np.ndarray,
    labels: np.ndarray,
    risk_records: list[dict],
    *,
    metric: str = "es",
    lats_1d: np.ndarray | None = None,
    lons_1d: np.ndarray | None = None,
    n_lat: int | None = None,
    n_lon: int | None = None,
    land_idx: np.ndarray | None = None,
    extent: tuple[float, float, float, float] = (3, 57, 28, 67),
    title: str = "Expected Shortfall ES$_{95}$",
    save_path: str | None = None,
    dpi: int = 300,
) -> plt.Figure:
    """Filled choropleth colouring cells by a per-cluster risk metric.

    Parameters
    ----------
    lat, lon : 1-D arrays
    labels : 1-D array
        Cluster label per cell.
    risk_records : list of dict
        Each dict has keys ``cluster``, ``var``, ``es``.
    metric : str
        Key in *risk_records* to plot (``"var"`` or ``"es"``).
    lats_1d, lons_1d, n_lat, n_lon, land_idx :
        Grid metadata for pcolormesh (optional — scatter fallback).
    extent, title, save_path, dpi :
        Display parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import cartopy.crs as ccrs

    rmap = np.zeros(len(labels), dtype=float)
    for r in risk_records:
        rmap[labels == r["cluster"]] = r[metric]

    fig = plt.figure(figsize=(10, 7))
    ax = _base_ax(fig, extent)

    _has_grid = all(x is not None for x in (lats_1d, lons_1d, n_lat, n_lon, land_idx))
    if _has_grid:
        grid, lat_e, lon_e = _to_grid(
            rmap, lats_1d, lons_1d, n_lat, n_lon, land_idx)
        mesh = _filled_plot(ax, grid, lat_e, lon_e, cmap="YlOrRd")
    else:
        mesh = ax.scatter(
            lon, lat, c=rmap, cmap="YlOrRd",
            s=45, edgecolors="k", linewidths=0.3,
            transform=ccrs.PlateCarree(), zorder=5,
        )

    cb_label = ("ES$_{95}$" if metric == "es" else "VaR$_{95}$") + "  (Fr\u00e9chet scale)"
    plt.colorbar(mesh, ax=ax, shrink=0.7, label=cb_label)
    ax.set_title(title, fontsize=12, fontweight="bold")
    _add_gridlines(ax)

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_summary_panel(
    lat: np.ndarray,
    lon: np.ndarray,
    labels_lec: np.ndarray,
    labels_edc: np.ndarray,
    k_lec: int,
    k_edc: int,
    smoothed: np.ndarray,
    risk_lec: list[dict],
    risk_edc: list[dict],
    *,
    lats_1d: np.ndarray | None = None,
    lons_1d: np.ndarray | None = None,
    n_lat: int | None = None,
    n_lon: int | None = None,
    land_idx: np.ndarray | None = None,
    extent: tuple[float, float, float, float] = (3, 57, 28, 67),
    suptitle: str = "Climate Risk — LEC / EDC on CPC Precipitation  (2000–2019)",
    save_path: str | None = None,
    dpi: int = 300,
) -> plt.Figure:
    """2×3 summary panel: clusters, parameters, risk.

    Layout::

        LEC clusters | EDC clusters | semi-minor a
        anisotropy b | ES₉₅  LEC   | ES₉₅  EDC

    Parameters
    ----------
    lat, lon : 1-D arrays, geographic coordinates.
    labels_lec, labels_edc : 1-D int arrays, cluster assignments.
    k_lec, k_edc : int, cluster counts.
    smoothed : (n, 3) array of (a, b, γ) per cell.
    risk_lec, risk_edc : lists of per-cluster risk dicts.
    lats_1d, lons_1d, n_lat, n_lon, land_idx :
        Grid metadata for pcolormesh (optional — scatter fallback).
    extent, suptitle, save_path, dpi : display parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    _has_grid = all(x is not None for x in (lats_1d, lons_1d, n_lat, n_lon, land_idx))

    fig = plt.figure(figsize=(20, 12))
    ext = list(extent)

    def _ax(pos):
        ax = fig.add_subplot(2, 3, pos, projection=ccrs.PlateCarree())
        ax.set_extent(ext, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.OCEAN, facecolor="lightskyblue", alpha=0.3)
        ax.add_feature(cfeature.LAND, facecolor="wheat", alpha=0.10)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")
        return ax

    def _draw_clusters(ax, labels_arr, k_val):
        cm = _qualitative_cmap(k_val)
        nm = BoundaryNorm(np.arange(0.5, k_val + 1.5), cm.N)
        if _has_grid:
            grid, lat_e, lon_e = _to_grid(
                labels_arr.astype(float), lats_1d, lons_1d, n_lat, n_lon, land_idx)
            _filled_plot(ax, grid, lat_e, lon_e, cmap=cm, norm=nm)
        else:
            ax.scatter(lon, lat, c=labels_arr, cmap=cm, norm=nm, s=22,
                       edgecolors="k", linewidths=0.2,
                       transform=ccrs.PlateCarree(), zorder=5)

    def _draw_field(ax, vals, cmap_name, label_text, vmin=None, vmax=None):
        if _has_grid:
            grid, lat_e, lon_e = _to_grid(
                vals, lats_1d, lons_1d, n_lat, n_lon, land_idx)
            sc = _filled_plot(ax, grid, lat_e, lon_e, cmap=cmap_name, vmin=vmin, vmax=vmax)
        else:
            sc = ax.scatter(lon, lat, c=vals, cmap=cmap_name, s=22,
                            edgecolors="k", linewidths=0.2, vmin=vmin, vmax=vmax,
                            transform=ccrs.PlateCarree(), zorder=5)
        plt.colorbar(sc, ax=ax, shrink=0.55, label=label_text)

    def _draw_risk(ax, labels_arr, risk_list, label_text):
        rmap = np.zeros(len(labels_arr), dtype=float)
        for r in risk_list:
            rmap[labels_arr == r["cluster"]] = r["es"]
        if _has_grid:
            grid, lat_e, lon_e = _to_grid(
                rmap, lats_1d, lons_1d, n_lat, n_lon, land_idx)
            sc = _filled_plot(ax, grid, lat_e, lon_e, cmap="YlOrRd")
        else:
            sc = ax.scatter(lon, lat, c=rmap, cmap="YlOrRd", s=22,
                            edgecolors="k", linewidths=0.2,
                            transform=ccrs.PlateCarree(), zorder=5)
        plt.colorbar(sc, ax=ax, shrink=0.55, label=label_text)

    # (1) LEC clusters
    ax = _ax(1)
    _draw_clusters(ax, labels_lec, k_lec)
    ax.set_title(f"LEC Clusters (k={k_lec})", fontsize=11, fontweight="bold")

    # (2) EDC clusters
    ax = _ax(2)
    _draw_clusters(ax, labels_edc, k_edc)
    ax.set_title(f"EDC Clusters (k={k_edc})", fontsize=11, fontweight="bold")

    # (3) parameter a
    ax = _ax(3)
    _draw_field(ax, smoothed[:, 0], "viridis", "a  (dependence range)")
    ax.set_title("Dependence Range (a)", fontsize=11, fontweight="bold")

    # (4) parameter b
    ax = _ax(4)
    _draw_field(ax, smoothed[:, 1], "magma", "b  (anisotropy)")
    ax.set_title("Anisotropy (b)", fontsize=11, fontweight="bold")

    # (5) ES₉₅ per LEC cluster
    ax = _ax(5)
    _draw_risk(ax, labels_lec, risk_lec, "ES$_{95}$ (Fr\u00e9chet)")
    ax.set_title("Tail-Risk ES$_{95}$ (LEC)", fontsize=11, fontweight="bold")

    # (6) ES₉₅ per EDC cluster
    ax = _ax(6)
    _draw_risk(ax, labels_edc, risk_edc, "ES$_{95}$ (Fr\u00e9chet)")
    ax.set_title("Tail-Risk ES$_{95}$ (EDC)", fontsize=11, fontweight="bold")

    plt.suptitle(suptitle, fontsize=14, fontweight="bold", y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fig

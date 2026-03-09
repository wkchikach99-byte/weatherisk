"""Grid construction and coordinate-index conversion utilities.

Provides the :class:`Grid` class for creating regular 2-D grids, converting
between ``(i, j)`` matrix indices, linear indices, and ``(x, y)`` coordinates.
Mirrors the R functions ``grid_number``, ``number_grid``, ``koord_num``, and
``number_koord`` from *functions.R*.
"""

from __future__ import annotations

import numpy as np


def rad(degrees: float) -> float:
    """Convert degrees to radians."""
    return degrees * np.pi / 180.0


def deg(radians: float) -> float:
    """Convert radians to degrees."""
    return 180.0 * radians / np.pi


class Grid:
    """A regular 2-D grid with coordinate and index helpers.

    Parameters
    ----------
    x_range : tuple[float, float]
        ``(x_min, x_max)`` extent of the horizontal axis.
    y_range : tuple[float, float]
        ``(y_min, y_max)`` extent of the vertical axis.
    resolution : int
        Number of grid points along each axis.
    """

    def __init__(
        self,
        x_range: tuple[float, float] = (-5, 5),
        y_range: tuple[float, float] = (-5, 5),
        resolution: int = 10,
    ) -> None:
        self.resolution = resolution
        self.x_range = x_range
        self.y_range = y_range

        # Axes — x runs low→high, y runs high→low (R convention)
        self.x_ax: np.ndarray = np.linspace(x_range[0], x_range[1], resolution)
        self.y_ax: np.ndarray = -np.linspace(y_range[0], y_range[1], resolution)

        self.nrow: int = resolution
        self.ncol: int = resolution
        self.n_grid: int = self.nrow * self.ncol

        # Meshgrids (nrow × ncol)
        self.X: np.ndarray
        self.Y: np.ndarray
        self.X, self.Y = np.meshgrid(self.x_ax, self.y_ax)

    # ------------------------------------------------------------------
    # Index conversions (all 0-based, matching NumPy convention)
    # ------------------------------------------------------------------

    def grid_number(self, i: int, j: int) -> int:
        """Convert ``(row, col)`` indices to a linear index (0-based).

        Raises
        ------
        IndexError
            If *i* or *j* is out of bounds.
        """
        if i < 0 or j < 0 or i >= self.nrow or j >= self.ncol:
            raise IndexError(
                f"Index ({i}, {j}) out of bounds for grid "
                f"{self.nrow}×{self.ncol}"
            )
        return j * self.nrow + i

    def number_grid(self, n: int) -> tuple[int, int]:
        """Convert linear index *n* (0-based) to ``(row, col)``.

        Raises
        ------
        IndexError
            If *n* is out of bounds.
        """
        if n < 0 or n >= self.n_grid:
            raise IndexError(
                f"Index {n} out of bounds for grid size {self.n_grid}"
            )
        i = n % self.nrow
        j = n // self.nrow
        return i, j

    @property
    def X_flat(self) -> np.ndarray:
        """X coordinates flattened in column-major order (matching ``grid_number``)."""
        return self.X.flatten(order="F")

    @property
    def Y_flat(self) -> np.ndarray:
        """Y coordinates flattened in column-major order (matching ``grid_number``)."""
        return self.Y.flatten(order="F")

    def koord_num(self, x: float, y: float) -> int:
        """Return the linear index of the nearest grid point to ``(x, y)``.

        The returned index uses column-major ordering, consistent with
        ``grid_number`` and ``number_grid``.
        """
        dists = (self.X - x) ** 2 + (self.Y - y) ** 2
        row, col = np.unravel_index(int(np.argmin(dists)), self.X.shape)
        return self.grid_number(row, col)

    def number_koord(self, idx: int) -> tuple[float, float]:
        """Return the ``(x, y)`` coordinates of linear index *idx*.

        *idx* uses column-major ordering, consistent with
        ``grid_number`` and ``number_grid``.
        """
        i, j = self.number_grid(idx)
        return float(self.X[i, j]), float(self.Y[i, j])

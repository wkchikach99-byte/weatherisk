"""Memory-efficient estimation and clustering for large global grids.

Provides coarse-grid proxy clustering, chunk-based parallel estimation,
and checkpoint/resume for HPC job arrays.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def downsample_estimates(
    estimates: np.ndarray,
    fine_shape: tuple[int, int],
    coarse_shape: tuple[int, int],
) -> np.ndarray:
    """Downsample local estimates by spatial block averaging.

    Parameters
    ----------
    estimates : ndarray, shape (n_fine, n_params)
        Fine-grid estimates (e.g. (400, 3) for a 20x20 grid).
    fine_shape : tuple
        (n_rows, n_cols) of the fine grid.
    coarse_shape : tuple
        (n_rows, n_cols) of the coarse grid.

    Returns
    -------
    ndarray, shape (n_coarse, n_params)
    """
    n_params = estimates.shape[1]
    fine_r, fine_c = fine_shape
    coarse_r, coarse_c = coarse_shape

    # Reshape to 2-D grid
    est_grid = estimates.reshape(fine_r, fine_c, n_params)

    # Block sizes
    br = fine_r // coarse_r
    bc = fine_c // coarse_c

    coarse = np.empty((coarse_r, coarse_c, n_params))
    for i in range(coarse_r):
        for j in range(coarse_c):
            block = est_grid[i * br : (i + 1) * br, j * bc : (j + 1) * bc, :]
            coarse[i, j, :] = block.mean(axis=(0, 1))

    return coarse.reshape(coarse_r * coarse_c, n_params)


def propagate_cluster_labels(
    coarse_labels: np.ndarray,
    coarse_shape: tuple[int, int],
    fine_shape: tuple[int, int],
) -> np.ndarray:
    """Assign fine-grid cells to the cluster of their nearest coarse cell.

    Parameters
    ----------
    coarse_labels : 1-D array of length n_coarse
        Cluster labels on the coarse grid.
    coarse_shape : tuple
        (n_rows, n_cols) of the coarse grid.
    fine_shape : tuple
        (n_rows, n_cols) of the fine grid.

    Returns
    -------
    1-D array of length n_fine
    """
    coarse_r, coarse_c = coarse_shape
    fine_r, fine_c = fine_shape

    label_grid = coarse_labels.reshape(coarse_r, coarse_c)

    fine_labels = np.empty(fine_r * fine_c, dtype=coarse_labels.dtype)
    for fi in range(fine_r):
        for fj in range(fine_c):
            ci = min(int(fi * coarse_r / fine_r), coarse_r - 1)
            cj = min(int(fj * coarse_c / fine_c), coarse_c - 1)
            fine_labels[fi * fine_c + fj] = label_grid[ci, cj]

    return fine_labels


def chunk_indices(
    n_total: int,
    n_chunks: int,
) -> list[tuple[int, int]]:
    """Split a range into roughly equal chunks.

    Parameters
    ----------
    n_total : int
        Total number of items.
    n_chunks : int
        Number of chunks.

    Returns
    -------
    list[tuple[int, int]]
        List of (start, end) pairs.
    """
    chunk_size = n_total // n_chunks
    remainder = n_total % n_chunks
    chunks = []
    start = 0
    for i in range(n_chunks):
        end = start + chunk_size + (1 if i < remainder else 0)
        chunks.append((start, end))
        start = end
    return chunks


def save_chunk(
    data: np.ndarray,
    chunk_id: int,
    output_dir: str,
) -> None:
    """Save a chunk of estimation results to disk."""
    path = Path(output_dir) / f"chunk_{chunk_id:04d}.npy"
    np.save(path, data)


def load_chunk(
    chunk_id: int,
    output_dir: str,
) -> np.ndarray:
    """Load a chunk of estimation results from disk."""
    path = Path(output_dir) / f"chunk_{chunk_id:04d}.npy"
    return np.load(path)

"""Tests for weatherisk.scalable — memory-efficient estimation and clustering."""

import numpy as np
import pytest


class TestCoarseGridProxy:
    def test_downsample_preserves_shape(self):
        from weatherisk.scalable import downsample_estimates

        # 20x20 grid, 3 parameters per cell
        estimates = np.random.default_rng(0).random((400, 3))
        coarse = downsample_estimates(
            estimates, fine_shape=(20, 20), coarse_shape=(5, 5)
        )
        assert coarse.shape == (25, 3)

    def test_propagate_labels(self):
        from weatherisk.scalable import propagate_cluster_labels

        # 4x4 fine grid, 2x2 coarse grid with 2 clusters
        coarse_labels = np.array([0, 0, 1, 1])
        fine_labels = propagate_cluster_labels(
            coarse_labels, coarse_shape=(2, 2), fine_shape=(4, 4)
        )
        assert fine_labels.shape == (16,)
        assert set(fine_labels) == {0, 1}


class TestParallelEstimation:
    def test_chunk_indices(self):
        from weatherisk.scalable import chunk_indices

        chunks = chunk_indices(n_total=100, n_chunks=4)
        assert len(chunks) == 4
        all_indices = []
        for start, end in chunks:
            all_indices.extend(range(start, end))
        assert sorted(all_indices) == list(range(100))


class TestCheckpointing:
    def test_save_and_load_chunk(self, tmp_path):
        from weatherisk.scalable import save_chunk, load_chunk

        data = np.random.default_rng(0).random((50, 3))
        save_chunk(data, chunk_id=0, output_dir=str(tmp_path))
        loaded = load_chunk(chunk_id=0, output_dir=str(tmp_path))
        np.testing.assert_array_equal(data, loaded)

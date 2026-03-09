"""Integration test — end-to-end pipeline on a tiny grid."""

import numpy as np
import pytest


@pytest.mark.slow
class TestPipelineIntegration:
    def test_full_run_tiny_grid(self, tmp_path):
        """Run the full pipeline on a 5×5 grid with 5 simulations."""
        from weatherisk.pipeline import run_pipeline

        result = run_pipeline(
            resolution=5,
            n_sim=5,
            df=5,
            alpha=1.0,
            seed=42,
            output_dir=str(tmp_path),
        )
        assert result is not None
        assert "clusters" in result
        assert len(result["clusters"]) == 25

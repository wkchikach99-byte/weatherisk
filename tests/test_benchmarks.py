from pathlib import Path

from weatherisk.benchmarks import HotPathBenchmarkConfig, run_hotpath_benchmark


def test_hotpath_benchmark_runs_and_writes_markdown(tmp_path: Path):
    output = tmp_path / "benchmark_results.md"
    result = run_hotpath_benchmark(
        HotPathBenchmarkConfig(
            n_years=12,
            n_lat=4,
            n_lon=4,
            n_workers=2,
            neighbor_radius=2.5,
        ),
        markdown_path=output,
    )

    assert result["total_seconds"] > 0
    assert result["derived"]["n_valid_cells"] > 0
    assert output.exists()
    text = output.read_text(encoding="utf-8")
    assert "# Benchmark Results" in text
    assert "step3_local_mle" in text
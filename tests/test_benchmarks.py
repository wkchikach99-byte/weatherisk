from pathlib import Path

from weatherisk.benchmarks import (
    Figure9BenchmarkConfig,
    run_figure9_benchmark,
)


def test_figure9_benchmark_runs_and_writes_markdown(tmp_path: Path):
    output = tmp_path / "benchmark_results.md"
    result = run_figure9_benchmark(
        Figure9BenchmarkConfig(
            n_years=12,
            n_lat=4,
            n_lon=4,
            n_workers=2,
            generate_plots=False,
        ),
        markdown_path=output,
    )

    assert result["total_seconds"] > 0
    assert result["derived"]["n_valid_cells"] > 0
    assert result["script_entrypoint"] == "scripts.reproduce_fig9.main"
    assert result["pipeline_entrypoint"] == "weatherisk.cmip6_pipeline.run_cmip6_pipeline"
    assert result["memory_profile"]["metric"] == "peak_process_tree_rss"
    assert result["memory_profile"]["max_rss_bytes"] > 0
    assert output.exists()
    text = output.read_text(encoding="utf-8")
    assert "# Benchmark Results" in text
    assert "scripts.reproduce_fig9.main" in text
    assert "_run_local_estimation_cmip6" in text
    assert "Max memory" in text
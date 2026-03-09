from pathlib import Path

from weatherisk.benchmarks import (
    DecisionBenchmarkConfig,
    HotPathBenchmarkConfig,
    run_decision_benchmark,
    run_hotpath_benchmark,
)


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
    assert result["entrypoint"] == "weatherisk.cmip6_pipeline.run_cmip6_pipeline"
    assert result["memory_profile"]["metric"] == "peak_process_tree_rss"
    assert result["memory_profile"]["max_rss_bytes"] > 0
    assert output.exists()
    text = output.read_text(encoding="utf-8")
    assert "# Benchmark Results" in text
    assert "run_cmip6_pipeline" in text
    assert "_run_local_estimation_cmip6" in text
    assert "Max memory" in text


def test_decision_benchmark_runs_and_writes_markdown(tmp_path: Path):
    output = tmp_path / "decision_benchmark_results.md"
    result = run_decision_benchmark(
        DecisionBenchmarkConfig(
            benchmark=HotPathBenchmarkConfig(
                n_years=16,
                n_lat=6,
                n_lon=6,
                n_workers=2,
                neighbor_radius=2.5,
            ),
            warmup_runs=1,
            measured_runs=2,
        ),
        markdown_path=output,
    )

    assert result["benchmark_method"]["warmup_runs"] == 1
    assert result["benchmark_method"]["measured_runs"] == 2
    assert len(result["warmup_runs"]) == 1
    assert len(result["measured_runs"]) == 2
    assert result["summary"]["checks_stable"] is True
    assert result["summary"]["total_seconds"]["mean_seconds"] > 0
    assert result["summary"]["memory_profile"]["max_max_rss_bytes"] > 0
    assert output.exists()
    text = output.read_text(encoding="utf-8")
    assert "Decision Benchmark" in text
    assert "1 warmup + 2 measured runs" in text
    assert "Mean (s)" in text
"""Command-line interface entry point.

Provides the `weatherisk` CLI with subcommands:
- validate: run the method validation flow on synthetic data
- risk: run the real-data risk analysis flow
- risk-pipeline: lightweight post-processing of pre-computed risk maps
"""

from __future__ import annotations

import click


@click.group()
@click.version_option(package_name="weatherisk")
def main() -> None:
    """weatherisk -- Climate risk analysis via max-stable process clustering."""


@main.command()
@click.option("--params", "preset", default="stripes", help="Parameter preset name.")
@click.option("--resolution", default=10, help="Grid resolution.")
@click.option("--n-sim", default=10, help="Number of simulations.")
@click.option("--seed", default=42, help="Random seed.")
@click.option("--output-dir", default="output/validate", help="Output directory.")
def validate(preset: str, resolution: int, n_sim: int, seed: int, output_dir: str) -> None:
    """Flow 1: Method validation on synthetic data."""
    from weatherisk.parameters import get_preset
    from weatherisk.pipeline import run_pipeline

    p = get_preset(preset)
    click.echo(f"Running validation with preset '{preset}' (res={resolution}, n_sim={n_sim})")
    result = run_pipeline(
        resolution=resolution,
        n_sim=n_sim,
        df=p.df,
        alpha=p.alpha,
        seed=seed,
        output_dir=output_dir,
    )
    n_clusters = len(set(result["clusters"]))
    click.echo(f"Done. Found {n_clusters} clusters. Results saved to {output_dir}")


@main.command()
@click.option("--netcdf", "netcdf_path", default=None, help="Path to NetCDF file.")
@click.option("--hazard", default="heat", help="Hazard type (heat|cold).")
@click.option("-k", default=25, help="Target number of clusters.")
@click.option("--workers", default=1, help="Number of parallel workers.")
@click.option("--output-dir", default="output/risk", help="Output directory.")
def risk(netcdf_path: str | None, hazard: str, k: int, workers: int, output_dir: str) -> None:
    """Flow 2: Real-data risk analysis."""
    click.echo(f"Risk analysis: hazard={hazard}, k={k}, workers={workers}")
    if netcdf_path is None:
        click.echo("Error: --netcdf path is required.", err=True)
        raise SystemExit(1)
    click.echo("Not yet implemented for real data. Use 'validate' for synthetic.")


@main.command(name="risk-pipeline")
@click.option("--csv", "csv_path", default="data/risk_map_grid.csv", help="Path to risk map CSV.")
@click.option("--bands", default=6, help="Number of ES quantile bands.")
@click.option("--sigma", default=0.8, help="Gaussian smoothing sigma.")
@click.option("--min-patch", default=30, help="Minimum patch size.")
@click.option("--output-dir", default="output/risk_pipeline", help="Output directory.")
def risk_pipeline(
    csv_path: str, bands: int, sigma: float, min_patch: int, output_dir: str
) -> None:
    """Flow 3: Lightweight post-processing of pre-computed risk maps."""
    import os
    import numpy as np
    import pandas as pd
    from weatherisk.risk_pipeline import (
        load_and_grid,
        smooth_field,
        quantile_bands,
        connected_patches,
        merge_tiny_regions,
        remap_ids_to_sequential,
        compute_cluster_stats,
    )

    os.makedirs(output_dir, exist_ok=True)

    click.echo(f"Loading {csv_path}...")
    data = load_and_grid(csv_path)

    ES_s = smooth_field(data["ES"], sigma, data["land_mask"])
    es_band, _ = quantile_bands(ES_s, bands)
    cluster_id = connected_patches(es_band, min_patch)
    cluster_id = merge_tiny_regions(cluster_id, data["lon_grid"], data["lat_grid"])
    cluster_id, K = remap_ids_to_sequential(cluster_id)

    out_df = pd.DataFrame({
        "lat": data["lat_grid"].ravel(),
        "lon": data["lon_grid"].ravel(),
        "cluster": cluster_id.ravel(),
    })
    out_csv = os.path.join(output_dir, "clusters.csv")
    out_df.to_csv(out_csv, index=False)

    stats = compute_cluster_stats(out_df, data["df"])
    stats.to_csv(os.path.join(output_dir, "cluster_stats.csv"), index=False)

    click.echo(f"Done. {K} clusters. Output in {output_dir}")


@main.command()
@click.option(
    "--data-dir", default="data/netcdf",
    help="Directory containing CPC NetCDF files.",
)
@click.option(
    "--output-dir", default="docs/figures",
    help="Directory for output PDFs.",
)
@click.option(
    "--variable", default="precip",
    help="NetCDF variable name (precip | tmax).",
)
@click.option(
    "--file-prefix", default=None,
    help="File name prefix (default: same as --variable).",
)
@click.option("--year-start", default=2000, help="First year (inclusive).")
@click.option("--year-end", default=2020, help="Last year (exclusive).")
@click.option(
    "--lat-range", nargs=2, type=float, default=(30.0, 65.0),
    help="Latitude range (min max).",
)
@click.option(
    "--lon-range", nargs=2, type=float, default=(5.0, 55.0),
    help="Longitude range (min max).",
)
@click.option("--coarsen", default=4, help="Coarsening factor (every Nth cell).")
@click.option("--df", "degrees_of_freedom", default=5.0, help="ν (degrees of freedom).")
@click.option("--alpha", default=1.0, help="α (smoothness exponent).")
@click.option("--neighbor-radius", default=3.0, help="ε (normalised MLE neighbourhood).")
@click.option("--smoothing-radius", default=2.0, help="Spatial smoothing radius.")
@click.option("--quantile-threshold", default=0.30, help="Quantile for k-selection.")
@click.option("--risk-level", default=0.95, help="p for VaR / ES.")
@click.option("--dpi", default=300, help="Output figure DPI.")
@click.option(
    "--gdp-path", default=None,
    help="Path to Kummu et al. GDP PPP NetCDF for exposure weighting.",
)
@click.option("--gdp-year", default=2015, help="GDP snapshot year (1990–2015).")
@click.option("--no-plots", is_flag=True, help="Run pipeline only, skip map generation.")
@click.option("--quiet", is_flag=True, help="Suppress progress output.")
def maps(
    data_dir: str,
    output_dir: str,
    variable: str,
    file_prefix: str | None,
    year_start: int,
    year_end: int,
    lat_range: tuple[float, float],
    lon_range: tuple[float, float],
    coarsen: int,
    degrees_of_freedom: float,
    alpha: float,
    neighbor_radius: float,
    smoothing_radius: float,
    quantile_threshold: float,
    risk_level: float,
    dpi: int,
    gdp_path: str | None,
    gdp_year: int,
    no_plots: bool,
    quiet: bool,
) -> None:
    """Run the LEC/EDC pipeline on CPC climate data and generate Cartopy maps.

    This reproduces the Justus (2025, Extremes) methodology on real
    NOAA CPC daily data (precipitation or tmax), adding VaR/ES risk
    metrics as an extension.  Use ``--variable precip`` for precipitation
    or ``--variable tmax`` for daily maximum temperature.

    \b
    Example:
        weatherisk maps --variable precip
        weatherisk maps --variable tmax --file-prefix tmax
        weatherisk maps --year-start 2010 --year-end 2020 --coarsen 2
    """
    from weatherisk.cpc_pipeline import PipelineConfig, run_cpc_pipeline, generate_maps as gen

    cfg = PipelineConfig(
        data_dir=data_dir,
        output_dir=output_dir,
        variable=variable,
        file_prefix=file_prefix if file_prefix is not None else variable,
        year_start=year_start,
        year_end=year_end,
        lat_range=lat_range,
        lon_range=lon_range,
        coarsen=coarsen,
        df=degrees_of_freedom,
        alpha=alpha,
        neighbor_radius=neighbor_radius,
        smoothing_radius=smoothing_radius,
        quantile_threshold=quantile_threshold,
        risk_level=risk_level,
        dpi=dpi,
        gdp_path=gdp_path,
        gdp_year=gdp_year,
    )

    verbose = not quiet
    result = run_cpc_pipeline(cfg, verbose=verbose)

    if no_plots:
        click.echo(f"Pipeline done. LEC k={result['k_lec']}, "
                    f"EDC k={result['k_edc']}. Plots skipped.")
        return

    saved = gen(result, verbose=verbose)
    click.echo(f"\nDone. {len(saved)} maps saved to {output_dir}")

"""Command-line interface entry point.

Provides the `weatherisk` CLI with subcommands:
- validate: run the method validation flow on synthetic data
- maps: run the CPC real-data pipeline
- cmip6: run the CMIP6 Figure 9 reproduction pipeline
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


@main.command()
@click.option(
    "--data-dir", default=None,
    help="Directory containing AWI-ESM-1-1-LR pr NetCDF files.",
)
@click.option("--output-dir", default="output/cmip6_fig9", help="Output directory.")
@click.option("--workers", default=1, help="Number of parallel workers.")
@click.option("--year-start", default=1850, help="First year (inclusive).")
@click.option("--year-end", default=2005, help="Last year (inclusive).")
@click.option("--dpi", default=300, help="Figure DPI.")
@click.option("--no-plots", is_flag=True, help="Skip figure generation.")
@click.option("--quiet", is_flag=True, help="Suppress progress output.")
def cmip6(
    data_dir: str | None,
    output_dir: str,
    workers: int,
    year_start: int,
    year_end: int,
    dpi: int,
    no_plots: bool,
    quiet: bool,
) -> None:
    """Reproduce Figure 9: LEC/EDC clustering on AWI-ESM-1-1-LR precipitation.

    Downloads data automatically if not found locally or on HPC.
    Uses STL de-trending, annual maxima of monthly data, GEV → Fréchet.

    \b
    Paper parameters: ν=5, α=1, ε=5 (grid point distance)
    Expected: k_EDC ≈ 104, k_LEC ≈ 24

    \b
    Example:
        weatherisk cmip6
        weatherisk cmip6 --workers 16
        weatherisk cmip6 --data-dir /pool/data/.../pr/gn/
    """
    from weatherisk.cmip6_pipeline import CMIP6Config, run_cmip6_pipeline, plot_figure9
    from weatherisk.cmip6_data import DEFAULT_DATA_DIR

    cfg = CMIP6Config(
        data_dir=data_dir or DEFAULT_DATA_DIR,
        output_dir=output_dir,
        year_start=year_start,
        year_end=year_end,
        n_workers=workers,
    )

    verbose = not quiet
    result = run_cmip6_pipeline(cfg, verbose=verbose)

    if no_plots:
        click.echo(f"Pipeline done. LEC k={result['k_lec']}, "
                    f"EDC k={result['k_edc']}. Plots skipped.")
        return

    saved = plot_figure9(result, dpi=dpi, verbose=verbose)
    click.echo(f"\nDone. {len(saved)} figures saved to {output_dir}")

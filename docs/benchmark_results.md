# Benchmark Results

This file records hot-path benchmark runs for the CMIP6 Figure 9 pipeline.
Each run uses the real pipeline functions on deterministic synthetic data.

## Run 2026-03-09T15:23:18.957400+00:00

- Git revision: `0bc6035`
- Total time: `1.590s`
- Config: `{'seed': 12345, 'n_years': 24, 'n_lat': 8, 'n_lon': 8, 'n_workers': 4, 'df': 5.0, 'alpha': 1.0, 'neighbor_radius': 3.0, 'smoothing_radius': 2.0, 'mle_ensemble': 3, 'stl_period': 12}`
- Derived: `{'n_months': 288, 'n_cells': 64, 'n_valid_cells': 64, 'n_years_complete': 24}`

| Step | Seconds |
| --- | ---: |
| step1a_detrend | 0.000 |
| step1b_annual_max | 0.126 |
| step2_gev_frechet | 0.692 |
| step3_local_mle | 0.770 |
| step5b_edc_matrix | 0.002 |

- Checks: `{'frechet_min': 0.17523249288252296, 'frechet_max': 193.8100447345456, 'edc_trace': 0.0, 'est_mean_a': 0.023881369901431172, 'est_mean_b': 1.283395966261962, 'est_mean_gamma': 0.36477793020266147}`

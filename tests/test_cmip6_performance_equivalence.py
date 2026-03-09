import numpy as np
import pytest
from scipy.stats import t as t_dist

from weatherisk.cmip6_pipeline import (
    CMIP6Config,
    _edc_condensed_flat,
    _compute_frechet_global,
    _edc_matrix_flat,
    _grid_coords,
    _incluster_reestimate_cmip6,
)
from weatherisk.density import pairwise_density_summand


def _reference_pairwise_density_summand(z1, z2, x, y, df, alpha, a, b, g):
    from weatherisk.covariance import cov_fkt_2d
    from weatherisk.density import _dtdiff

    cv = cov_fkt_2d(x, y, alpha, a, b, g)
    c = np.sqrt(1 - cv * cv) / np.sqrt(df + 1)

    m1 = ((z2 / z1) ** (1.0 / df) - cv) / c
    m2 = ((z1 / z2) ** (1.0 / df) - cv) / c

    dt_m1 = t_dist.pdf(m1, df + 1)
    dt_m2 = t_dist.pdf(m2, df + 1)
    pt_m1 = t_dist.cdf(m1, df + 1)
    pt_m2 = t_dist.cdf(m2, df + 1)

    term1_a = -pt_m1 / (z1 * z1)
    term1_b = -dt_m1 * z2 ** (1.0 / df) * z1 ** (-1.0 / df - 2) / c / df
    term1_c = dt_m2 * z1 ** (1.0 / df - 1) * z2 ** (-1.0 / df - 1) / c / df
    factor1 = term1_a + term1_b + term1_c

    term2_a = -pt_m2 / (z2 * z2)
    term2_b = -dt_m2 * z1 ** (1.0 / df) * z2 ** (-1.0 / df - 2) / c / df
    term2_c = dt_m1 * z2 ** (1.0 / df - 1) * z1 ** (-1.0 / df - 1) / c / df
    factor2 = term2_a + term2_b + term2_c

    dtd_m1 = _dtdiff(m1, df + 1)
    dtd_m2 = _dtdiff(m2, df + 1)
    cross = (
        dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1)
        + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1)
        + dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1) / df
        + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1) / df
        + dtd_m1 * z1 ** (-2.0 / df - 2) * z2 ** (2.0 / df - 1) / c / df
        + dtd_m2 * z2 ** (-2.0 / df - 2) * z1 ** (2.0 / df - 1) / c / df
    ) / c / df

    v = pt_m1 / z1 + pt_m2 / z2
    return np.log(np.maximum(factor1 * factor2 + cross, 1e-300)) - v


def _reference_edc_matrix_flat(frechet: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    n_years, n_cells = frechet.shape
    ranks = np.empty((n_cells, n_years))
    for s in range(n_cells):
        ranks[s] = rankdata(frechet[:, s])

    ec = np.zeros((n_cells, n_cells))
    for i in range(n_cells - 1):
        diff = np.abs(ranks[i] - ranks[i + 1:])
        v = diff.mean(axis=1) / (2.0 * (n_years + 1))
        denom = 1.0 - 2.0 * v
        denom[denom <= 0] = 1e-12
        ec[i, i + 1:] = np.minimum(1.0, (1.0 + 2.0 * v) / denom - 1.0)
    return ec + ec.T


def test_pairwise_density_matches_stats_reference_scalar_and_vector():
    scalar = pairwise_density_summand(1.5, 2.0, 1.0, 0.5, 5.0, 1.0, 2.0, 0.5, 0.2)
    scalar_ref = _reference_pairwise_density_summand(1.5, 2.0, 1.0, 0.5, 5.0, 1.0, 2.0, 0.5, 0.2)
    assert scalar == pytest.approx(scalar_ref, rel=1e-12, abs=1e-12)

    z1 = np.array([1.1, 1.5, 2.0, 3.0])
    z2 = np.array([1.3, 2.1, 1.8, 2.7])
    x = np.array([0.2, 0.5, 1.0, 1.5])
    y = np.array([0.1, -0.2, 0.7, -1.1])
    vec = pairwise_density_summand(z1, z2, x, y, 5.0, 1.0, 2.0, 0.5, -0.3)
    vec_ref = _reference_pairwise_density_summand(z1, z2, x, y, 5.0, 1.0, 2.0, 0.5, -0.3)
    np.testing.assert_allclose(vec, vec_ref, rtol=1e-12, atol=1e-12)


def test_compute_frechet_parallel_matches_serial():
    rng = np.random.default_rng(123)
    annual_max = rng.gamma(shape=2.5, scale=1.2, size=(18, 4, 5))
    annual_max[:, 0, 0] = annual_max[:, 0, 1]

    fr_serial, idx_serial = _compute_frechet_global(annual_max, n_workers=1, verbose=False)
    fr_parallel, idx_parallel = _compute_frechet_global(annual_max, n_workers=2, verbose=False)

    np.testing.assert_array_equal(idx_serial, idx_parallel)
    np.testing.assert_allclose(fr_serial, fr_parallel, rtol=1e-10, atol=1e-10)


def test_edc_matrix_fast_matches_reference_loop():
    rng = np.random.default_rng(456)
    frechet = rng.lognormal(mean=0.3, sigma=0.5, size=(20, 9))
    fast = _edc_matrix_flat(frechet)
    ref = _reference_edc_matrix_flat(frechet)
    np.testing.assert_allclose(fast, ref, rtol=0, atol=1e-12)


def test_edc_condensed_matches_squareform_of_full_matrix():
    from scipy.spatial.distance import squareform

    rng = np.random.default_rng(654)
    frechet = rng.lognormal(mean=0.1, sigma=0.4, size=(18, 10))
    full = _edc_matrix_flat(frechet)
    condensed = _edc_condensed_flat(frechet)
    np.testing.assert_allclose(condensed, squareform(full), rtol=0, atol=1e-12)


def test_incluster_reestimate_parallel_matches_serial():
    rng = np.random.default_rng(789)
    frechet = rng.lognormal(mean=0.2, sigma=0.35, size=(18, 12))
    valid_idx = np.arange(12)
    grid_coords = _grid_coords(valid_idx, 3, 4)
    labels = np.array([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])

    cfg_serial = CMIP6Config(
        neighbor_radius=3.0,
        mle_ensemble=3,
        n_workers=1,
    )
    cfg_parallel = CMIP6Config(
        neighbor_radius=3.0,
        mle_ensemble=3,
        n_workers=2,
    )

    serial = _incluster_reestimate_cmip6(
        frechet, grid_coords, labels, cfg_serial, "TEST", verbose=False
    )
    parallel = _incluster_reestimate_cmip6(
        frechet, grid_coords, labels, cfg_parallel, "TEST", verbose=False
    )

    assert set(serial) == set(parallel)
    for cl in serial:
        np.testing.assert_allclose(serial[cl], parallel[cl], rtol=1e-10, atol=1e-10)
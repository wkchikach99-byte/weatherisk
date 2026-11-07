"""weatherisk.stl against frozen outputs of R's stats::stl.

R is never called. The references come from

* ``tests/reference_data/stl/`` (``tests/generate_stl_reference.R``), stored
  as raw float64 so no decimal rounding enters the comparison, and
* ``tests/reference_data/cmip6_mini/stl_trends.csv``
  (``tests/generate_cmip6_mini_reference.R``), the trends the mini Figure 9
  chain detrends with.

Tolerances. The port and R agree to rounding, not bit for bit: R's
compiled loess sums use fused multiply-adds and the port does not, and
fifteen robustness passes amplify those differences. Every array is
held to a relative error of 1e-12. The trend is compared elementwise.
Seasonal and remainder are compared against the series' largest
magnitude, and weights against their scale of 1, because all three have
values near zero (a bisquare weight near zero is 1 - u^2 squared with u
near 1, so its elementwise relative error is cancellation, not
disagreement). ``dry240`` and the CSV fixture need their own bounds,
explained at their tests.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import weatherisk.stl as stl_mod
from weatherisk.stl import r_stl_params, stl

REF = Path(__file__).parent / "reference_data"
STL_REF = REF / "stl"
RTOL = 1e-12

CASES = pd.read_csv(STL_REF / "cases.csv")
WELL_CONDITIONED = CASES[CASES.case != "dry240"]


def _load(case: pd.Series):
    n, m = int(case.n), int(case.m)
    raw = np.fromfile(STL_REF / f"{case.case}.f64", dtype="<f8")
    assert raw.size == 5 * n * m
    x, seasonal, trend, remainder, weights = (
        raw[k * n * m:(k + 1) * n * m].reshape(m, n).T for k in range(5))
    return x, seasonal, trend, remainder, weights


def _kwargs(case: pd.Series) -> dict:
    kw = {"robust": bool(case.robust)}
    if not pd.isna(case.s_degree_arg):
        kw["s_degree"] = int(case.s_degree_arg)
    if case.jump1:
        kw.update(s_jump=1, t_jump=1, l_jump=1)
    if not pd.isna(case.inner_arg):
        kw["inner"] = int(case.inner_arg)
    if not pd.isna(case.outer_arg):
        kw["outer"] = int(case.outer_arg)
    return kw


def _s_window(case: pd.Series):
    return "periodic" if case.s_window == "periodic" else int(case.s_window)


def _fit(case: pd.Series, x: np.ndarray):
    return stl(x, int(case.period), _s_window(case), **_kwargs(case))


def _rel_to_series_max(a, b, x):
    return np.abs(a - b).max(axis=0) / np.abs(x).max(axis=0)


@pytest.mark.parametrize("case", [c for _, c in CASES.iterrows()],
                         ids=list(CASES.case))
def test_parameters_resolve_like_r(case):
    p = r_stl_params(int(case.n), int(case.period), _s_window(case),
                     **_kwargs(case))
    got = (p.s_window, p.t_window, p.l_window, p.s_degree, p.t_degree,
           p.l_degree, p.s_jump, p.t_jump, p.l_jump, p.inner, p.outer)
    want = tuple(int(case[k]) for k in (
        "win_s", "win_t", "win_l", "deg_s", "deg_t", "deg_l",
        "jump_s", "jump_t", "jump_l", "inner", "outer"))
    assert got == want


@pytest.mark.parametrize("case", [c for _, c in WELL_CONDITIONED.iterrows()],
                         ids=list(WELL_CONDITIONED.case))
def test_components_match_r(case):
    x, seasonal, trend, remainder, weights = _load(case)
    f = _fit(case, x)
    np.testing.assert_allclose(f.trend, trend, rtol=RTOL, atol=0)
    assert np.abs(f.weights - weights).max() <= RTOL
    assert _rel_to_series_max(f.seasonal, seasonal, x).max() <= RTOL
    assert _rel_to_series_max(f.remainder, remainder, x).max() <= RTOL


def _case(name):
    return CASES[CASES.case == name].iloc[0]


def test_dry_series_reach_zero_weights_like_r(monkeypatch):
    """dry240 is 53-75% exact zeros. It reaches the Fortran ok = FALSE
    branch (a window whose robustness weights are all zero), which no
    other case does.

    It is not held to 1e-12 because R is not. Once more than half the
    residuals are near zero, the median absolute residual collapses (to
    3e-10 in one series) and the weights amplify rounding. Perturbing the
    nonzero inputs by one ulp moves R 4.5.2's own trends by 7.0e-10
    relative and its weights by 9.7e-6; the port differs from R by
    3.9e-10 and 9.5e-6 (all measured). The test therefore checks the
    discrete outcome exactly and the values at R's own sensitivity.
    """
    calls = {"not_ok": 0}
    est = stl_mod._est

    def counting_est(*args):
        ys, ok = est(*args)
        calls["not_ok"] += int((~ok).sum())
        return ys, ok

    monkeypatch.setattr(stl_mod, "_est", counting_est)
    case = _case("dry240")
    x, _, trend, _, weights = _load(case)
    f = _fit(case, x)
    assert calls["not_ok"] > 0
    assert np.array_equal(f.weights == 0, weights == 0)
    np.testing.assert_allclose(f.trend, trend, rtol=1e-9, atol=0)
    assert np.abs(f.weights - weights).max() <= 2e-5


def test_psort_differs_from_the_exact_median(monkeypatch):
    """R's rwts calls psort with descending indices, which for even n can
    return a neighbouring order statistic instead of the median. The
    robust240 fit makes such passes, so an exact median would not match."""
    calls = {"off": 0}
    psort = stl_mod._psort

    def checking_psort(values, ind):
        got = psort(values, ind)
        s = sorted(values)
        calls["off"] += got != [s[i - 1] for i in ind]
        return got

    monkeypatch.setattr(stl_mod, "_psort", checking_psort)
    case = _case("robust240")
    _fit(case, _load(case)[0])
    assert calls["off"] > 0


def test_psort_matches_full_sort_for_ascending_indices():
    rng = np.random.default_rng(0)
    for n in (2, 3, 11, 12, 240, 1871):
        v = rng.normal(size=n)
        v[: n // 4] = 0.0  # ties
        ind = tuple(sorted({1, (n + 1) // 2, n}))
        s = np.sort(v)
        assert stl_mod._psort(v.tolist(), ind) == [s[i - 1] for i in ind]


def test_single_series_equals_batch_column():
    case = _case("robust240")
    x = _load(case)[0]
    batch = stl(x, 12, robust=True)
    for k in (0, 7, 24):
        one = stl(x[:, k], 12, robust=True)
        for a, b in zip(one, batch):
            assert np.array_equal(a, b[:, k])


def test_rejects_non_finite_input():
    x = np.ones(240)
    x[3] = np.nan
    with pytest.raises(ValueError):
        stl(x, 12)


def test_mini_fixture_trends_match_r():
    """Robust periodic STL of the 25 mini cells against stl_trends.csv.

    write.csv rounded both the input and R's trends to 15 significant
    digits. Rerunning R 4.5.2 on the rounded input moves its trends by up
    to 5.7e-14 relative to the CSV (measured), so 1e-13 is the floor this
    fixture supports. The raw-float64 cases above carry the 1e-12 target.
    """
    mini = REF / "cmip6_mini"
    x = pd.read_csv(mini / "monthly_pr.csv",
                    float_precision="round_trip").to_numpy()
    want = pd.read_csv(mini / "stl_trends.csv",
                       float_precision="round_trip").to_numpy()
    got = stl(x, 12, robust=True).trend
    np.testing.assert_allclose(got, want, rtol=1e-13, atol=0)

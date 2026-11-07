"""Seasonal-trend decomposition by Loess (STL), as computed by R's ``stl()``.

A NumPy port of the Fortran STL of Cleveland, Cleveland, McRae and
Terpenning (1990), "STL: A Seasonal-Trend Decomposition Procedure Based
on Loess", Journal of Official Statistics 6(1), 3-73.

Sources followed
----------------
* The algorithm was ported from the original netlib code
  (https://netlib.org/a/stl, file ``stl.f``: subroutines ``stl``,
  ``onestp``, ``ss``, ``fts``, ``ma``, ``ess``, ``est``, ``rwts`` and
  ``psort``).
  That file carries no authorship or copyright notice.
* R's ``src/library/stats/src/stl.f`` was read only to confirm that R
  computes the same algorithm. R's changes are double precision
  (B. D. Ripley, 1999) and restructured control flow (M. Maechler, 2001),
  neither of which changes the arithmetic. This port computes in float64,
  like R.
* The parameter defaults and the ``s.window = "periodic"`` post-step
  follow R 4.5.2's ``stats::stl``, as printed by ``body(stats::stl)``.
  No R code was copied.

Floating point
--------------
Every sum runs sequentially in the Fortran loop order (``np.cumsum``
along the time axis). The port does not fuse multiply-adds. R's compiled
``stats.so`` on arm64 macOS does (``fmadd``/``fmsub`` in ``stlest_`` and
``stless_``), so results agree with R to rounding, not bit for bit.

Arrays
------
A batch of ``m`` series is a float64 array of shape ``(n, m)``, time
along axis 0. Every time-step operation runs over all ``m`` cells at
once. Robustness weights differ per cell, so the Fortran ``ok`` branches
(a window whose weights sum to zero) become per-cell masks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

__all__ = ["STLParams", "STLResult", "r_stl_params", "stl"]


@dataclass(frozen=True)
class STLParams:
    """The arguments R's ``stl()`` passes to the Fortran routine."""

    period: int
    s_window: int
    s_degree: int
    t_window: int
    t_degree: int
    l_window: int
    l_degree: int
    s_jump: int
    t_jump: int
    l_jump: int
    inner: int
    outer: int
    periodic: bool


class STLResult(NamedTuple):
    """Components with the shape of the input series."""

    seasonal: np.ndarray
    trend: np.ndarray
    remainder: np.ndarray
    weights: np.ndarray


def _nextodd(x: float) -> int:
    x = round(x)
    return int(x + 1 if x % 2 == 0 else x)


def r_stl_params(
    n: int,
    period: int,
    s_window: int | str = "periodic",
    *,
    s_degree: int = 0,
    t_window: int | None = None,
    t_degree: int = 1,
    l_window: int | None = None,
    l_degree: int | None = None,
    s_jump: int | None = None,
    t_jump: int | None = None,
    l_jump: int | None = None,
    robust: bool = False,
    inner: int | None = None,
    outer: int | None = None,
) -> STLParams:
    """Resolve parameters exactly as R 4.5.2's ``stl()`` defaults them."""
    if period < 2 or n <= 2 * period:
        raise ValueError("series is not periodic or has less than two periods")
    periodic = s_window == "periodic"
    if periodic:
        s_window = 10 * n + 1
        s_degree = 0
    elif isinstance(s_window, str):
        raise ValueError("unknown string value for s_window")
    for name, deg in (("s_degree", s_degree), ("t_degree", t_degree),
                      ("l_degree", l_degree)):
        if deg is not None and deg not in (0, 1):
            raise ValueError(f"{name} must be 0 or 1")
    if t_window is None:
        t_window = _nextodd(math.ceil(1.5 * period / (1 - 1.5 / s_window)))
    if l_window is None:
        l_window = _nextodd(period)
    return STLParams(
        period=period,
        s_window=int(s_window),
        s_degree=s_degree,
        t_window=int(t_window),
        t_degree=t_degree,
        l_window=int(l_window),
        l_degree=t_degree if l_degree is None else l_degree,
        s_jump=math.ceil(s_window / 10) if s_jump is None else s_jump,
        t_jump=math.ceil(t_window / 10) if t_jump is None else t_jump,
        l_jump=math.ceil(l_window / 10) if l_jump is None else l_jump,
        inner=(1 if robust else 2) if inner is None else inner,
        outer=(15 if robust else 0) if outer is None else outer,
        periodic=periodic,
    )


def stl(
    x: np.ndarray,
    period: int,
    s_window: int | str = "periodic",
    *,
    robust: bool = False,
    **kwargs,
) -> STLResult:
    """STL decomposition of one series ``(n,)`` or a batch ``(n, m)``.

    Equivalent to ``stl(ts(x, frequency = period), s.window, robust = robust,
    ...)`` in R, column by column. The first observation is cycle position 1.
    Keyword arguments are the remaining R arguments in snake case
    (``t_window``, ``t_jump``, ``inner``, ...); see :func:`r_stl_params`.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim not in (1, 2):
        raise ValueError("x must have shape (n,) or (n, m)")
    if not np.isfinite(x).all():
        raise ValueError("x must be finite (R's stl uses na.fail)")
    y = np.ascontiguousarray(x.reshape(x.shape[0], -1))
    p = r_stl_params(y.shape[0], period, s_window, robust=robust, **kwargs)

    seasonal, trend, weights = _stl_fortran(y, p)
    if p.periodic:
        seasonal = _cycle_means(seasonal, p.period)
    remainder = y - seasonal - trend
    return STLResult(*(a.reshape(x.shape)
                       for a in (seasonal, trend, remainder, weights)))


def _seqsum(a: np.ndarray) -> np.ndarray:
    # np.sum switches to pairwise summation along a contiguous axis; a
    # cumulative sum is sequential by construction, like the Fortran loops.
    return np.cumsum(a, axis=0)[-1]


def _cycle_means(seasonal: np.ndarray, period: int) -> np.ndarray:
    """R's ``tapply(seasonal, cycle(x), mean)[cycle(x)]``.

    R's ``mean`` adds a second pass that sums the deviations from the
    first-pass mean. Its accumulator is a C ``long double``, which is
    float64 on arm64 macOS but 80-bit on x86.
    """
    out = np.empty_like(seasonal)
    for c in range(period):
        s = seasonal[c::period]
        k = s.shape[0]
        mean = _seqsum(s) / k
        mean = mean + _seqsum(s - mean) / k
        out[c::period] = mean
    return out


def _stl_fortran(y: np.ndarray, p: STLParams):
    """Fortran ``stl``: inner loop, then ``outer`` robustness passes."""
    n, m = y.shape
    ns, nt, nl = (w + 1 if w % 2 == 0 else w
                  for w in (max(3, p.s_window), max(3, p.t_window),
                            max(3, p.l_window)))
    np_ = max(2, p.period)
    season = np.zeros((n, m))
    trend = np.zeros((n, m))
    rw = None
    k = 0
    while True:
        season, trend = _onestp(y, np_, ns, nt, nl, p, rw, season, trend)
        k += 1
        if k > p.outer:
            break
        rw = _rwts(y, trend + season)
    weights = np.ones((n, m)) if p.outer <= 0 else rw
    return season, trend, weights


def _onestp(y, np_, ns, nt, nl, p, rw, season, trend):
    n = y.shape[0]
    for _ in range(p.inner):
        cycle = _ss(y - trend, np_, ns, p.s_degree, p.s_jump, rw)
        low = _ma(_ma(_ma(cycle, np_), np_), 3)
        low = _ess(low, nl, p.l_degree, p.l_jump, None)
        season = cycle[np_:np_ + n] - low
        trend = _ess(y - season, nt, p.t_degree, p.t_jump, rw)
    return season, trend


def _ss(y, np_, ns, isdeg, nsjump, rw):
    """Smooth each cycle-subseries and extend it by one point at each end.

    Returns shape ``(n + 2*np_, m)``.
    """
    n, m = y.shape
    season = np.empty((n + 2 * np_, m))
    for j in range(np_):
        sub = np.ascontiguousarray(y[j::np_])
        subw = None if rw is None else np.ascontiguousarray(rw[j::np_])
        k = sub.shape[0]
        work = np.empty((k + 2, m))
        work[1:k + 1] = _ess(sub, ns, isdeg, nsjump, subw)
        ys, ok = _est(sub, ns, isdeg, 0.0, 1, min(ns, k), subw)
        work[0] = np.where(ok, ys, work[1])
        ys, ok = _est(sub, ns, isdeg, float(k + 1), max(1, k - ns + 1), k,
                      subw)
        work[k + 1] = np.where(ok, ys, work[k])
        season[j::np_] = work
    return season


def _ma(x, length):
    """Running mean of ``length`` points, updated as Fortran ``ma`` does."""
    newn = x.shape[0] - length + 1
    ave = np.empty((newn, x.shape[1]))
    flen = float(length)
    v = _seqsum(x[:length])
    ave[0] = v / flen
    for j in range(1, newn):
        v = v - x[j - 1] + x[j - 1 + length]
        ave[j] = v / flen
    return ave


def _ess(y, length, ideg, njump, rw):
    """Loess smoothing of every column of ``y``, evaluated every ``njump``
    points and linearly interpolated in between."""
    n = y.shape[0]
    ys = np.empty_like(y)
    if n < 2:
        ys[0] = y[0]
        return ys

    def fit(i, nleft, nright):
        v, ok = _est(y, length, ideg, float(i), nleft, nright, rw)
        ys[i - 1] = np.where(ok, v, y[i - 1])

    newnj = min(njump, n - 1)
    nleft = nright = 0
    if length >= n:
        nleft, nright = 1, n
        for i in range(1, n + 1, newnj):
            fit(i, nleft, nright)
    elif newnj == 1:
        nsh = (length + 1) // 2
        nleft, nright = 1, length
        for i in range(1, n + 1):
            if i > nsh and nright != n:
                nleft += 1
                nright += 1
            fit(i, nleft, nright)
    else:
        nsh = (length + 1) // 2
        for i in range(1, n + 1, newnj):
            if i < nsh:
                nleft, nright = 1, length
            elif i >= n - nsh + 1:
                nleft, nright = n - length + 1, n
            else:
                nleft, nright = i - nsh + 1, length + i - nsh
            fit(i, nleft, nright)

    if newnj != 1:
        steps = np.arange(1, newnj, dtype=np.float64)[:, None]
        for i in range(1, n - newnj + 1, newnj):
            delta = (ys[i + newnj - 1] - ys[i - 1]) / float(newnj)
            ys[i:i + newnj - 1] = ys[i - 1] + delta * steps
        k = ((n - 1) // newnj) * newnj + 1
        if k != n:
            # Fortran reuses the window of the last evaluated point here.
            fit(n, nleft, nright)
            if k != n - 1:
                delta = (ys[n - 1] - ys[k - 1]) / float(n - k)
                tail = np.arange(1, n - k, dtype=np.float64)[:, None]
                ys[k:n - 1] = ys[k - 1] + delta * tail
    return ys


def _est(y, length, ideg, xs, nleft, nright, rw):
    """Local fit at abscissa ``xs`` from rows ``nleft..nright`` (1-based).

    Returns the fitted value per column and the per-column ``ok`` flag,
    False where the weights in the window sum to zero.
    """
    n = y.shape[0]
    rng = float(n) - 1.0
    h = max(xs - nleft, nright - xs)
    if length > n:
        h += (length - n) // 2
    j = np.arange(nleft, nright + 1, dtype=np.float64)
    r = np.abs(j - xs)
    w = np.zeros_like(r)
    near = r <= 0.001 * h
    mid = (r <= 0.999 * h) & ~near
    q = r[mid] / h
    t = 1.0 - q * (q * q)
    w[mid] = (t * t) * t
    w[near] = 1.0
    w = w[:, None]
    if rw is not None:
        w = rw[nleft - 1:nright] * w

    a = _seqsum(w)
    ok = a > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        w = w / a
        if h > 0.0 and ideg > 0:
            jj = j[:, None]
            a = _seqsum(w * jj)
            b = xs - a
            d = jj - a
            c = _seqsum(w * (d * d))
            steep = np.sqrt(c) > 0.001 * rng
            b = b / c
            w = np.where(steep, w * (b * d + 1.0), w)
        ys = _seqsum(w * y[nleft - 1:nright])
    return ys, ok


def _rwts(y, fit):
    """Tukey bisquare robustness weights with scale ``cmad``.

    ``cmad`` is meant to be 6 * median |residual|, but Fortran ``rwts``
    passes ``psort`` the indices ``(n/2+1, n-n/2)`` in descending order.
    For even ``n``, ``psort`` then sometimes leaves a neighbouring order
    statistic at position ``n/2``, and R's weights use that value. The
    exact median would not reproduce R, so ``psort`` is emulated per cell.
    """
    n = y.shape[0]
    r = np.abs(y - fit)
    mid1 = n // 2 + 1
    mid2 = n - mid1 + 1
    cmad = np.array([3.0 * sum(_psort(col, (mid1, mid2)))
                     for col in r.T.tolist()])
    with np.errstate(divide="ignore", invalid="ignore"):
        u = r / cmad
        t = 1.0 - u * u
        rw = np.where(r <= 0.999 * cmad, t * t, 0.0)
    return np.where(r <= 0.001 * cmad, 1.0, rw)


def _psort(values: list[float], ind: tuple[int, ...]) -> list[float]:
    """Fortran ``psort``: partial quicksort meant to put ``a(ind(.))`` in
    sorted position. Returns ``a(ind(.))`` after the sort.

    A line-by-line port, 1-based, of the netlib routine: a quicksort with
    median-of-three pivots, a 16-level segment stack, insertion sort for
    segments of at most 11 points, and bookkeeping that skips segments
    holding no requested index. The bookkeeping assumes ``ind`` is
    ascending; it is reproduced as written, not corrected.
    """
    a = [0.0, *values]
    n, ni = len(values), len(ind)
    idx = (0, *ind)
    if n < 2 or ni == 0:
        return [a[i] for i in ind]
    il, iu, indl, indu = [0] * 17, [0] * 17, [0] * 17, [0] * 17
    jl, ju = 1, ni
    indl[1], indu[1] = 1, ni
    i, j, m = 1, n, 1
    label = 161
    while True:
        if label == 161:
            label = 10 if i < j else 166
        elif label == 166:
            m -= 1
            if m == 0:
                return [a[i] for i in ind]
            i, j, jl, ju = il[m], iu[m], indl[m], indu[m]
            label = 173 if jl <= ju else 166
        elif label == 173:
            label = 10 if j - i > 10 else 174
        elif label == 10:
            k = i
            ij = (i + j) // 2
            t = a[ij]
            if a[i] > t:
                a[ij], a[i] = a[i], t
                t = a[ij]
            lo = j
            if a[j] < t:
                a[ij], a[j] = a[j], t
                t = a[ij]
                if a[i] > t:
                    a[ij], a[i] = a[i], t
                    t = a[ij]
            while True:
                lo -= 1
                if a[lo] <= t:
                    tt = a[lo]
                    k += 1
                    while not a[k] >= t:
                        k += 1
                    if k > lo:
                        break
                    a[lo] = a[k]
                    a[k] = tt
            indl[m], indu[m] = jl, ju
            p = m
            m += 1
            label = 173
            if lo - i <= j - k:
                il[p], iu[p] = k, j
                j = lo
                while True:
                    if jl > ju:
                        label = 166
                        break
                    if idx[ju] > j:
                        ju -= 1
                    else:
                        indl[p] = ju + 1
                        break
            else:
                il[p], iu[p] = i, lo
                i = k
                while True:
                    if jl > ju:
                        label = 166
                        break
                    if idx[jl] < i:
                        jl += 1
                    else:
                        indu[p] = jl - 1
                        break
        else:  # 174: insertion sort, which needs a(i-1) as a sentinel
            if i == 1:
                label = 161
                continue
            i -= 1
            while True:
                i += 1
                if i == j:
                    break
                t = a[i + 1]
                if a[i] > t:
                    k = i
                    while True:
                        a[k + 1] = a[k]
                        k -= 1
                        if t >= a[k]:
                            break
                    a[k + 1] = t
            label = 166

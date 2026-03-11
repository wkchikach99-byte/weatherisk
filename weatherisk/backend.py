"""Backend dispatch: Rust (weatherisk_core) or pure-Python fallback.

Import this module wherever the hot kernels are needed.  It exposes a
uniform API regardless of whether the Rust extension is installed:

    from weatherisk.backend import neg_log_likelihood_sum, calc_distance_ellipses

The active backend can be forced via the ``WEATHERISK_BACKEND`` env var:

    WEATHERISK_BACKEND=python   — always use pure-Python
    WEATHERISK_BACKEND=rust     — require Rust (error if missing)
    (unset)                     — auto-detect (Rust if available)
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

# ── Resolve backend ──────────────────────────────────────────────────────

_BACKEND_ENV = os.environ.get("WEATHERISK_BACKEND", "").lower()

_USE_RUST = False

if _BACKEND_ENV == "python":
    _USE_RUST = False
elif _BACKEND_ENV == "rust":
    try:
        import weatherisk_core as _rc  # noqa: F401

        _USE_RUST = True
    except ImportError as exc:
        raise ImportError(
            "WEATHERISK_BACKEND=rust but weatherisk_core is not installed. "
            "Build with: cd crates/weatherisk_core && maturin develop --release"
        ) from exc
else:
    # Auto-detect
    try:
        import weatherisk_core as _rc  # noqa: F401

        _USE_RUST = True
    except ImportError:
        _USE_RUST = False

BACKEND: str = "rust" if _USE_RUST else "python"

_REPO_ROOT = Path(__file__).resolve().parents[1]
_R_LOCAL_MLE_CACHE: dict[str, np.ndarray] = {}


def _rscript_path() -> str:
        candidate = shutil.which("Rscript")
        if candidate is None:
                raise RuntimeError(
                        "Rscript not found on PATH. "
                        "R is required for exact numerical parity."
                )
        return candidate


def _local_mle_cache_key(
        zi: np.ndarray,
        zj: np.ndarray,
        xl: np.ndarray,
        yl: np.ndarray,
        df: float,
        alpha: float,
        ensemble: int,
        seed: int,
        max_boundary_retries: int,
    lower_bounds: tuple[float, float],
    upper_bounds: tuple[float, float],
) -> str:
        hasher = hashlib.sha256()
        for arr in (zi, zj, xl, yl):
                arr64 = np.ascontiguousarray(arr, dtype=np.float64)
                hasher.update(arr64.tobytes())
        hasher.update(
            (
                f"{df:.17g}|{alpha:.17g}|{ensemble}|{seed}|{max_boundary_retries}|"
                f"{lower_bounds[0]:.17g}|{lower_bounds[1]:.17g}|"
                f"{upper_bounds[0]:.17g}|{upper_bounds[1]:.17g}"
            ).encode()
        )
        return hasher.hexdigest()


def _optimize_local_mle_via_r(
        zi: np.ndarray,
        zj: np.ndarray,
        xl: np.ndarray,
        yl: np.ndarray,
        df: float,
        alpha: float,
        ensemble: int,
        seed: int,
        max_boundary_retries: int,
    lower_bounds: tuple[float, float],
    upper_bounds: tuple[float, float],
) -> np.ndarray:
        rscript = _rscript_path()

        cache_key = _local_mle_cache_key(
        zi, zj, xl, yl, df, alpha, ensemble, seed, max_boundary_retries,
        lower_bounds, upper_bounds,
        )
        cached = _R_LOCAL_MLE_CACHE.get(cache_key)
        if cached is not None:
                return cached.copy()

        arrays = np.column_stack([
                np.asarray(zi, dtype=np.float64),
                np.asarray(zj, dtype=np.float64),
                np.asarray(xl, dtype=np.float64),
                np.asarray(yl, dtype=np.float64),
        ])
        n_rows = arrays.shape[0]
        arr_c = np.ascontiguousarray(arrays, dtype=np.float64)

        r_code = (
            """
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == '--args') {
    args <- args[-1]
}
input_path <- args[1]
df_val <- as.numeric(args[2])
alpha_val <- as.numeric(args[3])
ensemble_val <- as.integer(args[4])
seed_val <- as.integer(args[5])
max_boundary_retries_val <- as.integer(args[6])
lower_a <- as.numeric(args[7])
lower_b <- as.numeric(args[8])
upper_a <- as.numeric(args[9])
upper_b <- as.numeric(args[10])

suppressMessages(library(lhs))
source(file.path(getwd(), 'r_code', 'functions.R'))

"""
            + f'raw <- readBin(input_path, what="double", n={n_rows * 4}, size=8, endian="little")\n'
            + f"dat <- matrix(raw, nrow={n_rows}, ncol=4, byrow=TRUE)\n"
            + """
zlist1 <- dat[,1]
zlist2 <- dat[,2]
Xlist <- dat[,3]
Ylist <- dat[,4]

llh <- function(par) {
    sum(pairwise_density_summand(
        zlist1, zlist2, Xlist, Ylist,
        df_val, alpha_val, par[1], par[2], par[3]
    ))
}

parameters_lower_bound <- c(lower_a, lower_b, -pi / 2)
parameters_upper_bound <- c(upper_a, upper_b, pi / 2)
parscale <- (parameters_upper_bound - parameters_lower_bound) / 100

set.seed(seed_val)
starts <- matrix(
    parameters_lower_bound,
    ensemble_val + max_boundary_retries_val,
    length(parameters_lower_bound),
    byrow = TRUE
) + maximinLHS(
    ensemble_val + max_boundary_retries_val,
    length(parameters_lower_bound)
) %*% diag(parameters_upper_bound - parameters_lower_bound)

num_calc <- 1L
num_more_calc_done <- 0L
fit <- NULL
while (num_calc <= ensemble_val) {
    fit_ <- optim(
        starts[num_calc + num_more_calc_done, ],
        fn = llh,
        method = 'L-BFGS-B',
        lower = parameters_lower_bound,
        upper = parameters_upper_bound,
        control = list(fnscale = -1, parscale = parscale, maxit = 10000)
    )
    if (abs(fit_$par[3]) == pi / 2) {
        fit_ <- optim(
            c(fit_$par[1], fit_$par[2], -fit_$par[3]),
            fn = llh,
            method = 'L-BFGS-B',
            lower = parameters_lower_bound,
            upper = parameters_upper_bound,
            control = list(fnscale = -1, parscale = parscale, maxit = 10000)
        )
    }
    if (is.null(fit) || (-fit_$value < -fit$value)) {
        fit <- fit_
    }
    if (
        min(abs(fit_$par - parameters_lower_bound)) < 0.01 ||
        min(abs(fit_$par - parameters_upper_bound)) < 0.01
    ) {
        if (num_more_calc_done < max_boundary_retries_val) {
            num_more_calc_done <- num_more_calc_done + 1L
            next
        }
    }
    num_calc <- num_calc + 1L
}

cat(sprintf('%.17g,%.17g,%.17g\n', fit$par[1], fit$par[2], fit$par[3]))
"""
        )

        try:
                with tempfile.TemporaryDirectory() as tmpdir:
                        input_path = Path(tmpdir) / "local_mle_inputs.bin"
                        arr_c.tofile(input_path)
                        result = subprocess.run(
                                [
                                        rscript,
                                        "-e",
                                        r_code,
                                "--args",
                                        str(input_path),
                                        f"{df:.17g}",
                                        f"{alpha:.17g}",
                                        str(int(ensemble)),
                                        str(int(seed)),
                                        str(int(max_boundary_retries)),
                                        f"{lower_bounds[0]:.17g}",
                                        f"{lower_bounds[1]:.17g}",
                                        f"{upper_bounds[0]:.17g}",
                                        f"{upper_bounds[1]:.17g}",
                                ],
                                check=True,
                                capture_output=True,
                                text=True,
                                cwd=_REPO_ROOT,
                        )
        except (OSError, subprocess.SubprocessError, ValueError) as exc:
                raise RuntimeError(f"R local MLE optimizer subprocess failed: {exc}") from exc

        output = result.stdout.strip().splitlines()
        if not output:
                raise RuntimeError("R local MLE optimizer returned no output")
        try:
                est = np.fromstring(output[-1], sep=",", dtype=np.float64)
        except ValueError as exc:
                raise RuntimeError(f"R local MLE optimizer returned unparseable output: {output[-1]}") from exc
        if est.shape != (3,):
                raise RuntimeError(f"R local MLE optimizer returned {est.shape[0]} params, expected 3")

        _R_LOCAL_MLE_CACHE[cache_key] = est.copy()
        return est


# ── Neg log-likelihood sum ───────────────────────────────────────────────


def neg_log_likelihood_sum(
    z1: np.ndarray,
    z2: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    df: float,
    alpha: float,
    a: float,
    b: float,
    g: float,
) -> float:
    """Compute -sum(pairwise_density_summand(...)) over paired arrays."""
    if _USE_RUST:
        return _rc.neg_log_likelihood_sum(z1, z2, x, y, df, alpha, a, b, g)
    else:
        from weatherisk.density import pairwise_density_summand

        return -float(
            np.sum(pairwise_density_summand(z1, z2, x, y, df, alpha, a, b, g))
        )


# ── LEC dissimilarity matrix ────────────────────────────────────────────


def calc_distance_ellipses(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Jaccard-like ellipse overlap dissimilarity matrix (n×n, 0–100)."""
    if _USE_RUST:
        return _rc.calc_distance_ellipses(estimates, res)
    else:
        from weatherisk.clustering import (
            calc_distance_ellipses as _py_calc_distance_ellipses,
        )

        return _py_calc_distance_ellipses(estimates, res=res, chunk_size=chunk_size)


def calc_distance_ellipses_condensed(
    estimates: np.ndarray,
    res: int = 21,
    chunk_size: int | None = None,
) -> np.ndarray:
    """Condensed upper-triangle ellipse dissimilarity (n*(n-1)/2, 0–100).

    Uses the Rust backend's native condensed function when available,
    avoiding the full n×n matrix allocation entirely.
    """
    if _USE_RUST:
        return _rc.calc_distance_ellipses_condensed(estimates, res)
    else:
        from weatherisk.clustering import (
            calc_distance_ellipses_condensed as _py_condensed,
        )

        return _py_condensed(estimates, res=res, chunk_size=chunk_size)


# ── Full optimizer loops ─────────────────────────────────────────────────


def _nll_and_grad_factory(zi, zj, xl, yl, df, alpha, parscale):
    """Return an (f, jac_flag) callable for SciPy minimize.

    Parameters operate in *scaled* space (par_scaled = par / parscale).
    Both Rust and Python paths return (value, gradient) so that
    ``jac=True`` is used uniformly — ensuring identical optimization
    trajectories regardless of backend.

    Parameters
    ----------
    parscale : ndarray, shape (3,)
        Scaling factors; real params = par_scaled * parscale.
    """
    if _USE_RUST:
        def f_and_grad(p_scaled):
            p = p_scaled * parscale
            fval, grad_raw = _rc.nll_with_gradient(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            return fval, grad_raw * parscale
        return f_and_grad, True
    else:
        # Pure-Python: forward-difference gradient (matches Rust semantics)
        _EPS_FD = 1e-5

        def f_and_grad(p_scaled):
            p = p_scaled * parscale
            f0 = neg_log_likelihood_sum(
                zi, zj, xl, yl, df, alpha, p[0], p[1], p[2],
            )
            if not np.isfinite(f0):
                f0 = 1e20
            grad = np.empty(3)
            for k in range(3):
                p_k = p.copy()
                p_k[k] += _EPS_FD
                fk = neg_log_likelihood_sum(
                    zi, zj, xl, yl, df, alpha, p_k[0], p_k[1], p_k[2],
                )
                if not np.isfinite(fk):
                    fk = 1e20
                grad[k] = (fk - f0) / _EPS_FD
            # Chain-rule for parscale
            grad *= parscale
            return f0, grad
        return f_and_grad, True


def optimize_pairwise_density(
    z: np.ndarray,
    df: float,
    alpha: float,
    X: np.ndarray,
    Y: np.ndarray,
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
    ensemble: int = 3,
    max_dist: float = 0.0,
    seed: int = 42,
) -> np.ndarray:
    """Global MLE of (a, b, g) via multi-start L-BFGS-B.

    When Rust is available, each iteration does only 1 FFI crossing
    (NLL + gradient computed together in Rust).
    """
    from weatherisk.density import pairwise_density_optim as _py_optim

    return _py_optim(
        z, df, alpha, X, Y,
        lower_bounds=lower_bounds, upper_bounds=upper_bounds,
        ensemble=ensemble, max_dist=max_dist,
    )


def optimize_local_mle(
    zi: np.ndarray,
    zj: np.ndarray,
    xl: np.ndarray,
    yl: np.ndarray,
    df: float,
    alpha: float,
    ensemble: int = 5,
    seed: int = 42,
    max_boundary_retries: int = 5,
    lower_bounds: tuple[float, float] = (0.01, 0.01),
    upper_bounds: tuple[float, float] = (15.0, 15.0),
) -> np.ndarray:
    """Local MLE of (a, b, γ) from pre-built pair arrays.

    Calls R's optim() via subprocess for exact R parity.
    """
    return _optimize_local_mle_via_r(
        zi, zj, xl, yl, df, alpha, ensemble, seed, max_boundary_retries,
        lower_bounds, upper_bounds,
    )


def clustering_via_r(
    dist_matrix: np.ndarray,
    method: str = "average",
) -> np.ndarray:
    """Call R's hclust via subprocess for exact R parity.

    Uses binary (raw bytes) data transfer to avoid any floating-point
    precision loss from decimal CSV round-trips.
    """
    rscript = _rscript_path()

    n = dist_matrix.shape[0]
    arr = np.ascontiguousarray(dist_matrix, dtype=np.float64)

    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        tmp_in = f.name
        arr.tofile(f)

    r_code = f"""
raw <- readBin("{tmp_in}", what="double", n={n * n}, size=8, endian="little")
m <- matrix(raw, nrow={n}, ncol={n}, byrow=TRUE)
d <- as.dist(t(m), diag=TRUE)
hc <- hclust(d, method="{method}")
cat(hc$merge[,1], "\\n")
cat(hc$merge[,2], "\\n")
cat(sprintf("%.18e", hc$height), "\\n")
"""
    try:
        result = subprocess.run(
            [rscript, "-e", r_code],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.SubprocessError) as exc:
        os.unlink(tmp_in)
        raise RuntimeError(f"R hclust subprocess failed: {exc}") from exc

    os.unlink(tmp_in)
    lines = result.stdout.strip().split("\n")
    if len(lines) < 3:
        raise RuntimeError(
            f"R hclust returned unexpected output ({len(lines)} lines)"
        )

    merge_i = np.array(lines[0].split(), dtype=float)
    merge_j = np.array(lines[1].split(), dtype=float)
    heights = np.array(lines[2].split(), dtype=float)

    # Convert R's merge format to scipy linkage format
    n_merges = len(heights)
    Z = np.zeros((n_merges, 4))
    cluster_sizes = np.ones(n + n_merges)
    for k in range(n_merges):
        ri, rj = int(merge_i[k]), int(merge_j[k])
        # R: negative = singleton (1-indexed), positive = cluster
        si = (-ri - 1) if ri < 0 else (n + ri - 1)
        sj = (-rj - 1) if rj < 0 else (n + rj - 1)
        Z[k, 0] = min(si, sj)
        Z[k, 1] = max(si, sj)
        Z[k, 2] = heights[k]
        Z[k, 3] = cluster_sizes[si] + cluster_sizes[sj]
        cluster_sizes[n + k] = Z[k, 3]

    return Z

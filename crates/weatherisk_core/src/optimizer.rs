/// Full L-BFGS-B optimizer loop in Rust.
///
/// Mirrors the Python `pairwise_density_optim` and `_local_mle_one_cmip6`
/// functions, but runs the entire optimizer inside compiled code so there
/// is only ONE PyO3 boundary crossing per cell.
///
/// Uses `lbfgsb-rs-pure`, a Rust port of the original Fortran L-BFGS-B v3.0
/// (the same algorithm SciPy wraps).

use crate::density::neg_log_likelihood_sum;

// ── Latin Hypercube Sampling ─────────────────────────────────────────────

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

/// Generate `n` Latin Hypercube samples in `d` dimensions, values in [0, 1].
///
/// Simple stratified approach: divide [0,1] into n strata per dimension,
/// randomly permute the strata, then sample uniformly within each stratum.
fn latin_hypercube(n: usize, d: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut samples = vec![vec![0.0; d]; n];

    for dim in 0..d {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.gen_range(0..=i);
            perm.swap(i, j);
        }

        for i in 0..n {
            let stratum = perm[i] as f64;
            let u: f64 = rng.gen();
            samples[i][dim] = (stratum + u) / n as f64;
        }
    }

    samples
}

/// Scale Latin Hypercube samples from [0,1]^d to [lo, hi]^d.
fn scale_samples(samples: &[Vec<f64>], lo: &[f64], hi: &[f64]) -> Vec<Vec<f64>> {
    samples
        .iter()
        .map(|s| {
            s.iter()
                .enumerate()
                .map(|(j, &v)| lo[j] + v * (hi[j] - lo[j]))
                .collect()
        })
        .collect()
}

// ── L-BFGS-B wrapper ────────────────────────────────────────────────────

use lbfgsb_rs_pure::LBFGSB;

/// Result of an L-BFGS-B optimization run.
struct OptResult {
    x: Vec<f64>,
    fun: f64,
}

/// Run L-BFGS-B on the NLL objective with box constraints.
///
/// Internally computes finite-difference gradients (central differences),
/// matching SciPy's `approx_grad=True` behavior.
fn run_lbfgsb(
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    zi: &[f64],
    zj: &[f64],
    x_arr: &[f64],
    y_arr: &[f64],
    df: f64,
    alpha: f64,
    maxiter: usize,
) -> OptResult {
    let n = x0.len();
    let eps = 1e-8; // finite-difference step

    let mut x = x0.to_vec();

    // Build the closure that returns (f, grad)
    let mut f_and_grad = |x_current: &[f64]| -> (f64, Vec<f64>) {
        let fval = neg_log_likelihood_sum(
            zi, zj, x_arr, y_arr, df, alpha,
            x_current[0], x_current[1], x_current[2],
        );

        // Return large value with zero gradient for non-finite objective
        if !fval.is_finite() {
            return (1e20, vec![0.0; n]);
        }

        let mut grad = vec![0.0; n];
        for i in 0..n {
            let h = eps * x_current[i].abs().max(1.0);

            let mut x_plus = x_current.to_vec();
            let mut x_minus = x_current.to_vec();
            x_plus[i] = (x_current[i] + h).min(hi[i]);
            x_minus[i] = (x_current[i] - h).max(lo[i]);

            let fp = neg_log_likelihood_sum(
                zi, zj, x_arr, y_arr, df, alpha,
                x_plus[0], x_plus[1], x_plus[2],
            );
            let fm = neg_log_likelihood_sum(
                zi, zj, x_arr, y_arr, df, alpha,
                x_minus[0], x_minus[1], x_minus[2],
            );

            let actual_h = x_plus[i] - x_minus[i];
            grad[i] = if actual_h.abs() > 0.0 {
                let g = (fp - fm) / actual_h;
                if g.is_finite() { g } else { 0.0 }
            } else {
                0.0
            };
        }

        (fval, grad)
    };

    let mut solver = LBFGSB::new(10) // m=10 (SciPy default)
        .with_max_iter(maxiter)
        .with_pgtol(1e-10);

    let result = solver.minimize(
        &mut x,
        lo,
        hi,
        &mut f_and_grad,
    );

    let fval = match result {
        Ok(sol) => sol.f,
        Err(_) => {
            neg_log_likelihood_sum(zi, zj, x_arr, y_arr, df, alpha, x[0], x[1], x[2])
        }
    };

    OptResult { x, fun: fval }
}

// ── Public API ───────────────────────────────────────────────────────────

/// Build pair arrays from observation matrix z and coordinates X, Y.
///
/// Mirrors the pair assembly in `pairwise_density_optim`:
///   ilist, jlist = np.triu_indices(n_grid, k=1)
///   Xlist = np.repeat(X[ilist] - X[jlist], n_sim)
///   ...
fn build_pair_arrays(
    z: &[f64],    // row-major: z[i * n_sim + t]
    n_grid: usize,
    n_sim: usize,
    x_coords: &[f64],
    y_coords: &[f64],
    max_dist: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut zi = Vec::new();
    let mut zj = Vec::new();
    let mut xl = Vec::new();
    let mut yl = Vec::new();

    let max_dist_sq = if max_dist > 0.0 {
        max_dist * max_dist
    } else {
        f64::INFINITY
    };

    for i in 0..n_grid {
        for j in (i + 1)..n_grid {
            let dx = x_coords[i] - x_coords[j];
            let dy = y_coords[i] - y_coords[j];
            let dist_sq = dx * dx + dy * dy;

            if dist_sq <= max_dist_sq {
                for t in 0..n_sim {
                    zi.push(z[i * n_sim + t]);
                    zj.push(z[j * n_sim + t]);
                    xl.push(dx);
                    yl.push(dy);
                }
            }
        }
    }

    (zi, zj, xl, yl)
}

/// Full multi-start L-BFGS-B optimizer for the global pairwise density MLE.
///
/// Mirrors `density.pairwise_density_optim` — pair assembly,
/// multi-start LHS, gamma-wrapping retry — all in compiled Rust.
pub fn optimize_pairwise_density(
    z: &[f64],
    n_grid: usize,
    n_sim: usize,
    df: f64,
    alpha: f64,
    x_coords: &[f64],
    y_coords: &[f64],
    lower_a: f64,
    lower_b: f64,
    upper_a: f64,
    upper_b: f64,
    ensemble: usize,
    max_dist: f64,
    seed: u64,
) -> [f64; 3] {
    let (zi, zj, xl, yl) =
        build_pair_arrays(z, n_grid, n_sim, x_coords, y_coords, max_dist);

    if zi.is_empty() {
        return [0.0, 0.0, 0.0];
    }

    let lo = [lower_a, lower_b, -std::f64::consts::FRAC_PI_2];
    let hi = [upper_a, upper_b, std::f64::consts::FRAC_PI_2];

    // Latin Hypercube starting points
    let raw_samples = latin_hypercube(ensemble, 3, seed);
    let starts = scale_samples(&raw_samples, &lo, &hi);

    let mut best_val = f64::INFINITY;
    let mut best_par = [starts[0][0], starts[0][1], starts[0][2]];

    for i in 0..ensemble {
        let start = &starts[i];
        let result = run_lbfgsb(start, &lo, &hi, &zi, &zj, &xl, &yl, df, alpha, 10000);

        let mut par = result.x;
        let mut fval = result.fun;

        // Gamma-wrapping retry: if gamma hits ±π/2, flip and re-run
        if (par[2].abs() - std::f64::consts::FRAC_PI_2).abs() < 1e-10 {
            let retry_start = [par[0], par[1], -par[2]];
            let result2 =
                run_lbfgsb(&retry_start, &lo, &hi, &zi, &zj, &xl, &yl, df, alpha, 10000);
            par = result2.x;
            fval = result2.fun;
        }

        if fval < best_val {
            best_val = fval;
            best_par = [par[0], par[1], par[2]];
        }
    }

    best_par
}

/// Full multi-start L-BFGS-B optimizer for the local MLE.
///
/// Takes pre-built pair arrays (zi, zj, xl, yl) and runs the optimizer.
/// Mirrors `_local_mle_one_cmip6` in cmip6_pipeline.py.
pub fn optimize_local_mle(
    zi: &[f64],
    zj: &[f64],
    xl: &[f64],
    yl: &[f64],
    df: f64,
    alpha: f64,
    lower_a: f64,
    lower_b: f64,
    _lower_g: f64,
    upper_a: f64,
    upper_b: f64,
    _upper_g: f64,
    ensemble: usize,
    seed: u64,
) -> [f64; 3] {
    if zi.is_empty() {
        return [1.0, 0.0, 0.0];
    }

    let lo = [lower_a, lower_b, -std::f64::consts::FRAC_PI_2];
    let hi = [upper_a, upper_b, std::f64::consts::FRAC_PI_2];

    let n_starts = ensemble.max(5);
    let raw_samples = latin_hypercube(n_starts, 3, seed);
    let starts = scale_samples(&raw_samples, &lo, &hi);

    let mut best_val = f64::INFINITY;
    let mut best_par = [1.0, 0.0, 0.0];

    for s in 0..ensemble {
        let start = &starts[s];
        let result = run_lbfgsb(start, &lo, &hi, zi, zj, xl, yl, df, alpha, 10000);

        if result.fun < best_val {
            best_val = result.fun;
            best_par = [result.x[0], result.x[1], result.x[2]];
        }
    }

    best_par
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_latin_hypercube_bounds() {
        let samples = latin_hypercube(10, 3, 42);
        assert_eq!(samples.len(), 10);
        for s in &samples {
            assert_eq!(s.len(), 3);
            for &v in s {
                assert!(v >= 0.0 && v <= 1.0, "LHS value out of [0,1]: {}", v);
            }
        }
    }

    #[test]
    fn test_scale_samples() {
        let samples = latin_hypercube(5, 2, 42);
        let lo = vec![1.0, -1.0];
        let hi = vec![3.0, 1.0];
        let scaled = scale_samples(&samples, &lo, &hi);
        for s in &scaled {
            assert!(s[0] >= 1.0 && s[0] <= 3.0);
            assert!(s[1] >= -1.0 && s[1] <= 1.0);
        }
    }

    #[test]
    fn test_build_pair_arrays() {
        let z = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 0.0, 0.0];

        let (zi, zj, xl, yl) = build_pair_arrays(&z, 3, 2, &x, &y, 0.0);
        assert_eq!(zi.len(), 6);
        assert_eq!(zj.len(), 6);
        assert_eq!(xl.len(), 6);
        assert_eq!(yl.len(), 6);
    }

    #[test]
    fn test_build_pair_arrays_with_max_dist() {
        let z = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x = vec![0.0, 1.0, 10.0];
        let y = vec![0.0, 0.0, 0.0];

        let (zi, _zj, _xl, _yl) = build_pair_arrays(&z, 3, 2, &x, &y, 2.0);
        assert_eq!(zi.len(), 2);
    }

    #[test]
    fn test_optimizer_produces_finite() {
        let n_grid = 4;
        let n_sim = 10;
        let mut z = vec![0.0; n_grid * n_sim];
        for i in 0..z.len() {
            z[i] = 0.5 + (i as f64 * 1.7 % 5.0);
        }
        let x = vec![0.0, 1.0, 0.0, 1.0];
        let y = vec![0.0, 0.0, 1.0, 1.0];

        let result = optimize_pairwise_density(
            &z, n_grid, n_sim, 5.0, 1.0, &x, &y, 0.01, 0.01, 15.0, 15.0, 3, 0.0, 42,
        );

        assert!(result[0].is_finite());
        assert!(result[1].is_finite());
        assert!(result[2].is_finite());
        assert!(result[0] >= 0.01);
        assert!(result[0] <= 15.0);
    }

    #[test]
    fn test_local_optimizer_produces_finite() {
        let zi = vec![1.5, 2.0, 1.2, 3.0, 2.5];
        let zj = vec![2.0, 1.5, 1.8, 2.5, 3.0];
        let xl = vec![0.5, 0.3, 0.7, -0.5, -0.3];
        let yl = vec![0.3, 0.5, 0.2, -0.3, -0.5];

        let result = optimize_local_mle(
            &zi, &zj, &xl, &yl,
            5.0, 1.0,
            0.01, 0.0, -std::f64::consts::FRAC_PI_2,
            15.0, 15.0, std::f64::consts::FRAC_PI_2,
            3, 42,
        );

        assert!(result[0].is_finite());
        assert!(result[1].is_finite());
        assert!(result[2].is_finite());
    }
}

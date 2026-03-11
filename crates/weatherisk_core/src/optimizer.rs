/// Full L-BFGS-B optimizer loop in Rust.
///
/// Mirrors the Python `pairwise_density_optim` and `_local_mle_one_cmip6`
/// functions, but runs the entire optimizer inside compiled code so there
/// is only ONE PyO3 boundary crossing per cell.
///
/// Uses `rustimization`, a wrapper around the original Fortran L-BFGS-B.

use crate::density::neg_log_likelihood_sum;
#[cfg(test)]
use crate::grid::{grid_number, koord_num, number_grid};

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

use rustimization::lbfgsb_minimizer::Lbfgsb;
use std::cell::RefCell;

/// Result of an L-BFGS-B optimization run.
struct OptResult {
    x: Vec<f64>,
    fun: f64,
    success: bool,
    status: String,
    iterations: usize,
}

fn eval_nll(
    x_current: &[f64],
    parscale: Option<&[f64]>,
    zi: &[f64],
    zj: &[f64],
    x_arr: &[f64],
    y_arr: &[f64],
    df: f64,
    alpha: f64,
) -> f64 {
    let mut par = [x_current[0], x_current[1], x_current[2]];
    if let Some(scale) = parscale {
        for i in 0..3 {
            par[i] *= scale[i];
        }
    }

    neg_log_likelihood_sum(
        zi, zj, x_arr, y_arr, df, alpha,
        par[0], par[1], par[2],
    )
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
    run_lbfgsb_with_scale(
        x0, lo, hi, None,
        zi, zj, x_arr, y_arr, df, alpha, maxiter,
    )
}

fn run_lbfgsb_with_scale(
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    parscale: Option<&[f64]>,
    zi: &[f64],
    zj: &[f64],
    x_arr: &[f64],
    y_arr: &[f64],
    df: f64,
    alpha: f64,
    maxiter: usize,
) -> OptResult {
    run_lbfgsb_generic(x0, lo, hi, maxiter, |x_current| {
        eval_nll(x_current, parscale, zi, zj, x_arr, y_arr, df, alpha)
    })
}

fn run_lbfgsb_generic<F>(
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    maxiter: usize,
    objective: F,
) -> OptResult
where
    F: FnMut(&[f64]) -> f64,
{
    run_lbfgsb_generic_with_eps(x0, lo, hi, maxiter, 1e-3, objective)
}

fn run_lbfgsb_generic_with_eps<F>(
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    maxiter: usize,
    eps: f64,
    mut objective: F,
) -> OptResult
where
    F: FnMut(&[f64]) -> f64,
{
    run_lbfgsb_generic_with_fgrad(x0, lo, hi, maxiter, |x_current| {
        let fval = objective(x_current);

        if !fval.is_finite() {
            return (1e20, vec![0.0; x_current.len()]);
        }

        let mut grad = vec![0.0; x_current.len()];
        for i in 0..x_current.len() {
            let h = eps;

            let mut x_plus = x_current.to_vec();
            x_plus[i] = (x_current[i] + h).min(hi[i]);
            let mut x_minus = x_current.to_vec();
            x_minus[i] = (x_current[i] - h).max(lo[i]);

            let fp = objective(&x_plus);
            let fm = objective(&x_minus);

            let actual_h = x_plus[i] - x_minus[i];
            grad[i] = if actual_h.abs() > 0.0 {
                let g = (fp - fm) / actual_h;
                if g.is_finite() { g } else { 0.0 }
            } else {
                0.0
            };
        }

        (fval, grad)
    })
}

fn run_lbfgsb_generic_with_fgrad<F>(
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    maxiter: usize,
    mut f_and_grad: F,
) -> OptResult
where
    F: FnMut(&[f64]) -> (f64, Vec<f64>),
{
    let mut x = x0.to_vec();

    struct CachedEval<'a, F>
    where
        F: FnMut(&[f64]) -> (f64, Vec<f64>),
    {
        evaluator: RefCell<&'a mut F>,
        last_x: RefCell<Option<Vec<f64>>>,
        last_fg: RefCell<Option<(f64, Vec<f64>)>>,
    }

    impl<'a, F> CachedEval<'a, F>
    where
        F: FnMut(&[f64]) -> (f64, Vec<f64>),
    {
        fn new(evaluator: &'a mut F) -> Self {
            Self {
                evaluator: RefCell::new(evaluator),
                last_x: RefCell::new(None),
                last_fg: RefCell::new(None),
            }
        }

        fn eval(&self, x: &[f64]) -> (f64, Vec<f64>) {
            let needs_refresh = self
                .last_x
                .borrow()
                .as_ref()
                .map(|last| last.as_slice() != x)
                .unwrap_or(true);

            if needs_refresh {
                let fg = (self.evaluator.borrow_mut())(x);
                *self.last_x.borrow_mut() = Some(x.to_vec());
                *self.last_fg.borrow_mut() = Some(fg);
            }

            self.last_fg
                .borrow()
                .as_ref()
                .cloned()
                .expect("cached optimizer evaluation missing")
        }
    }

    let cache = CachedEval::new(&mut f_and_grad);
    let value = |x_vec: &Vec<f64>| -> f64 { cache.eval(x_vec.as_slice()).0 };
    let gradient = |x_vec: &Vec<f64>| -> Vec<f64> { cache.eval(x_vec.as_slice()).1 };

    let n_params = x.len();
    let mut solver = Lbfgsb::new(&mut x, &value, &gradient);
    solver.set_matric_correction(5);
    solver.set_termination_tolerance(1.0e7);
    solver.set_tolerance(0.0);
    solver.max_iteration(maxiter as u32);
    solver.set_verbosity(-1);
    for i in 0..n_params {
        solver.set_lower_bound(i, lo[i]);
        solver.set_upper_bound(i, hi[i]);
    }
    solver.minimize();

    let fval = cache.eval(&x).0;
    let success = fval.is_finite();
    let status = if success {
        "rustimization-finished".to_string()
    } else {
        "rustimization-nonfinite".to_string()
    };
    let iterations = 0;

    OptResult {
        x,
        fun: fval,
        success,
        status,
        iterations,
    }
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
    lower_g: f64,
    upper_a: f64,
    upper_b: f64,
    upper_g: f64,
    ensemble: usize,
    seed: u64,
) -> [f64; 3] {
    let max_boundary_retries = 5usize;
    let total_starts = ensemble + max_boundary_retries;
    let raw_samples = latin_hypercube(total_starts, 3, seed);
    optimize_local_mle_from_unit_starts(
        zi,
        zj,
        xl,
        yl,
        df,
        alpha,
        lower_a,
        lower_b,
        lower_g,
        upper_a,
        upper_b,
        upper_g,
        ensemble,
        &raw_samples,
    )
}

fn optimize_local_mle_from_unit_starts(
    zi: &[f64],
    zj: &[f64],
    xl: &[f64],
    yl: &[f64],
    df: f64,
    alpha: f64,
    lower_a: f64,
    lower_b: f64,
    lower_g: f64,
    upper_a: f64,
    upper_b: f64,
    upper_g: f64,
    ensemble: usize,
    unit_starts: &[Vec<f64>],
) -> [f64; 3] {
    if zi.is_empty() {
        return [1.0, 0.0, 0.0];
    }

    let lo = [lower_a, lower_b, lower_g];
    let hi = [upper_a, upper_b, upper_g];
    let parscale = [
        (hi[0] - lo[0]) / 100.0,
        (hi[1] - lo[1]) / 100.0,
        (hi[2] - lo[2]) / 100.0,
    ];
    let max_boundary_retries = 5usize;
    let lo_scaled = [lo[0] / parscale[0], lo[1] / parscale[1], lo[2] / parscale[2]];
    let hi_scaled = [hi[0] / parscale[0], hi[1] / parscale[1], hi[2] / parscale[2]];
    let starts = scale_samples(unit_starts, &lo_scaled, &hi_scaled);

    let mut best_val = f64::INFINITY;
    let mut best_par = [1.0, 0.0, 0.0];
    let mut runs_completed = 0usize;
    let mut boundary_retries = 0usize;
    let mut start_idx = 0usize;

    while runs_completed < ensemble && start_idx < starts.len() {
        let start = &starts[start_idx];
        let result = run_lbfgsb_with_scale(
            start,
            &lo_scaled,
            &hi_scaled,
            Some(&parscale),
            zi,
            zj,
            xl,
            yl,
            df,
            alpha,
            10000,
        );

        let mut par = [
            result.x[0] * parscale[0],
            result.x[1] * parscale[1],
            result.x[2] * parscale[2],
        ];
        let mut fval = result.fun;

        if (par[2].abs() - upper_g.abs()).abs() < 1e-10 {
            let retry_start = [
                par[0] / parscale[0],
                par[1] / parscale[1],
                -par[2] / parscale[2],
            ];
            let result2 = run_lbfgsb_with_scale(
                &retry_start,
                &lo_scaled,
                &hi_scaled,
                Some(&parscale),
                zi,
                zj,
                xl,
                yl,
                df,
                alpha,
                10000,
            );
            par = [
                result2.x[0] * parscale[0],
                result2.x[1] * parscale[1],
                result2.x[2] * parscale[2],
            ];
            fval = result2.fun;
        }

        if fval < best_val {
            best_val = fval;
            best_par = par;
        }

        let near_lower = (0..3).any(|i| (par[i] - lo[i]).abs() < 0.01);
        let near_upper = (0..3).any(|i| (par[i] - hi[i]).abs() < 0.01);
        if (near_lower || near_upper) && boundary_retries < max_boundary_retries {
            boundary_retries += 1;
            start_idx += 1;
            continue;
        }

        runs_completed += 1;
        start_idx += 1;
    }

    best_par
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_csv_fixture(csv_text: &str) -> Vec<Vec<String>> {
        csv_text
            .lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split(',')
                    .map(|field| field.trim().trim_matches('"').to_string())
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn load_axis_fixture(csv_text: &str, column_name: &str) -> Vec<f64> {
        let rows = parse_csv_fixture(csv_text);
        let header = csv_text.lines().next().expect("fixture header missing");
        let columns: Vec<&str> = header
            .split(',')
            .map(|field| field.trim().trim_matches('"'))
            .collect();
        let index = columns
            .iter()
            .position(|&name| name == column_name)
            .expect("fixture column missing");
        rows.into_iter()
            .map(|row| row[index].parse::<f64>().expect("invalid fixture value"))
            .collect()
    }

    fn load_named_scalar(csv_text: &str, parameter_name: &str) -> f64 {
        for row in parse_csv_fixture(csv_text) {
            if row[0] == parameter_name {
                return row[1].parse::<f64>().expect("invalid scalar value");
            }
        }
        panic!("missing scalar parameter");
    }

    fn load_simulation_fixture(csv_text: &str) -> Vec<Vec<f64>> {
        parse_csv_fixture(csv_text)
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .skip(1)
                    .map(|value| value.parse::<f64>().expect("invalid simulation value"))
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn load_upper_bounds_from_parameters(csv_text: &str) -> (f64, f64) {
        let rows = parse_csv_fixture(csv_text);
        let mut max_a = f64::NEG_INFINITY;
        let mut max_b = f64::NEG_INFINITY;
        for row in rows {
            let a = row[1].parse::<f64>().expect("invalid a value");
            let b = row[2].parse::<f64>().expect("invalid b value");
            max_a = max_a.max(a);
            max_b = max_b.max(b);
        }
        (max_a + 5.0, 2.0 * max_b)
    }

    fn build_local_pair_arrays(
        sim_data: &[Vec<f64>],
        x_ax: &[f64],
        y_ax: &[f64],
        x: f64,
        y: f64,
        abstand: usize,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let nrow = y_ax.len();
        let ncol = x_ax.len();
        let xy_pos = koord_num(x, y, x_ax, y_ax).expect("grid point missing");
        let (row0, col0) = number_grid(xy_pos, nrow, ncol).expect("index inverse failed");
        let n_sim = sim_data[xy_pos].len();

        let mut sel_grid = Vec::new();
        let radius_sq = (abstand * abstand) as isize;
        for dx in -(abstand as isize)..=(abstand as isize) {
            for dy in -(abstand as isize)..=(abstand as isize) {
                let dist_sq = dx * dx + dy * dy;
                if dist_sq == 0 || dist_sq > radius_sq {
                    continue;
                }

                let sel_col = col0 as isize + dx;
                let sel_row = row0 as isize + dy;
                if sel_row < 0 || sel_row >= nrow as isize || sel_col < 0 || sel_col >= ncol as isize {
                    continue;
                }

                let idx = grid_number(sel_row as usize, sel_col as usize, nrow, ncol)
                    .expect("grid number failed");
                sel_grid.push(idx);
            }
        }

        let mut zi = Vec::new();
        let mut zj = Vec::new();
        let mut xl = Vec::new();
        let mut yl = Vec::new();

        for &neighbor_idx in &sel_grid {
            let dx = x_ax[neighbor_idx / nrow] - x_ax[col0];
            let dy = y_ax[neighbor_idx % nrow] - y_ax[row0];
            for t in 0..n_sim {
                zi.push(sim_data[neighbor_idx][t]);
                zj.push(sim_data[xy_pos][t]);
                xl.push(dx);
                yl.push(dy);
            }
        }

        (zi, zj, xl, yl)
    }

    fn micro_case_objective(case_id: &str, par: &[f64]) -> f64 {
        match case_id {
            "interior_scaled" => {
                ((par[0] - 1.25) / 0.35).powi(2)
                    + ((par[1] - 3.4) / 1.8).powi(2)
                    + ((par[2] + 0.42) / 0.12).powi(2)
            }
            "upper_boundary_b" => {
                (par[0] - 0.45).powi(2)
                    + (par[1] - 10.8).powi(2)
                    + 0.2 * (par[2] + 0.3).powi(2)
            }
            "tilted_valley" => {
                let dx = par[0] - 0.55;
                let dy = par[1] - 6.0;
                let dz = par[2] + 0.65;
                dx.powi(2) + 2.0 * dy.powi(2) + 0.5 * dz.powi(2) + 0.6 * dx * dy
            }
            _ => panic!("unknown micro case"),
        }
    }

    fn micro_case_gradient(case_id: &str, par: &[f64]) -> [f64; 3] {
        match case_id {
            "interior_scaled" => [
                2.0 * (par[0] - 1.25) / 0.35_f64.powi(2),
                2.0 * (par[1] - 3.4) / 1.8_f64.powi(2),
                2.0 * (par[2] + 0.42) / 0.12_f64.powi(2),
            ],
            "upper_boundary_b" => [
                2.0 * (par[0] - 0.45),
                2.0 * (par[1] - 10.8),
                0.4 * (par[2] + 0.3),
            ],
            "tilted_valley" => {
                let dx = par[0] - 0.55;
                let dy = par[1] - 6.0;
                let dz = par[2] + 0.65;
                [
                    2.0 * dx + 0.6 * dy,
                    4.0 * dy + 0.6 * dx,
                    dz,
                ]
            }
            _ => panic!("unknown micro case"),
        }
    }

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

    #[test]
    fn test_local_optimizer_matches_r_min_fixture() {
        let input_rows = parse_csv_fixture(include_str!(
            "../../../tests/reference_data/local_optimizer_min_inputs.csv"
        ));
        let output_rows = parse_csv_fixture(include_str!(
            "../../../tests/reference_data/local_optimizer_min_outputs.csv"
        ));
        let unit_start_rows = parse_csv_fixture(include_str!(
            "../../../tests/reference_data/maximin_lhs_10x3_seed42.csv"
        ));

        assert_eq!(output_rows.len(), 1, "expected exactly one optimizer output row");

        let mut zi = Vec::with_capacity(input_rows.len());
        let mut zj = Vec::with_capacity(input_rows.len());
        let mut xl = Vec::with_capacity(input_rows.len());
        let mut yl = Vec::with_capacity(input_rows.len());

        for row in input_rows {
            zi.push(row[4].parse::<f64>().expect("invalid zi"));
            zj.push(row[5].parse::<f64>().expect("invalid zj"));
            xl.push(row[6].parse::<f64>().expect("invalid xl"));
            yl.push(row[7].parse::<f64>().expect("invalid yl"));
        }

        let output = &output_rows[0];
        let lower_a = output[3].parse::<f64>().expect("invalid lower_a");
        let lower_b = output[4].parse::<f64>().expect("invalid lower_b");
        let lower_g = output[5].parse::<f64>().expect("invalid lower_g");
        let upper_a = output[6].parse::<f64>().expect("invalid upper_a");
        let upper_b = output[7].parse::<f64>().expect("invalid upper_b");
        let upper_g = output[8].parse::<f64>().expect("invalid upper_g");
        let ensemble = output[9].parse::<usize>().expect("invalid ensemble");
        let expected = [
            output[11].parse::<f64>().expect("invalid a_est"),
            output[12].parse::<f64>().expect("invalid b_est"),
            output[13].parse::<f64>().expect("invalid g_est"),
        ];

        let unit_starts = unit_start_rows
            .into_iter()
            .map(|row| {
                vec![
                    row[0].parse::<f64>().expect("invalid u1"),
                    row[1].parse::<f64>().expect("invalid u2"),
                    row[2].parse::<f64>().expect("invalid u3"),
                ]
            })
            .collect::<Vec<_>>();

        let actual = optimize_local_mle_from_unit_starts(
            &zi,
            &zj,
            &xl,
            &yl,
            5.0,
            1.0,
            lower_a,
            lower_b,
            lower_g,
            upper_a,
            upper_b,
            upper_g,
            ensemble,
            &unit_starts,
        );

        let expected_nll = neg_log_likelihood_sum(
            &zi,
            &zj,
            &xl,
            &yl,
            5.0,
            1.0,
            expected[0],
            expected[1],
            expected[2],
        );
        let actual_nll = neg_log_likelihood_sum(
            &zi,
            &zj,
            &xl,
            &yl,
            5.0,
            1.0,
            actual[0],
            actual[1],
            actual[2],
        );

        for i in 0..3 {
            assert!(
                (actual[i] - expected[i]).abs() <= 5e-2,
                "min-fixture mismatch param {i}: actual={}, expected={}, actual_nll={}, expected_nll={}, nll_diff={}",
                actual[i],
                expected[i],
                actual_nll,
                expected_nll,
                actual_nll - expected_nll
            );
        }
    }

    #[test]
    fn test_lbfgsb_matches_r_micro_cases() {
        let rows = parse_csv_fixture(include_str!(
            "../../../tests/reference_data/optim_lbfgsb_micro_cases.csv"
        ));

        for row in rows {
            let case_id = row[0].as_str();
            let x0_real = [
                row[1].parse::<f64>().expect("invalid start_a"),
                row[2].parse::<f64>().expect("invalid start_b"),
                row[3].parse::<f64>().expect("invalid start_g"),
            ];
            let lo_real = [
                row[4].parse::<f64>().expect("invalid lower_a"),
                row[5].parse::<f64>().expect("invalid lower_b"),
                row[6].parse::<f64>().expect("invalid lower_g"),
            ];
            let hi_real = [
                row[7].parse::<f64>().expect("invalid upper_a"),
                row[8].parse::<f64>().expect("invalid upper_b"),
                row[9].parse::<f64>().expect("invalid upper_g"),
            ];
            let expected = [
                row[10].parse::<f64>().expect("invalid opt_a"),
                row[11].parse::<f64>().expect("invalid opt_b"),
                row[12].parse::<f64>().expect("invalid opt_g"),
            ];
            let expected_value = row[13].parse::<f64>().expect("invalid opt_value");

            let parscale = [
                (hi_real[0] - lo_real[0]) / 100.0,
                (hi_real[1] - lo_real[1]) / 100.0,
                (hi_real[2] - lo_real[2]) / 100.0,
            ];
            let x0_scaled = [
                x0_real[0] / parscale[0],
                x0_real[1] / parscale[1],
                x0_real[2] / parscale[2],
            ];
            let lo_scaled = [
                lo_real[0] / parscale[0],
                lo_real[1] / parscale[1],
                lo_real[2] / parscale[2],
            ];
            let hi_scaled = [
                hi_real[0] / parscale[0],
                hi_real[1] / parscale[1],
                hi_real[2] / parscale[2],
            ];

            let result = run_lbfgsb_generic(&x0_scaled, &lo_scaled, &hi_scaled, 10000, |scaled| {
                let real = [
                    scaled[0] * parscale[0],
                    scaled[1] * parscale[1],
                    scaled[2] * parscale[2],
                ];
                micro_case_objective(case_id, &real)
            });
            let actual = [
                result.x[0] * parscale[0],
                result.x[1] * parscale[1],
                result.x[2] * parscale[2],
            ];
            let unscaled_result = run_lbfgsb_generic(&x0_real, &lo_real, &hi_real, 10000, |real| {
                micro_case_objective(case_id, real)
            });
            let unscaled_actual = [unscaled_result.x[0], unscaled_result.x[1], unscaled_result.x[2]];
            let rlike_scaled_result = run_lbfgsb_generic_with_eps(
                &x0_scaled,
                &lo_scaled,
                &hi_scaled,
                10000,
                1e-3,
                |scaled| {
                    let real = [
                        scaled[0] * parscale[0],
                        scaled[1] * parscale[1],
                        scaled[2] * parscale[2],
                    ];
                    micro_case_objective(case_id, &real)
                },
            );
            let rlike_scaled_actual = [
                rlike_scaled_result.x[0] * parscale[0],
                rlike_scaled_result.x[1] * parscale[1],
                rlike_scaled_result.x[2] * parscale[2],
            ];
            let exact_grad_scaled_result = run_lbfgsb_generic_with_fgrad(
                &x0_scaled,
                &lo_scaled,
                &hi_scaled,
                10000,
                |scaled| {
                    let real = [
                        scaled[0] * parscale[0],
                        scaled[1] * parscale[1],
                        scaled[2] * parscale[2],
                    ];
                    let grad_real = micro_case_gradient(case_id, &real);
                    (
                        micro_case_objective(case_id, &real),
                        vec![
                            grad_real[0] * parscale[0],
                            grad_real[1] * parscale[1],
                            grad_real[2] * parscale[2],
                        ],
                    )
                },
            );
            let exact_grad_scaled_actual = [
                exact_grad_scaled_result.x[0] * parscale[0],
                exact_grad_scaled_result.x[1] * parscale[1],
                exact_grad_scaled_result.x[2] * parscale[2],
            ];

            for i in 0..3 {
                assert!(
                    (actual[i] - expected[i]).abs() <= 1e-5,
                    "micro-case {case_id} param {i} mismatch: scaled_actual={}, expected={}, scaled_value={}, expected_value={}, scaled_success={}, scaled_status={}, scaled_iterations={}, unscaled_actual={:?}, unscaled_value={}, unscaled_success={}, unscaled_status={}, unscaled_iterations={}, rlike_scaled_actual={:?}, rlike_scaled_value={}, rlike_scaled_success={}, rlike_scaled_status={}, rlike_scaled_iterations={}, exact_grad_scaled_actual={:?}, exact_grad_scaled_value={}, exact_grad_scaled_success={}, exact_grad_scaled_status={}, exact_grad_scaled_iterations={}",
                    actual[i],
                    expected[i],
                    result.fun,
                    expected_value,
                    result.success,
                    result.status,
                    result.iterations,
                    unscaled_actual,
                    unscaled_result.fun,
                    unscaled_result.success,
                    unscaled_result.status,
                    unscaled_result.iterations,
                    rlike_scaled_actual,
                    rlike_scaled_result.fun,
                    rlike_scaled_result.success,
                    rlike_scaled_result.status,
                    rlike_scaled_result.iterations,
                    exact_grad_scaled_actual,
                    exact_grad_scaled_result.fun,
                    exact_grad_scaled_result.success,
                    exact_grad_scaled_result.status,
                    exact_grad_scaled_result.iterations,
                );
            }

            assert!(
                (result.fun - expected_value).abs() <= 1e-8,
                "micro-case {case_id} value mismatch: actual={}, expected={}",
                result.fun,
                expected_value
            );
        }
    }
}

#![allow(clippy::many_single_char_names)]
//! High-level LBFGS-B solver implementation using the helper modules.
//!
//! This module provides `LBFGSB`, a safe, idiomatic Rust implementation of the
//! high-level driver that uses the ported BLAS, line-search and subalgorithms
//! modules. It relies on the crate-level `Status` and `Solution` types for
//! results and termination codes.

use crate::blas;
use crate::linesearch;
use crate::subalgorithms;
use std::f64;

/// Information passed to iteration callbacks.
#[derive(Debug, Clone)]
pub struct IterationInfo {
    pub iteration: usize,
    pub f: f64,
    pub proj_grad_norm: f64,
    pub n_func_evals: usize,
    pub n_segments: usize,
    pub n_skipped: usize,
    pub n_active: usize,
}

/// Control action returned by iteration callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IterationControl {
    /// Continue optimization
    Continue,
    /// Stop with convergence status
    StopConverged,
    /// Stop with custom reason (will return MaxIter status)
    StopCustom,
}

/// Compact L-BFGS-B solver.
///
/// This struct stores limited-memory correction pairs and options for the
/// optimization run. The implementation here is a high-level orchestrator that
/// uses the other modules for lower-level linear algebra and line-search
/// operations. It aims to follow the reference algorithm behavior while being
/// safe and idiomatic in Rust.
pub struct LBFGSB {
    /// maximum number of corrections to keep (m)
    m: usize,
    /// maximum iterations
    max_iter: usize,
    /// projected gradient tolerance
    pgtol: f64,
    /// whether to emit verbose output
    verbose: bool,
    /// stored s vectors (oldest -> newest)
    s: Vec<Vec<f64>>,
    /// stored y vectors (oldest -> newest)
    y: Vec<Vec<f64>>,
    /// rho = 1 / (y^T s) for each stored pair
    rho: Vec<f64>,

    // Column-major compact storage used by the faithful reference pipeline.
    // These are lazily initialized when a problem dimension `n` is known.
    /// WS: n x m column-major storage for S (each column contiguous)
    ws: Vec<f64>,
    /// WY: n x m column-major storage for Y (each column contiguous)
    wy: Vec<f64>,
    /// SY: m x m column-major storage for S'Y
    sy: Vec<f64>,
    /// SS: m x m column-major storage for S'S
    ss: Vec<f64>,
    /// WT: m x m column-major storage for triangular factor
    wt: Vec<f64>,

    /// current head index (0-based) into circular WS/WY storage
    head: usize,
    /// current itail index (0-based) into circular WS/WY storage
    itail: usize,
    /// number of updates performed so far (1-based count like reference)
    iupdat: usize,
    /// current number of stored columns (<= m)
    col: usize,
    /// current theta scaling
    theta: f64,
}

impl LBFGSB {
    /// Create a new solver storing up to `m` corrections.
    pub fn new(m: usize) -> Self {
        let m = if m == 0 { 1 } else { m };
        LBFGSB {
            m,
            max_iter: 1000,
            pgtol: 1e-6,
            verbose: false,
            s: Vec::new(),
            y: Vec::new(),
            rho: Vec::new(),
            // compact buffers left empty; will be lazily allocated when `minimize` sees n
            ws: Vec::new(),
            wy: Vec::new(),
            sy: Vec::new(),
            ss: Vec::new(),
            wt: Vec::new(),
            head: 0,
            itail: 0,
            iupdat: 0,
            col: 0,
            theta: 1.0,
        }
    }

    /// Set maximum iterations.
    pub fn with_max_iter(mut self, it: usize) -> Self {
        self.max_iter = it;
        self
    }

    /// Set projected gradient tolerance.
    pub fn with_pgtol(mut self, tol: f64) -> Self {
        self.pgtol = tol;
        self
    }

    /// Enable verbose output.
    pub fn with_verbose(mut self, v: bool) -> Self {
        self.verbose = v;
        self
    }

    /// Minimize with iteration callback for custom control.
    ///
    /// The callback receives `IterationInfo` at each iteration and returns
    /// `IterationControl` to continue or stop the optimization.
    ///
    /// This is useful for:
    /// - Custom stopping criteria (like driver2.c)
    /// - Time-based termination (like driver3.c)
    /// - Custom logging or monitoring
    ///
    /// Returns a `crate::Solution` describing the final point, objective value,
    /// iteration count and termination `crate::Status`.
    pub fn minimize_with_callback<F, C>(
        &mut self,
        x: &mut [f64],
        lower: &[f64],
        upper: &[f64],
        f_and_grad: &mut F,
        callback: &mut C,
    ) -> Result<crate::Solution, &'static str>
    where
        F: FnMut(&[f64]) -> (f64, Vec<f64>),
        C: FnMut(&IterationInfo, &[f64]) -> IterationControl,
    {
        let n = x.len();
        if lower.len() != n || upper.len() != n {
            return Err("length mismatch between x and bounds");
        }

        // Validate bounds and project initial x
        for i in 0..n {
            if lower[i] > upper[i] {
                return Err("lower > upper for some index");
            }
            if x[i] < lower[i] {
                x[i] = lower[i];
            } else if x[i] > upper[i] {
                x[i] = upper[i];
            }
        }

        // Build nbd flags similar to the reference (0=unbounded,1=lower,2=both,3=upper)
        let mut nbd = vec![0i32; n];
        for i in 0..n {
            let lo_inf = lower[i].is_infinite() && lower[i].is_sign_negative();
            let hi_inf = upper[i].is_infinite() && upper[i].is_sign_positive();
            nbd[i] = if lo_inf && hi_inf {
                0
            } else if !lo_inf && hi_inf {
                1
            } else if lo_inf && !hi_inf {
                3
            } else {
                2
            };
        }

        // Initialize iwhere using the lightweight cmprlb helper so cauchy/form* can use it.
        let mut iwhere = vec![0_i32; n];
        subalgorithms::cmprlb(x, lower, upper, &nbd, &mut iwhere).map_err(|_| "cmprlb failed")?;

        // Wrap the user function to count evaluations and to be usable by line-search wrappers.
        let eval_count = std::rc::Rc::new(std::cell::Cell::new(0usize));
        let eval_count_cloned = eval_count.clone();
        let mut func = |xx: &[f64]| -> (f64, Vec<f64>) {
            eval_count_cloned.set(eval_count_cloned.get() + 1);
            f_and_grad(xx)
        };

        // initial function/gradient evaluation (through the wrapper)
        let (mut f, mut grad) = func(x);
        if grad.len() != n {
            return Err("gradient length mismatch");
        }

        // Lazily allocate column-major compact buffers for the reference-style pipeline.
        // We allocate WS/WY as n-by-m and SY/SS/WT as m-by-m.
        if self.ws.len() != n * self.m {
            self.ws = vec![0.0; n * self.m];
            self.wy = vec![0.0; n * self.m];
            self.sy = vec![0.0; self.m * self.m];
            self.ss = vec![0.0; self.m * self.m];
            self.wt = vec![0.0; self.m * self.m];
            self.head = 0;
            self.itail = 0;
            self.iupdat = 0;
            self.col = 0;
            self.theta = 1.0;
        }

        // projected gradient and stopping check (use reference projgr)
        let mut pg_norm =
            subalgorithms::projgr(x, lower, upper, &nbd, &grad).map_err(|_| "projgr failed")?;
        // keep a projected-gradient vector available for fallback direction
        let mut pg = vec![0.0; n];

        if self.verbose {
            eprintln!("iter=0 f={:.6e} ||proj_grad||_inf={:.3e}", f, pg_norm);
        }

        if pg_norm <= self.pgtol {
            return Ok(crate::Solution {
                x: x.to_vec(),
                f,
                iterations: 0,
                status: crate::Status::Converged,
            });
        }

        // Preallocate cauchy / subspace workspace once to avoid per-iteration allocations.
        let mut iorder = vec![0i32; n];
        let mut t_work = vec![0.0f64; n];
        let mut d_work = vec![0.0f64; n];
        let mut xcp = vec![0.0f64; n];
        let max_col = self.m; // allocate for worst-case number of columns
        let mut p_work = vec![0.0f64; 2 * max_col];
        let mut c_work = vec![0.0f64; 2 * max_col];
        let mut wbp = vec![0.0f64; 2 * max_col];
        let mut vwrk = vec![0.0f64; 2 * max_col];
        let _kmat = vec![0.0f64; max_col * max_col];
        let _tmat = vec![0.0f64; max_col * max_col];
        let _rvec = vec![0.0f64; max_col];
        let _zvec = vec![0.0f64; max_col];

        // main iteration loop
        let mut tnint: usize = 0;
        let mut nskip: usize = 0;
        let mut nact: usize = 0;
        for iter in 1..=self.max_iter {
            // two-loop recursion to compute H * (-grad)
            let mut q = grad.iter().map(|v| -*v).collect::<Vec<f64>>();
            let col = self.s.len();
            let mut alpha: Vec<f64> = vec![0.0; col];

            // first loop: newest -> oldest
            for i in (0..col).rev() {
                alpha[i] = self.rho[i] * blas::ddot(&self.s[i], &q);
                // q = q - alpha[i] * y[i]
                for (q_item, &y_item) in q.iter_mut().zip(self.y[i].iter()) {
                    *q_item -= alpha[i] * y_item;
                }
            }

            // apply initial H0 scaling
            if col > 0 {
                let last = col - 1;
                let sy = blas::ddot(&self.s[last], &self.y[last]);
                let yy = blas::ddot(&self.y[last], &self.y[last]);
                let gamma = if yy > 0.0 { sy / yy } else { 1.0 };
                for v in q.iter_mut() {
                    *v *= gamma.max(1e-20);
                }
            }

            // second loop: oldest -> newest
            for ((&rho_i, &alpha_i), (y_i, s_i)) in self
                .rho
                .iter()
                .zip(alpha.iter())
                .zip(self.y.iter().zip(self.s.iter()))
                .take(col)
            {
                let beta = rho_i * blas::ddot(y_i, &q);
                for (q_item, &s_item) in q.iter_mut().zip(s_i.iter()) {
                    *q_item += s_item * (alpha_i - beta);
                }
            }

            // By default, pick q as approximate H*(-g)
            let mut d = q;

            // If we have compact memory, compute the Generalized Cauchy Point (GCP)
            // and then run the faithful subspace reconstruction (subsm_full) to get xp.
            if self.col > 0 {
                let col_use = self.col;
                // reuse preallocated buffers (only use the leading slices)
                let p_slice = &mut p_work[0..(2 * col_use)];
                let c_slice = &mut c_work[0..(2 * col_use)];
                let wbp_slice = &mut wbp[0..(2 * col_use)];
                let v_slice = &mut vwrk[0..(2 * col_use)];
                let mut nseg: i32 = 0;
                let mut info: i32 = 0;

                // Call cauchy to compute xcp and c = W'(xcp-x)
                let cauchy_res = subalgorithms::cauchy(
                    n,
                    x,
                    lower,
                    upper,
                    &nbd,
                    &grad,
                    &mut iorder,
                    &mut iwhere,
                    &mut t_work,
                    &mut d_work,
                    &mut xcp,
                    self.m,
                    &self.wy,
                    &self.ws,
                    &self.sy,
                    &self.wt,
                    self.theta,
                    col_use,
                    self.head,
                    p_slice,
                    c_slice,
                    wbp_slice,
                    v_slice,
                    &mut nseg,
                    if self.verbose { 0 } else { -1 },
                    pg_norm,
                    &mut info,
                    f64::EPSILON,
                );

                if cauchy_res.is_ok() {
                    // accumulate segments explored by cauchy for trace comparisons
                    tnint = tnint.saturating_add(nseg as usize);

                    // count active bounds at GCP (iwhere > 0)
                    nact = iwhere.iter().filter(|&&w| w > 0).count();

                    if self.verbose {
                        eprintln!(
                            "iter {}: cauchy nseg={} info={} nact={}",
                            iter, nseg, info, nact
                        );
                    }

                    // Reconstruct full-space subspace minimizer xp using faithful wrapper.
                    let mut xp = vec![0.0f64; n];
                    match subalgorithms::subsm_full(
                        n, self.m, &xcp, lower, upper, &nbd, &grad, &iwhere, &self.ws, &self.wy,
                        &self.sy, &self.wt, self.theta, col_use, self.head, &mut xp,
                    ) {
                        Ok(iword) => {
                            // xp returned; form direction d = xp - x and use it if descent
                            let mut d_xp = vec![0.0f64; n];
                            for i in 0..n {
                                d_xp[i] = xp[i] - x[i];
                            }
                            let dd = blas::ddot(&d_xp, &grad);
                            if dd < 0.0 {
                                // use subspace xp direction
                                d = d_xp;
                                if self.verbose {
                                    eprintln!(
                                        "iter {}: using subspace xp (iword={}) as direction",
                                        iter, iword
                                    );
                                }
                            } else {
                                // fallback to GCP (xcp - x) if it is descent
                                let mut d_gcp = vec![0.0f64; n];
                                for i in 0..n {
                                    d_gcp[i] = xcp[i] - x[i];
                                }
                                if blas::ddot(&d_gcp, &grad) < 0.0 {
                                    d = d_gcp;
                                    if self.verbose {
                                        eprintln!(
                                            "iter {}: subspace xp not descent, using GCP",
                                            iter
                                        );
                                    }
                                } else {
                                    // keep default q-based direction (already in d)
                                    if self.verbose {
                                        eprintln!(
                                            "iter {}: neither xp nor GCP are descent, keep q-direction",
                                            iter
                                        );
                                    }
                                }
                            }
                        }
                        Err(_) => {
                            // subsm_full failed - fallback to using GCP if available
                            for i in 0..n {
                                d[i] = xcp[i] - x[i];
                            }
                            if self.verbose {
                                eprintln!(
                                    "iter {}: subsm_full failed, using GCP as fallback",
                                    iter
                                );
                            }
                        }
                    }
                } else {
                    // cauchy failed: leave d as the two-loop q direction
                    if self.verbose {
                        eprintln!("iter {}: cauchy failed, skipping subspace", iter);
                    }
                }
            } // if self.col > 0

            // ensure descent: d^T grad < 0
            let d_dot_grad = blas::ddot(&d, &grad);
            if d_dot_grad >= 0.0 {
                // fallback to negative projected gradient
                if self.verbose {
                    eprintln!(
                        "iter {}: computed direction not descent (d.g={:.3e}), using -proj_grad",
                        iter, d_dot_grad
                    );
                }
                for i in 0..n {
                    if (x[i] <= lower[i] && grad[i] >= 0.0) || (x[i] >= upper[i] && grad[i] <= 0.0)
                    {
                        pg[i] = 0.0;
                    } else {
                        pg[i] = -grad[i];
                    }
                }
                // copy projected-gradient into d (avoid moving `pg`)
                d[..].copy_from_slice(&pg);
            }

            // perform line-search using direction d
            match linesearch::lnsrlb_search(x, &d, &mut f, &mut grad, lower, upper, &mut func) {
                Ok((x_new, f_new, g_new)) => {
                    let s_vec = sub_vecs(&x_new, x);
                    let y_vec = sub_vecs(&g_new, &grad);
                    let sty = blas::ddot(&s_vec, &y_vec);
                    if sty > 1e-12 {
                        // keep history in the classic data structures
                        self.push_correction(s_vec.clone(), y_vec.clone());
                        // update reference compact matrices (column-major)
                        self.iupdat += 1;
                        let rr = blas::ddot(&y_vec, &y_vec);
                        let dr = blas::ddot(&y_vec, &s_vec);
                        let stp = 1.0f64;
                        let dtd = blas::ddot(&s_vec, &s_vec);
                        let _ = subalgorithms::matupd(
                            n,
                            self.m,
                            &mut self.ws,
                            &mut self.wy,
                            &mut self.sy,
                            &mut self.ss,
                            &s_vec,
                            &y_vec,
                            &mut self.itail,
                            self.iupdat,
                            &mut self.col,
                            &mut self.head,
                            &mut self.theta,
                            rr,
                            dr,
                            stp,
                            dtd,
                        );
                    } else {
                        // record skipped BFGS update for trace parity
                        nskip = nskip.saturating_add(1);
                    }
                    x.copy_from_slice(&x_new);
                    f = f_new;
                    grad = g_new;
                }
                Err(status) => {
                    return Ok(crate::Solution {
                        x: x.to_vec(),
                        f,
                        iterations: iter - 1,
                        status,
                    });
                }
            }

            // stopping check using projected gradient (use reference projgr)
            pg_norm =
                subalgorithms::projgr(x, lower, upper, &nbd, &grad).map_err(|_| "projgr failed")?;

            if self.verbose {
                eprintln!(
                    "{:5} {:5} {:5} {:5} {:5} {:5} {:12.5e} {:12.5e}",
                    n,
                    iter,
                    eval_count.get(),
                    tnint,
                    nskip,
                    nact,
                    pg_norm,
                    f
                );
            }

            // Call user callback for custom control
            let info = IterationInfo {
                iteration: iter,
                f,
                proj_grad_norm: pg_norm,
                n_func_evals: eval_count.get(),
                n_segments: tnint,
                n_skipped: nskip,
                n_active: nact,
            };

            match callback(&info, x) {
                IterationControl::Continue => {
                    // Check default convergence
                    if pg_norm <= self.pgtol {
                        return Ok(crate::Solution {
                            x: x.to_vec(),
                            f,
                            iterations: iter,
                            status: crate::Status::Converged,
                        });
                    }
                }
                IterationControl::StopConverged => {
                    return Ok(crate::Solution {
                        x: x.to_vec(),
                        f,
                        iterations: iter,
                        status: crate::Status::Converged,
                    });
                }
                IterationControl::StopCustom => {
                    return Ok(crate::Solution {
                        x: x.to_vec(),
                        f,
                        iterations: iter,
                        status: crate::Status::MaxIter,
                    });
                }
            }
        }

        // max iterations reached
        Ok(crate::Solution {
            x: x.to_vec(),
            f,
            iterations: self.max_iter,
            status: crate::Status::MaxIter,
        })
    }

    /// Minimize the function starting from `x` with box bounds.
    ///
    /// - `x` is updated in-place and contains the initial guess.
    /// - `lower` / `upper` are bound vectors (use +/-INFINITY for unbounded).
    /// - `f_and_grad` is a callback returning (f, grad_vec) for any `x`.
    ///
    /// Returns a `crate::Solution` describing the final point, objective value,
    /// iteration count and termination `crate::Status`.
    pub fn minimize<F>(
        &mut self,
        x: &mut [f64],
        lower: &[f64],
        upper: &[f64],
        f_and_grad: &mut F,
    ) -> Result<crate::Solution, &'static str>
    where
        F: FnMut(&[f64]) -> (f64, Vec<f64>),
    {
        let n = x.len();
        if lower.len() != n || upper.len() != n {
            return Err("length mismatch between x and bounds");
        }

        // Validate bounds and project initial x
        for i in 0..n {
            if lower[i] > upper[i] {
                return Err("lower > upper for some index");
            }
            if x[i] < lower[i] {
                x[i] = lower[i];
            } else if x[i] > upper[i] {
                x[i] = upper[i];
            }
        }

        // Build nbd flags similar to the reference (0=unbounded,1=lower,2=both,3=upper)
        let mut nbd = vec![0i32; n];
        for i in 0..n {
            let lo_inf = lower[i].is_infinite() && lower[i].is_sign_negative();
            let hi_inf = upper[i].is_infinite() && upper[i].is_sign_positive();
            nbd[i] = if lo_inf && hi_inf {
                0
            } else if !lo_inf && hi_inf {
                1
            } else if lo_inf && !hi_inf {
                3
            } else {
                2
            };
        }

        // Initialize iwhere using the lightweight cmprlb helper so cauchy/form* can use it.
        let mut iwhere = vec![0_i32; n];
        subalgorithms::cmprlb(x, lower, upper, &nbd, &mut iwhere).map_err(|_| "cmprlb failed")?;

        // Wrap the user function to count evaluations and to be usable by line-search wrappers.
        let eval_count = std::rc::Rc::new(std::cell::Cell::new(0usize));
        let eval_count_cloned = eval_count.clone();
        let mut func = |xx: &[f64]| -> (f64, Vec<f64>) {
            eval_count_cloned.set(eval_count_cloned.get() + 1);
            f_and_grad(xx)
        };

        // initial function/gradient evaluation (through the wrapper)
        let (mut f, mut grad) = func(x);
        if grad.len() != n {
            return Err("gradient length mismatch");
        }

        // Lazily allocate column-major compact buffers for the reference-style pipeline.
        // We allocate WS/WY as n-by-m and SY/SS/WT as m-by-m.
        if self.ws.len() != n * self.m {
            self.ws = vec![0.0; n * self.m];
            self.wy = vec![0.0; n * self.m];
            self.sy = vec![0.0; self.m * self.m];
            self.ss = vec![0.0; self.m * self.m];
            self.wt = vec![0.0; self.m * self.m];
            self.head = 0;
            self.itail = 0;
            self.iupdat = 0;
            self.col = 0;
            self.theta = 1.0;
        }

        // projected gradient and stopping check (use reference projgr)
        let mut pg_norm =
            subalgorithms::projgr(x, lower, upper, &nbd, &grad).map_err(|_| "projgr failed")?;
        // keep a projected-gradient vector available for fallback direction
        let mut pg = vec![0.0; n];

        if self.verbose {
            eprintln!("iter=0 f={:.6e} ||proj_grad||_inf={:.3e}", f, pg_norm);
        }

        if pg_norm <= self.pgtol {
            return Ok(crate::Solution {
                x: x.to_vec(),
                f,
                iterations: 0,
                status: crate::Status::Converged,
            });
        }

        // Preallocate cauchy / subspace workspace once to avoid per-iteration allocations.
        let mut iorder = vec![0i32; n];
        let mut t_work = vec![0.0f64; n];
        let mut d_work = vec![0.0f64; n];
        let mut xcp = vec![0.0f64; n];
        let max_col = self.m; // allocate for worst-case number of columns
        let mut p_work = vec![0.0f64; 2 * max_col];
        let mut c_work = vec![0.0f64; 2 * max_col];
        let mut wbp = vec![0.0f64; 2 * max_col];
        let mut vwrk = vec![0.0f64; 2 * max_col];
        let _kmat = vec![0.0f64; max_col * max_col];
        let _tmat = vec![0.0f64; max_col * max_col];
        let _rvec = vec![0.0f64; max_col];
        let _zvec = vec![0.0f64; max_col];

        // main iteration loop
        // Counters matching the C reference semantics:
        // - Tit : total number of (outer) iterations (we use `iter`)
        // - Tnf : total number of function/gradient evaluations (tracked via eval_count)
        // - Tnint: total number of segments explored during cauchy (accumulated)
        // - Skip : number of BFGS updates skipped (when s'y is too small)
        // - Nact : number of active bounds at the current GCP (computed after cauchy)
        let mut tnint: usize = 0;
        let mut nskip: usize = 0;
        let mut nact: usize = 0;
        for iter in 1..=self.max_iter {
            // two-loop recursion to compute H * (-grad)
            let mut q = grad.iter().map(|v| -*v).collect::<Vec<f64>>();
            let col = self.s.len();
            let mut alpha: Vec<f64> = vec![0.0; col];

            // first loop: newest -> oldest
            for i in (0..col).rev() {
                alpha[i] = self.rho[i] * blas::ddot(&self.s[i], &q);
                // q = q - alpha[i] * y[i]
                for (q_item, &y_item) in q.iter_mut().zip(self.y[i].iter()) {
                    *q_item -= alpha[i] * y_item;
                }
            }

            // apply initial H0 scaling
            if col > 0 {
                let last = col - 1;
                let sy = blas::ddot(&self.s[last], &self.y[last]);
                let yy = blas::ddot(&self.y[last], &self.y[last]);
                let gamma = if yy > 0.0 { sy / yy } else { 1.0 };
                for v in q.iter_mut() {
                    *v *= gamma.max(1e-20);
                }
            }

            // second loop: oldest -> newest
            for ((&rho_i, &alpha_i), (y_i, s_i)) in self
                .rho
                .iter()
                .zip(alpha.iter())
                .zip(self.y.iter().zip(self.s.iter()))
                .take(col)
            {
                let beta = rho_i * blas::ddot(y_i, &q);
                for (q_item, &s_item) in q.iter_mut().zip(s_i.iter()) {
                    *q_item += s_item * (alpha_i - beta);
                }
            }

            // By default, pick q as approximate H*(-g)
            let mut d = q;

            // If we have compact memory, compute the Generalized Cauchy Point (GCP)
            // and then run the faithful subspace reconstruction (subsm_full) to get xp.
            if self.col > 0 {
                let col_use = self.col;
                // reuse preallocated buffers (only use the leading slices)
                let p_slice = &mut p_work[0..(2 * col_use)];
                let c_slice = &mut c_work[0..(2 * col_use)];
                let wbp_slice = &mut wbp[0..(2 * col_use)];
                let v_slice = &mut vwrk[0..(2 * col_use)];
                let mut nseg: i32 = 0;
                let mut info: i32 = 0;

                // Call cauchy to compute xcp and c = W'(xcp-x)
                let cauchy_res = subalgorithms::cauchy(
                    n,
                    x,
                    lower,
                    upper,
                    &nbd,
                    &grad,
                    &mut iorder,
                    &mut iwhere,
                    &mut t_work,
                    &mut d_work,
                    &mut xcp,
                    self.m,
                    &self.wy,
                    &self.ws,
                    &self.sy,
                    &self.wt,
                    self.theta,
                    col_use,
                    self.head,
                    p_slice,
                    c_slice,
                    wbp_slice,
                    v_slice,
                    &mut nseg,
                    if self.verbose { 0 } else { -1 },
                    pg_norm,
                    &mut info,
                    f64::EPSILON,
                );

                if cauchy_res.is_ok() {
                    // accumulate segments explored by cauchy for trace comparisons
                    tnint = tnint.saturating_add(nseg as usize);

                    // count active bounds at GCP (iwhere > 0)
                    nact = iwhere.iter().filter(|&&w| w > 0).count();

                    if self.verbose {
                        eprintln!(
                            "iter {}: cauchy nseg={} info={} nact={}",
                            iter, nseg, info, nact
                        );
                    }

                    // Reconstruct full-space subspace minimizer xp using faithful wrapper.
                    let mut xp = vec![0.0f64; n];
                    match subalgorithms::subsm_full(
                        n, self.m, &xcp, lower, upper, &nbd, &grad, &iwhere, &self.ws, &self.wy,
                        &self.sy, &self.wt, self.theta, col_use, self.head, &mut xp,
                    ) {
                        Ok(iword) => {
                            // xp returned; form direction d = xp - x and use it if descent
                            let mut d_xp = vec![0.0f64; n];
                            for i in 0..n {
                                d_xp[i] = xp[i] - x[i];
                            }
                            let dd = blas::ddot(&d_xp, &grad);
                            if dd < 0.0 {
                                // use subspace xp direction
                                d = d_xp;
                                if self.verbose {
                                    eprintln!(
                                        "iter {}: using subspace xp (iword={}) as direction",
                                        iter, iword
                                    );
                                }
                            } else {
                                // fallback to GCP (xcp - x) if it is descent
                                let mut d_gcp = vec![0.0f64; n];
                                for i in 0..n {
                                    d_gcp[i] = xcp[i] - x[i];
                                }
                                if blas::ddot(&d_gcp, &grad) < 0.0 {
                                    d = d_gcp;
                                    if self.verbose {
                                        eprintln!(
                                            "iter {}: subspace xp not descent, using GCP",
                                            iter
                                        );
                                    }
                                } else {
                                    // keep default q-based direction (already in d)
                                    if self.verbose {
                                        eprintln!(
                                            "iter {}: neither xp nor GCP are descent, keep q-direction",
                                            iter
                                        );
                                    }
                                }
                            }
                        }
                        Err(_) => {
                            // subsm_full failed - fallback to using GCP if available
                            for i in 0..n {
                                d[i] = xcp[i] - x[i];
                            }
                            if self.verbose {
                                eprintln!(
                                    "iter {}: subsm_full failed, using GCP as fallback",
                                    iter
                                );
                            }
                        }
                    }
                } else {
                    // cauchy failed: leave d as the two-loop q direction
                    if self.verbose {
                        eprintln!("iter {}: cauchy failed, skipping subspace", iter);
                    }
                }
            } // if self.col > 0

            // ensure descent: d^T grad < 0
            let d_dot_grad = blas::ddot(&d, &grad);
            if d_dot_grad >= 0.0 {
                // fallback to negative projected gradient
                if self.verbose {
                    eprintln!(
                        "iter {}: computed direction not descent (d.g={:.3e}), using -proj_grad",
                        iter, d_dot_grad
                    );
                }
                for i in 0..n {
                    if (x[i] <= lower[i] && grad[i] >= 0.0) || (x[i] >= upper[i] && grad[i] <= 0.0)
                    {
                        pg[i] = 0.0;
                    } else {
                        pg[i] = -grad[i];
                    }
                }
                // copy projected-gradient into d (avoid moving `pg`)
                d[..].copy_from_slice(&pg);
            }

            // perform line-search using direction d
            match linesearch::lnsrlb_search(x, &d, &mut f, &mut grad, lower, upper, &mut func) {
                Ok((x_new, f_new, g_new)) => {
                    let s_vec = sub_vecs(&x_new, x);
                    let y_vec = sub_vecs(&g_new, &grad);
                    let sty = blas::ddot(&s_vec, &y_vec);
                    if sty > 1e-12 {
                        // keep history in the classic data structures
                        self.push_correction(s_vec.clone(), y_vec.clone());
                        // update reference compact matrices (column-major)
                        self.iupdat += 1;
                        let rr = blas::ddot(&y_vec, &y_vec);
                        let dr = blas::ddot(&y_vec, &s_vec);
                        let stp = 1.0f64;
                        let dtd = blas::ddot(&s_vec, &s_vec);
                        let _ = subalgorithms::matupd(
                            n,
                            self.m,
                            &mut self.ws,
                            &mut self.wy,
                            &mut self.sy,
                            &mut self.ss,
                            &s_vec,
                            &y_vec,
                            &mut self.itail,
                            self.iupdat,
                            &mut self.col,
                            &mut self.head,
                            &mut self.theta,
                            rr,
                            dr,
                            stp,
                            dtd,
                        );
                    } else {
                        // record skipped BFGS update for trace parity
                        nskip = nskip.saturating_add(1);
                    }
                    x.copy_from_slice(&x_new);
                    f = f_new;
                    grad = g_new;
                }
                Err(status) => {
                    return Ok(crate::Solution {
                        x: x.to_vec(),
                        f,
                        iterations: iter - 1,
                        status,
                    });
                }
            }

            // stopping check using projected gradient (use reference projgr)
            pg_norm =
                subalgorithms::projgr(x, lower, upper, &nbd, &grad).map_err(|_| "projgr failed")?;

            if self.verbose {
                // Print a C-reference-like per-iteration summary line for easy comparison:
                // N Tit Tnf Tnint Skip Nact Projg F
                // where:
                //   N    : problem dimension
                //   Tit  : iteration number (iter)
                //   Tnf  : function/gradient evals (eval_count)
                //   Tnint: accumulated cauchy segments (tnint)
                //   Skip : skipped updates (nskip)
                //   Nact : active bounds at GCP (nact)
                //   Projg: infinity norm of projected gradient
                //   F    : current objective value
                eprintln!(
                    "{:5} {:5} {:5} {:5} {:5} {:5} {:12.5e} {:12.5e}",
                    n,
                    iter,
                    eval_count.get(),
                    tnint,
                    nskip,
                    nact,
                    pg_norm,
                    f
                );
            }

            if pg_norm <= self.pgtol {
                return Ok(crate::Solution {
                    x: x.to_vec(),
                    f,
                    iterations: iter,
                    status: crate::Status::Converged,
                });
            }
        }

        // max iterations reached
        Ok(crate::Solution {
            x: x.to_vec(),
            f,
            iterations: self.max_iter,
            status: crate::Status::MaxIter,
        })
    }

    /// Push a new correction pair (s, y) into memory; evict oldest if capacity reached.
    fn push_correction(&mut self, s_vec: Vec<f64>, y_vec: Vec<f64>) {
        let sty = blas::ddot(&s_vec, &y_vec);
        if sty == 0.0 {
            return;
        }
        let rho_val = 1.0 / sty;
        if self.s.len() == self.m {
            self.s.remove(0);
            self.y.remove(0);
            self.rho.remove(0);
        }
        self.s.push(s_vec);
        self.y.push(y_vec);
        self.rho.push(rho_val);
    }
}

// ---------------- helper functions ----------------

// `projected_gradient` removed — use `subalgorithms::projgr(x, l, u, &nbd, &g)`
// which returns the projected-gradient infinity norm. The previous helper
// populated a vector with the projected components; the solver now computes
// the norm directly using the reference-style `projgr`.

/// Vector difference a - b
fn sub_vecs(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai - bi).collect()
}

/// Infinity norm of vector (internal helper; intentionally unused in some builds)
fn _inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, &x| m.max(x.abs()))
}

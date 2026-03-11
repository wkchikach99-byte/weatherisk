//! Driver 3: Time-controlled driver for L-BFGS-B
//!
//! This example shows how to terminate a run after some prescribed CPU time
//! has elapsed, similar to driver3.c from the original L-BFGS-B package.
//!
//! Features demonstrated:
//! - Time-based termination
//! - Custom iteration monitoring
//! - Handling early termination gracefully

use lbfgsb_rs_pure::{IterationControl, IterationInfo, LBFGSB};
use std::time::Instant;

fn main() {
    println!("     Solving sample problem.");
    println!("      (f = 0.0 at the optimal solution.)");

    let n = 1000;
    let m = 10;

    // Set bounds: odd-numbered variables (0, 2, 4, ...) have bounds [1, 100]
    // even-numbered variables (1, 3, 5, ...) have bounds [-100, 100]
    let mut lower = vec![0.0; n];
    let mut upper = vec![0.0; n];

    for i in (0..n).step_by(2) {
        lower[i] = 1.0;
        upper[i] = 100.0;
    }

    for i in (1..n).step_by(2) {
        lower[i] = -100.0;
        upper[i] = 100.0;
    }

    // Starting point: all variables at 3.0
    let mut x = vec![3.0; n];

    // Create solver with custom parameters
    // We suppress the default stopping tests by setting pgtol to 0
    // and providing our own criteria via callback
    let mut solver = LBFGSB::new(m).with_pgtol(0.0).with_verbose(false);

    // Set time limit (in seconds) - matching driver3.c's 0.2 seconds
    let time_limit_secs = 0.2;
    let start_time = Instant::now();

    // Define custom stopping criteria
    let max_func_evals = 900;
    let custom_gradient_tolerance = 1e-10;

    // Chained Rosenbrock function (as in C driver3.c)
    let sol = solver
        .minimize_with_callback(
            &mut x,
            &lower,
            &upper,
            &mut |x: &[f64]| {
                // f = (x[0] - 1)^2 + 4 * sum_{i=2}^{n} (x[i-1] - x[i-2]^2)^2
                let mut f = (x[0] - 1.0) * (x[0] - 1.0);
                for i in 1..n {
                    let diff = x[i] - x[i - 1] * x[i - 1];
                    f += 4.0 * diff * diff;
                }

                // Compute gradient
                let mut g = vec![0.0; n];
                let mut t1 = x[1] - x[0] * x[0];
                g[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * t1;

                for i in 1..(n - 1) {
                    let t2 = t1;
                    t1 = x[i + 1] - x[i] * x[i];
                    g[i] = 8.0 * t2 - 16.0 * x[i] * t1;
                }
                g[n - 1] = 8.0 * t1;

                (f, g)
            },
            &mut |info: &IterationInfo, x: &[f64]| {
                let elapsed = start_time.elapsed().as_secs_f64();

                // Check if time limit exceeded first
                if elapsed > time_limit_secs {
                    println!("STOP: CPU EXCEEDING THE TIME LIMIT.");
                    println!(
                        "At latest iterate   f = {:.2e}, |proj g| = {:.2e}",
                        info.f, info.proj_grad_norm
                    );

                    // Print some of the final X values
                    if n < 30 {
                        print!("Latest iterate X = ");
                        for &xi in x.iter() {
                            print!("{:.2e} ", xi);
                        }
                        println!();
                    } else {
                        print!("Latest iterate X = ");
                        for &xi in x.iter().take(10) {
                            print!("{:.2e} ", xi);
                        }
                        println!("...");
                        println!("...(suppressing some output)...");
                        print!("... ");
                        for &xi in x.iter().skip(n - 10) {
                            print!("{:.2e} ", xi);
                        }
                        println!();
                    }

                    return IterationControl::StopCustom;
                }

                // Print iteration information
                println!(
                    "Iterate {:5}, nfg = {:5}, f = {:.2e}, |proj g| = {:.2e}",
                    info.iteration, info.n_func_evals, info.f, info.proj_grad_norm
                );

                // Custom stopping test 1: Terminate if total number of f and g
                // evaluations exceeds limit
                if info.n_func_evals >= max_func_evals {
                    println!("STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT");
                    print_final_x(x, n);
                    return IterationControl::StopCustom;
                }

                // Custom stopping test 2: Terminate if |proj g|/(1+|f|) < tolerance
                if info.proj_grad_norm <= (info.f.abs() + 1.0) * custom_gradient_tolerance {
                    println!("STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL");
                    print_final_x(x, n);
                    return IterationControl::StopConverged;
                }

                IterationControl::Continue
            },
        )
        .unwrap();

    println!("\n           * * *");
    println!("Final status: {:?}", sol.status);
    println!("Iterations: {}", sol.iterations);
    println!("F(x) = {:.9e}", sol.f);
    println!(
        "Total elapsed time: {:.3} seconds",
        start_time.elapsed().as_secs_f64()
    );
}

fn print_final_x(x: &[f64], n: usize) {
    if n < 30 {
        print!("Final X = ");
        for &xi in x.iter() {
            print!("{:.2e} ", xi);
        }
        println!();
    } else {
        print!("Final X = ");
        for &xi in x.iter().take(10) {
            print!("{:.2e} ", xi);
        }
        println!("...");
        println!("...(suppressing some output)...");
        print!("... ");
        for &xi in x.iter().skip(n - 10) {
            print!("{:.2e} ", xi);
        }
        println!();
    }
}

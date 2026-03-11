use lbfgsb_rs_pure::LBFGSB;

fn main() {
    println!("     Solving sample problem (Rosenbrock test fcn).");
    println!("      (f = 0.0 at the optimal solution.)");

    let n = 25;
    let m = 5;

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

    // Create solver with m=5 limited memory corrections
    // factr = 1e7 in C code, which is used as: factr * machine_eps
    // machine_eps ~= 2.22e-16, so factr*eps ~= 2.22e-9
    // This is used for function value convergence check
    // We'll use pgtol for projected gradient tolerance
    let mut solver = LBFGSB::new(m)
        .with_pgtol(1e-5)
        .with_verbose(true)
        .with_max_iter(1000);

    // Chained Rosenbrock function (as in C driver1.c)
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
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
        })
        .unwrap();

    println!("\n           * * *");
    println!("Tit   = total number of iterations");
    println!("Tnf   = total number of function evaluations");
    println!("Tnint = total number of segments explored during Cauchy searches");
    println!("Skip  = number of BFGS updates skipped");
    println!("Nact  = number of active bounds at final generalized Cauchy point");
    println!("Projg = norm of the final projected gradient");
    println!("F     = final function value");
    println!("           * * *");
    println!(
        "Status: {:?}, Iterations: {}, F(x) = {:.9e}",
        sol.status, sol.iterations, sol.f
    );
}

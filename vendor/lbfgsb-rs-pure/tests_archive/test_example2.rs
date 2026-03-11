use lbfgsb_rs::{LBFGSB, Status};

fn main() {
    // Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
    // Minimum at (1, 1)
    let mut x = vec![0.0, 0.0];
    let lower = vec![-2.0, -2.0];
    let upper = vec![2.0, 2.0];

    let mut solver = LBFGSB::new(5).with_verbose(true);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            let a = 1.0 - x[0];
            let b = x[1] - x[0] * x[0];
            let f = a * a + 100.0 * b * b;
            let g = vec![
                -2.0 * a - 400.0 * x[0] * b,
                200.0 * b,
            ];
            (f, g)
        })
        .unwrap();

    println!("\nFinal Status: {:?}", sol.status);
    println!("x: {:?}", sol.x);
    println!("f: {}", sol.f);
    println!("iterations: {}", sol.iterations);
    
    if sol.status == Status::Converged {
        println!("Test passed!");
    } else {
        println!("Test failed - did not converge");
    }
}

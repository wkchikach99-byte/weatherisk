use lbfgsb_rs::{LBFGSB, Status};

fn main() {
    // Simple linear function with a bound constraint
    // f(x) = x, x in [0, 2]
    // Minimum at x = 0
    let mut x = vec![1.0];  // Start at 1.0
    let lower = vec![0.0];
    let upper = vec![2.0];

    let mut solver = LBFGSB::new(5).with_verbose(true).with_max_iter(10);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            (x[0], vec![1.0])  // f=x, grad=1
        })
        .unwrap();

    println!("\nFinal Status: {:?}", sol.status);
    println!("x: {:?}", sol.x);
    println!("f: {}", sol.f);
    println!("iterations: {}", sol.iterations);
    
    if sol.status == Status::Converged && (sol.x[0] - 0.0).abs() < 1e-4 {
        println!("Test passed!");
    } else {
        println!("Test result: x should be near 0.0");
    }
}

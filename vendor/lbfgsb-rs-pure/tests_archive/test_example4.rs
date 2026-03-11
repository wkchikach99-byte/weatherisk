use lbfgsb_rs::{LBFGSB, Status};

fn main() {
    // Unconstrained quadratic
    let mut x = vec![5.0];
    let lower = vec![f64::NEG_INFINITY];
    let upper = vec![f64::INFINITY];

    let mut solver = LBFGSB::new(5).with_verbose(true).with_max_iter(20);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            let f = x[0] * x[0];
            let g = vec![2.0 * x[0]];
            (f, g)
        })
        .unwrap();

    println!("\nFinal Status: {:?}", sol.status);
    println!("x: {:?}", sol.x);
    println!("f: {}", sol.f);
    println!("iterations: {}", sol.iterations);
    
    if sol.status == Status::Converged {
        println!("Test passed!");
    }
}

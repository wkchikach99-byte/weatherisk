use lbfgsb_rs::{LBFGSB, Status};

fn main() {
    // Simple quadratic: f(x) = (x - 0.5)^2, minimum at x = 0.5
    let mut x = vec![0.0];
    let lower = vec![-1.0];
    let upper = vec![1.0];

    let mut solver = LBFGSB::new(5);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            let dx = x[0] - 0.5;
            (dx * dx, vec![2.0 * dx])
        })
        .unwrap();

    println!("Status: {:?}", sol.status);
    println!("x: {:?}", sol.x);
    println!("f: {}", sol.f);
    println!("iterations: {}", sol.iterations);
    
    assert_eq!(sol.status, Status::Converged);
    assert!((sol.x[0] - 0.5).abs() < 1e-6);
    println!("\nTest passed!");
}

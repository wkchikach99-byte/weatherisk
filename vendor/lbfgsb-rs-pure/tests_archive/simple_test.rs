use lbfgsb_rs::LBFGSB;

fn main() {
    // Simple 2D quadratic: f(x,y) = (x-1)^2 + (y-2)^2
    // Minimum at (1, 2) with f=0
    let mut x = vec![0.0, 0.0];
    let lower = vec![-10.0, -10.0];
    let upper = vec![10.0, 10.0];

    let mut solver = LBFGSB::new(5).with_pgtol(1e-8).with_verbose(true);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            let f = (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2);
            let g = vec![2.0 * (x[0] - 1.0), 2.0 * (x[1] - 2.0)];
            (f, g)
        })
        .unwrap();

    println!("\nFinal: x=[{:.10}, {:.10}], f={:.15e}", sol.x[0], sol.x[1], sol.f);
    println!("Expected: x=[1.0, 2.0], f=0.0");
    println!("Error: dx={:.3e}, dy={:.3e}", (sol.x[0]-1.0).abs(), (sol.x[1]-2.0).abs());
}

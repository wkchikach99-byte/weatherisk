use lbfgsb_rs::{LBFGSB, Status};

fn main() {
    let mut x = vec![0.0, 0.0];
    let lower = vec![-2.0, -2.0];
    let upper = vec![2.0, 2.0];

    let mut call_count = 0;
    let mut solver = LBFGSB::new(5).with_verbose(false).with_max_iter(2);
    let sol = solver
        .minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
            let a = 1.0 - x[0];
            let b = x[1] - x[0] * x[0];
            let f = a * a + 100.0 * b * b;
            let g = vec![
                -2.0 * a - 400.0 * x[0] * b,
                200.0 * b,
            ];
            call_count += 1;
            println!("Call {}: x={:?}, f={}, g={:?}", call_count, x, f, g);
            (f, g)
        })
        .unwrap();

    println!("\nFinal: {:?}, x={:?}", sol.status, sol.x);
}

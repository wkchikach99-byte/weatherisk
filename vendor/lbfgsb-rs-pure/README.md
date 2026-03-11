# lbfgsb-rs-pure

A safe, idiomatic Rust port of the L-BFGS-B optimization algorithm (version 3.0).

## Overview

L-BFGS-B is a limited-memory quasi-Newton optimization algorithm for solving large-scale bound-constrained optimization problems. This crate provides a faithful Rust implementation of the algorithm, ported from the C reference implementation.

The algorithm is particularly well-suited for:
- Large-scale optimization (thousands to millions of variables)
- Problems with simple box constraints (lower/upper bounds on variables)
- Cases where computing the full Hessian matrix is impractical

## Features

- ✅ **Complete implementation** - All core functionality ported and verified
- ✅ **Safe Rust** - No unsafe code, bounds-checked operations
- ✅ **Idiomatic API** - Builder pattern, closures, Result types
- ✅ **Iteration callbacks** - Custom stopping criteria and monitoring
- ✅ **Well-tested** - Unit tests and examples validate correctness
- ✅ **Documented** - Comprehensive documentation and examples

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
lbfgsb-rs-pure = "0.1.0"
```

### Basic Usage

```rust
use lbfgsb_rs_pure::LBFGSB;

fn main() {
    let n = 25; // problem size
    let m = 5;  // number of corrections to approximate Hessian

    // Set bounds: [1, 100] for even indices, [-100, 100] for odd
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

    // Starting point
    let mut x = vec![3.0; n];

    // Create solver
    let mut solver = LBFGSB::new(m)
        .with_pgtol(1e-5)
        .with_max_iter(1000);

    // Minimize (Rosenbrock function)
    let sol = solver.minimize(&mut x, &lower, &upper, &mut |x: &[f64]| {
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
    }).unwrap();

    println!("Status: {:?}", sol.status);
    println!("Iterations: {}", sol.iterations);
    println!("F(x) = {:.9e}", sol.f);
}
```

### Advanced: Custom Stopping Criteria

```rust
use lbfgsb_rs::{LBFGSB, IterationControl};

let mut solver = LBFGSB::new(5).with_pgtol(0.0); // Disable default stopping

let sol = solver.minimize_with_callback(
    &mut x,
    &lower,
    &upper,
    &mut |x| (f, grad), // Your objective function
    &mut |info, x| {
        // Custom stopping: max 100 evaluations or relative gradient < 1e-10
        if info.n_func_evals >= 100 {
            return IterationControl::StopCustom;
        }
        if info.proj_grad_norm <= (info.f.abs() + 1.0) * 1e-10 {
            return IterationControl::StopConverged;
        }
        IterationControl::Continue
    },
).unwrap();
```

### Time-Based Termination

```rust
use std::time::Instant;
use lbfgsb_rs::IterationControl;

let start = Instant::now();
let time_limit = 0.5; // seconds

let sol = solver.minimize_with_callback(
    &mut x, &lower, &upper,
    &mut |x| (f, grad),
    &mut |info, _| {
        if start.elapsed().as_secs_f64() > time_limit {
            IterationControl::StopCustom
        } else {
            IterationControl::Continue
        }
    },
).unwrap();
```

## Examples

Three example drivers demonstrate different features:

```bash
# Basic usage
cargo run --example driver1

# Custom stopping criteria (max evaluations, relative gradient)
cargo run --example driver2

# Time-controlled optimization
cargo run --example driver3

# Run all examples
./scripts/run_all_examples.sh
```

See [`docs/EXAMPLES.md`](docs/EXAMPLES.md) for detailed documentation.

## API Documentation

### `LBFGSB`

Main solver struct with builder pattern configuration:

- `new(m: usize)` - Create solver with `m` limited-memory corrections
- `with_pgtol(tol: f64)` - Set projected gradient tolerance (default: 1e-6)
- `with_max_iter(n: usize)` - Set maximum iterations (default: 1000)
- `with_verbose(v: bool)` - Enable verbose output (default: false)

### Methods

- **`minimize`** - Simple API with default stopping criteria
- **`minimize_with_callback`** - Advanced API with iteration callbacks

### `IterationInfo`

Information passed to callbacks:

```rust
pub struct IterationInfo {
    pub iteration: usize,        // Current iteration
    pub f: f64,                   // Objective value
    pub proj_grad_norm: f64,      // ||projected gradient||_inf
    pub n_func_evals: usize,      // Total f/g evaluations
    pub n_segments: usize,        // Cauchy search segments
    pub n_skipped: usize,         // Skipped BFGS updates
    pub n_active: usize,          // Active bounds at GCP
}
```

### `IterationControl`

Control flow returned by callbacks:

- `Continue` - Continue optimization
- `StopConverged` - Stop with converged status
- `StopCustom` - Stop with custom reason

### `Solution`

Result returned by solver:

```rust
pub struct Solution {
    pub x: Vec<f64>,        // Final point
    pub f: f64,              // Final objective value
    pub iterations: usize,   // Number of iterations
    pub status: Status,      // Termination status
}
```

### `Status`

Termination status:

- `Converged` - Gradient tolerance satisfied
- `MaxIter` - Maximum iterations reached
- `LineSearchFailure` - Line search failed
- `NumericalFailure` - Numerical issues detected

## Testing

```bash
# Run unit tests
cargo test

# Run with verbose output
cargo test -- --nocapture

# Build release version
cargo build --release
```

All 13 unit tests pass, covering BLAS operations, line search, linear algebra, and subalgorithms.

## Project Structure

```
lbfgsb-rs-pure/
├── src/
│   ├── lib.rs              # Public API and types
│   ├── solver.rs           # High-level LBFGS-B solver
│   ├── blas.rs             # Level-1 BLAS operations
│   ├── linpack.rs          # Cholesky factorization
│   ├── linesearch.rs       # Line search routines
│   └── subalgorithms.rs    # Core L-BFGS-B subalgorithms
├── examples/
│   ├── driver1.rs          # Basic usage example
│   ├── driver2.rs          # Custom stopping criteria
│   └── driver3.rs          # Time-controlled optimization
├── docs/
│   ├── EXAMPLES.md         # Example documentation
│   ├── COMPARISON.md       # C vs Rust comparison
│   ├── PORTING_STATUS.md   # Porting details
│   └── DRIVER_PORT_SUMMARY.md
├── scripts/
│   ├── run_all_examples.sh # Run all examples
│   └── verify.sh           # Verification script
├── tests_archive/          # Archived development tests
└── L-BFGS-B-C/             # Original C reference
```

## Algorithm References

- R. H. Byrd, P. Lu and J. Nocedal. *A Limited Memory Algorithm for Bound Constrained Optimization* (1995), SIAM Journal on Scientific and Statistical Computing, 16, 5, pp. 1190-1208.
- C. Zhu, R. H. Byrd and J. Nocedal. *L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization* (1997), ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550-560.
- J.L. Morales and J. Nocedal. *L-BFGS-B: Remark on Algorithm 778* (2011), ACM Transactions on Mathematical Software, Vol 38, Num. 1.

More information: [L-BFGS-B Wikipedia](http://en.wikipedia.org/wiki/L-BFGS-B:_Optimization_subject_to_simple_bounds)

## Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- **[EXAMPLES.md](docs/EXAMPLES.md)** - Detailed guide to all example drivers
- **[COMPARISON.md](docs/COMPARISON.md)** - C vs Rust implementation comparison
- **[PORTING_STATUS.md](docs/PORTING_STATUS.md)** - Porting progress and notes
- **[DRIVER_PORT_SUMMARY.md](docs/DRIVER_PORT_SUMMARY.md)** - Driver port technical details

## License

This crate and the repository are provided under the BSD 3-Clause License. See the top-level `LICENSE` file for the full license text and copyright notice. The original L-BFGS-B C reference implementation is published under the BSD 3-Clause License as well.

The original L-BFGS-B algorithm was developed by:
- Ciyou Zhu (1994, with revisions in 1996)
- In collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal
- Version 3.0 algorithmic updates (2011) by J. L. Morales

This Rust port maintains algorithmic fidelity while providing a modern, safe API.

## Citation

If you use this software in research, please cite the original L-BFGS-B papers listed above.

## Contributing

Contributions are welcome! Areas for potential improvement:

- Performance optimization and BLAS integration
- Additional examples and benchmarks
- Extended test coverage
- Documentation improvements

## Acknowledgments

This Rust port is based on the [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C) version by Stephen Becker, which itself is based on the original Fortran code.

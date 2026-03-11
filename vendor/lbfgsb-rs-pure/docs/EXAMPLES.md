# Example Drivers

This directory contains three example drivers that demonstrate different features of the L-BFGS-B Rust implementation.

## driver1.rs - Basic Usage

The simplest example demonstrating basic usage of the L-BFGS-B solver.

**Features demonstrated:**
- Basic solver setup with default convergence criteria
- Simple bounds on variables
- Extended Rosenbrock test function (25 variables)
- Standard convergence with gradient tolerance

**Run:**
```bash
cargo run --example driver1
```

This example uses the high-level `minimize` API with default stopping criteria (projected gradient tolerance).

## driver2.rs - Custom Stopping Criteria

Shows how to implement custom termination conditions and monitor optimization progress.

**Features demonstrated:**
- Custom stopping criteria via iteration callbacks
- Monitoring function evaluation count
- Custom gradient tolerance: `|proj g| / (1 + |f|) < 1e-10`
- Maximum function evaluation limit (99 evaluations)
- Iteration-by-iteration output

**Run:**
```bash
cargo run --example driver2
```

This example corresponds to `driver2.c` from the original L-BFGS-B-C package. It demonstrates:
- Suppressing default stopping tests (by setting `pgtol = 0`)
- Using `minimize_with_callback` for custom control
- Implementing relative gradient tolerance
- Printing custom iteration information

## driver3.rs - Time-Controlled Optimization

Demonstrates how to impose a time limit on optimization runs.

**Features demonstrated:**
- Time-based termination (0.2 second limit)
- Custom stopping criteria combined with time checks
- Larger problem size (1000 variables)
- Handling early termination gracefully
- Printing partial solution vectors for large problems

**Run:**
```bash
cargo run --example driver3
```

This example corresponds to `driver3.c` from the original L-BFGS-B-C package. It shows:
- Using `Instant::now()` to track elapsed time
- Checking time in the iteration callback
- Combining time-based and iteration-based stopping criteria
- Terminating optimization cleanly when time limit is exceeded

## API Comparison

### Simple API (driver1)

```rust
let sol = solver.minimize(&mut x, &lower, &upper, &mut |x| {
    // compute f and gradient
    (f, grad)
}).unwrap();
```

### Callback API (driver2, driver3)

```rust
let sol = solver.minimize_with_callback(
    &mut x,
    &lower,
    &upper,
    &mut |x| {
        // compute f and gradient
        (f, grad)
    },
    &mut |info: &IterationInfo, x: &[f64]| {
        // custom logic for monitoring and stopping
        // return Continue, StopConverged, or StopCustom
        IterationControl::Continue
    }
).unwrap();
```

## IterationInfo Structure

The callback receives an `IterationInfo` struct with the following fields:

- `iteration: usize` - Current iteration number
- `f: f64` - Current objective function value
- `proj_grad_norm: f64` - Infinity norm of projected gradient
- `n_func_evals: usize` - Total function/gradient evaluations
- `n_segments: usize` - Total segments explored in Cauchy searches
- `n_skipped: usize` - Number of BFGS updates skipped (due to low curvature)
- `n_active: usize` - Number of active bounds at generalized Cauchy point

## IterationControl Return Values

Callbacks must return one of:

- `IterationControl::Continue` - Continue optimization
- `IterationControl::StopConverged` - Stop with converged status
- `IterationControl::StopCustom` - Stop with custom reason (MaxIter status)

## Comparison with C Reference

All three drivers replicate the behavior of their C counterparts:

| Driver | Problem Size | Feature | C Version | Rust Version |
|--------|-------------|---------|-----------|--------------|
| driver1 | n=25 | Basic | `driver1.c` | `driver1.rs` |
| driver2 | n=25 | Custom stop | `driver2.c` | `driver2.rs` |
| driver3 | n=1000 | Time limit | `driver3.c` | `driver3.rs` |

The Rust versions produce comparable results to the C reference implementation, with minor numerical differences due to floating-point ordering and implementation details.

## Notes

- The callback API provides fine-grained control over the optimization process
- Custom stopping criteria can combine multiple conditions (evaluations, time, gradient, etc.)
- The `IterationInfo` struct provides comprehensive monitoring data
- Time-based termination is useful for real-time applications with hard deadlines
- All examples use the chained Rosenbrock function as the test problem
# Archived Test Files

This directory contains legacy test files that were used during the initial development and porting process. These files are kept for historical reference but are not part of the active test suite.

## Files

### simple_test.rs
Early simple test case used during initial development.

### test_debug.rs
Debug-focused test case used for troubleshooting during the port.

### test_example.rs, test_example2.rs, test_example3.rs, test_example4.rs
Various test cases created during the development process to verify specific functionality.

## Note

These files are archived because:
- They were standalone test binaries rather than proper Rust tests
- The functionality they tested is now covered by:
  - Unit tests in `src/` modules
  - Example drivers in `examples/`
  - Cargo test suite (`cargo test`)

## Active Testing

For current testing, use:

```bash
# Run all unit tests
cargo test

# Run examples
cargo run --example driver1
cargo run --example driver2
cargo run --example driver3

# Or use the convenience script
../scripts/run_all_examples.sh
```

## Restoration

If you need to restore any of these tests:
1. Move the `.rs` file to the project root or `examples/` directory
2. Compile with `rustc <file>.rs` or add to `examples/` in `Cargo.toml`
3. Or convert to proper integration tests in a `tests/` directory
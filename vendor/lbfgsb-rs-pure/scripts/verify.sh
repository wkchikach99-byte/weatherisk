#!/bin/bash
echo "=== L-BFGS-B Rust Port Verification ==="
echo ""
echo "1. Building release version..."
cargo build --release 2>&1 | grep -E "(Finished|error)"
echo ""
echo "2. Running unit tests..."
cargo test --release --quiet 2>&1 | grep "test result"
echo ""
echo "3. Running example (Driver1)..."
cargo run --release --example driver1 2>&1 | tail -5
echo ""
echo "=== Verification Complete ==="

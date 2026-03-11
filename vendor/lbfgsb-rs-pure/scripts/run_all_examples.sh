#!/bin/bash
# Script to run all three L-BFGS-B example drivers

set -e

echo "=================================="
echo "Running L-BFGS-B Rust Examples"
echo "=================================="
echo ""

echo "1. Running driver1 (basic usage)..."
echo "-----------------------------------"
cargo run --release --example driver1
echo ""
echo ""

echo "2. Running driver2 (custom stopping criteria)..."
echo "-------------------------------------------------"
cargo run --release --example driver2
echo ""
echo ""

echo "3. Running driver3 (time-controlled optimization)..."
echo "----------------------------------------------------"
cargo run --release --example driver3
echo ""
echo ""

echo "=================================="
echo "All examples completed successfully!"
echo "=================================="

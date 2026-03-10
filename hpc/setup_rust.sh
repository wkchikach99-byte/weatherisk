#!/bin/bash
# ============================================================
#  Install Rust toolchain + build weatherisk_core on albedo HPC
#
#  Run this ONCE on a login node (not via SLURM):
#
#    ssh albedo
#    cd /albedo/work/user/wikchi001/weatherisk
#    bash hpc/setup_rust.sh
#
#  After this, the Rust backend is available automatically —
#  backend.py auto-detects `import weatherisk_core`.
#
#  To rebuild after code changes:
#    source ~/.cargo/env
#    source .venv/bin/activate
#    maturin develop --release -m crates/weatherisk_core/Cargo.toml
# ============================================================

set -euo pipefail

WORK_DIR="/albedo/work/user/wikchi001/weatherisk"
cd "${WORK_DIR}"

echo "============================================"
echo "  Rust backend setup for weatherisk"
echo "  $(date)"
echo "============================================"

# ── 1. Install Rust (userspace, no root needed) ──────────────
if command -v rustc &>/dev/null; then
    echo "  Rust already installed: $(rustc --version)"
else
    echo "  Installing Rust via rustup ..."
    # Use minimal profile + single-threaded extraction to avoid
    # resource limits on login nodes (RLIMIT_NPROC causes threadpool panic)
    export RUSTUP_IO_THREADS=1
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --profile minimal
    source "${HOME}/.cargo/env"
    echo "  Installed: $(rustc --version)"
fi

# Make sure cargo is on PATH for this session
source "${HOME}/.cargo/env"

# ── 2. Activate Python venv ─────────────────────────────────
source "${WORK_DIR}/.venv/bin/activate"
echo "  Python: $(python --version) ($(which python))"

# ── 3. Install maturin ──────────────────────────────────────
if python -m pip show maturin &>/dev/null; then
    echo "  maturin already installed"
else
    echo "  Installing maturin ..."
    python -m pip install --quiet maturin
fi

# ── 4. Verify crates/ directory exists ───────────────────────
if [ ! -f "${WORK_DIR}/crates/weatherisk_core/Cargo.toml" ]; then
    echo ""
    echo "  ERROR: crates/weatherisk_core/Cargo.toml not found!"
    echo "  Copy it from your local machine first:"
    echo ""
    echo "    scp -r crates albedo:/albedo/work/user/wikchi001/weatherisk/"
    echo ""
    exit 1
fi

# ── 5. Build the Rust extension ─────────────────────────────
echo ""
echo "  Building weatherisk_core (release mode) ..."
echo "  This takes ~1-2 minutes on the first build."
echo ""
maturin develop --release -m crates/weatherisk_core/Cargo.toml

# ── 6. Verify ────────────────────────────────────────────────
echo ""
echo "  Verifying import ..."
python -c "
import weatherisk_core as rc
print(f'  weatherisk_core loaded: {rc.__file__}')
print(f'  Functions: {[x for x in dir(rc) if not x.startswith(\"_\")]}')
"

echo ""
echo "============================================"
echo "  Setup complete! The Rust backend will be"
echo "  used automatically by weatherisk."
echo "============================================"

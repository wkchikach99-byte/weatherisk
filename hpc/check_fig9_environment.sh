#!/bin/bash

set -euo pipefail

WORK_DIR="${WORK_DIR:-/albedo/work/user/wikchi001/weatherisk}"

echo "============================================"
echo "  Figure 9 albedo preflight"
echo "  $(date)"
echo "============================================"

if [ ! -d "${WORK_DIR}" ]; then
    echo "ERROR: work directory not found: ${WORK_DIR}"
    exit 1
fi

cd "${WORK_DIR}"

echo "Host: $(hostname)"
echo "Work dir: ${WORK_DIR}"

module purge 2>/dev/null || true
module load python/3.11.7 2>/dev/null || module load python/3.11 2>/dev/null || module load python/3.10.4 2>/dev/null || true

if [ -f "${WORK_DIR}/.venv/bin/activate" ]; then
    source "${WORK_DIR}/.venv/bin/activate"
else
    echo "ERROR: virtual environment not found at ${WORK_DIR}/.venv"
    exit 1
fi

echo "Python: $(python3 --version) ($(which python3))"

if command -v git >/dev/null 2>&1 && [ -d .git ]; then
    echo "Git HEAD: $(git rev-parse --short HEAD)"
    if [ -n "$(git status --porcelain)" ]; then
        echo "WARNING: remote worktree has uncommitted changes"
    else
        echo "Git status: clean"
    fi
fi

python3 - <<'PY'
import importlib
import os
import sys

mods = ["weatherisk", "numpy", "scipy", "pandas", "statsmodels", "xarray", "netCDF4"]
for name in mods:
    module = importlib.import_module(name)
    print(f"import ok: {name} -> {getattr(module, '__file__', '<built-in>')}")

try:
    rc = importlib.import_module("weatherisk_core")
    print(f"import ok: weatherisk_core -> {rc.__file__}")
except Exception as exc:
    print(f"WARNING: weatherisk_core not importable: {exc}")

from weatherisk.backend import _USE_RUST
print(f"backend autodetect: {'Rust' if _USE_RUST else 'Python'}")

required = [
    "scripts/reproduce_fig9.py",
    "hpc/run_fig9_clean.slurm",
    "crates/weatherisk_core/Cargo.toml",
]
for path in required:
    if not os.path.exists(path):
        raise SystemExit(f"missing required path: {path}")
    print(f"path ok: {path}")
PY

python3 scripts/reproduce_fig9.py --help >/dev/null
echo "CLI contract: reproduce_fig9.py --help ok"

echo "============================================"
echo "  Preflight passed"
echo "============================================"
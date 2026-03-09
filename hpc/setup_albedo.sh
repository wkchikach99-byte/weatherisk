#!/bin/bash
# ============================================================
# Setup script for AWI Albedo HPC
# Run this ONCE after cloning the repo on Albedo.
#
# Usage:
#   ssh wikchi001@albedo0.awi.de
#   cd /work/wikchi001
#   git clone <repo-url> weatherisk
#   cd weatherisk
#   bash hpc/setup_albedo.sh
# ============================================================

set -euo pipefail

echo "======================================"
echo "  Setting up weatherisk on Albedo HPC"
echo "======================================"

WORK_DIR="${PWD}"
echo "  Working directory: ${WORK_DIR}"

# ── 1. Create directories ──────────────────────────────────────
mkdir -p "${WORK_DIR}/data/cmip6"
mkdir -p "${WORK_DIR}/output/cmip6_fig9"
mkdir -p "${WORK_DIR}/logs"

# ── 2. Check available modules ─────────────────────────────────
echo ""
echo "  Checking available modules ..."
echo "  ─────────────────────────────"

# Try common module names on Albedo
for mod in python python/3.10 python/3.11 python3 anaconda3 conda miniconda; do
    if module avail "${mod}" 2>&1 | grep -q "${mod}"; then
        echo "  ✅ Module available: ${mod}"
    fi
done

# ── 3. Load Python ─────────────────────────────────────────────
echo ""
echo "  Loading Python environment ..."

# Albedo typically uses module system. Try common options.
module purge 2>/dev/null || true

if module load python/3.11 2>/dev/null; then
    echo "  ✅ Loaded python/3.11"
elif module load python/3.10 2>/dev/null; then
    echo "  ✅ Loaded python/3.10"
elif module load anaconda3 2>/dev/null; then
    echo "  ✅ Loaded anaconda3"
elif module load conda 2>/dev/null; then
    echo "  ✅ Loaded conda"
else
    echo "  ⚠️  Could not load a Python module."
    echo "     Try: module avail | grep -i python"
    echo "     Then edit this script accordingly."
fi

# ── 4. Create virtual environment ──────────────────────────────
echo ""
if [ ! -d "${WORK_DIR}/.venv" ]; then
    echo "  Creating virtual environment ..."
    python3 -m venv "${WORK_DIR}/.venv"
    echo "  ✅ Created .venv"
else
    echo "  ⏭️  .venv already exists"
fi

source "${WORK_DIR}/.venv/bin/activate"
echo "  Activated: $(python3 --version) at $(which python3)"

# ── 5. Install dependencies ───────────────────────────────────
echo ""
echo "  Installing dependencies ..."
pip install --upgrade pip setuptools wheel --quiet
pip install -e ".[dev]" --quiet
echo "  ✅ Dependencies installed"

# ── 6. Verify imports ──────────────────────────────────────────
echo ""
echo "  Verifying critical imports ..."
python3 -c "
import numpy; print(f'  numpy     {numpy.__version__}')
import scipy; print(f'  scipy     {scipy.__version__}')
import xarray; print(f'  xarray    {xarray.__version__}')
import statsmodels; print(f'  statsmodels {statsmodels.__version__}')
import weatherisk; print(f'  weatherisk  OK')
" && echo "  ✅ All imports OK" || echo "  ⚠️ Some imports failed"

# ── 7. Check for pre-staged CMIP6 data ────────────────────────
echo ""
echo "  Searching for pre-staged AWI-ESM-1-1-LR data ..."
bash "${WORK_DIR}/hpc/find_cmip6_data.sh"

echo ""
echo "======================================"
echo "  Setup complete!"
echo ""
echo "  Activate environment:"
echo "    source ${WORK_DIR}/.venv/bin/activate"
echo ""
echo "  Submit Figure 9 job:"
echo "    sbatch hpc/run_fig9.slurm"
echo ""
echo "  Or run interactively:"
echo "    python scripts/reproduce_fig9.py --workers 16"
echo "======================================"

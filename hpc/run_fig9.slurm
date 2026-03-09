#!/bin/bash
#SBATCH --job-name=fig9_cmip6
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/fig9_%j.out
#SBATCH --error=logs/fig9_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wikchi001@awi.de
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h

# ============================================================
#  SLURM job: Reproduce Figure 9 from Contzen et al. (2025)
#
#  AWI-ESM-1-1-LR monthly precipitation → STL → annual maxima
#  → GEV → Fréchet → LEC/EDC clustering → cluster maps
#
#  Submit:
#    sbatch hpc/run_fig9.slurm
#
#  Monitor:
#    squeue -u wikchi001
#    tail -f logs/fig9_*.out
#
#  NOTE: Edit --account if required by your allocation.
#        Check with: sacctmgr show assoc user=wikchi001
#        Then add:   #SBATCH --account=<project_id>
# ============================================================

echo "============================================"
echo "  Job started: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================"

# ── Setup ──────────────────────────────────────────────────────
WORK_DIR="/albedo/work/user/wikchi001/weatherisk"
cd "${WORK_DIR}"

mkdir -p logs

# Load Python
module purge 2>/dev/null || true
module load python/3.11.7 2>/dev/null || module load python/3.11 2>/dev/null || \
    module load python/3.10.4 2>/dev/null || true

# Activate venv
if [ -f "${WORK_DIR}/.venv/bin/activate" ]; then
    source "${WORK_DIR}/.venv/bin/activate"
fi

echo "  Python: $(python3 --version) ($(which python3))"
echo "  Working dir: ${WORK_DIR}"
echo ""

# Set threading — each MLE worker gets 1 thread to avoid oversubscription
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Force unbuffered output so progress lines appear in real time
export PYTHONUNBUFFERED=1

# ── Run pipeline ──────────────────────────────────────────────
python3 scripts/reproduce_fig9.py \
    --data-dir /albedo/pool/CMIP6/data/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212 \
    --workers ${SLURM_CPUS_PER_TASK} \
    --output-dir output/cmip6_fig9 \
    --dpi 300 \
    "$@"

echo ""
echo "============================================"
echo "  Job finished: $(date)"
echo "============================================"

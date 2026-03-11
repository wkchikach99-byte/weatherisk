#!/bin/bash
#SBATCH --job-name=fig9_clean
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --time=12:00:00
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk/logs/fig9_clean_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk/logs/fig9_clean_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wikchi001@awi.de
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h

set -euo pipefail

WORK_DIR="/albedo/work/user/wikchi001/weatherisk"
DATA_DIR_DEFAULT="/albedo/pool/CMIP6/data/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212"

FIG9_RUN_ID="${FIG9_RUN_ID:-cmip6_fig9_clean_${SLURM_JOB_ID}}"
FIG9_OUTPUT_DIR="${FIG9_OUTPUT_DIR:-output/${FIG9_RUN_ID}}"
FIG9_DATA_DIR="${FIG9_DATA_DIR:-${DATA_DIR_DEFAULT}}"

cd "${WORK_DIR}"
mkdir -p logs "${FIG9_OUTPUT_DIR}"

module purge 2>/dev/null || true
module load python/3.11.7 2>/dev/null || module load python/3.11 2>/dev/null || module load python/3.10.4 2>/dev/null || true

if [ -f "${WORK_DIR}/.venv/bin/activate" ]; then
    source "${WORK_DIR}/.venv/bin/activate"
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1

echo "============================================"
echo "  Job started: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  CPUs: ${SLURM_CPUS_PER_TASK}"
echo "  Python: $(python3 --version) ($(which python3))"
echo "  Data dir: ${FIG9_DATA_DIR}"
echo "  Output dir: ${FIG9_OUTPUT_DIR}"
echo "============================================"

python3 scripts/reproduce_fig9.py \
    --data-dir "${FIG9_DATA_DIR}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --output-dir "${FIG9_OUTPUT_DIR}" \
    --dpi 300

echo ""
echo "============================================"
echo "  Job finished: $(date)"
echo "  Output dir: ${FIG9_OUTPUT_DIR}"
echo "============================================"
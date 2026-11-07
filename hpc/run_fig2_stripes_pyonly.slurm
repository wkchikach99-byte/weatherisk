#!/bin/bash
#SBATCH --job-name=fig2_stripes
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/fig2_stripes_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/fig2_stripes_%j.err
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h

set -euo pipefail

BASE="/albedo/work/user/wikchi001/weatherisk_pyonly"
SEED="${SEED:-42}"
OUT_DIR="${BASE}/output/fig2_stripes_seed${SEED}_${SLURM_JOB_ID}"

cd "${BASE}/src"
mkdir -p "${BASE}/logs" "${OUT_DIR}"

module purge 2>/dev/null || true
module load python/3.11.7
source "${BASE}/venv/bin/activate"

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1

if command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is on PATH; refusing to run the R-free reproduction" >&2
    exit 2
fi

echo "start $(date -u +%FT%TZ) node $(hostname) job ${SLURM_JOB_ID} seed ${SEED}"
echo "git $(cat "${BASE}/src/GIT_COMMIT" 2>/dev/null || echo unknown)"

/usr/bin/time -v python scripts/reproduce_fig3.py \
    --seed "${SEED}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --output "${OUT_DIR}/fig2_stripes.png" \
    --output-dir "${OUT_DIR}" ${EXTRA_ARGS:-}

echo "end $(date -u +%FT%TZ) output ${OUT_DIR}"

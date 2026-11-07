#!/bin/bash
#SBATCH --job-name=fig9_pyonly
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --time=12:00:00
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/fig9_pyonly_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/fig9_pyonly_%j.err
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h

set -euo pipefail

BASE="/albedo/work/user/wikchi001/weatherisk_pyonly"
DATA_DIR="${FIG9_DATA_DIR:-/albedo/pool/CMIP6/data/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212}"
OUT_DIR="${FIG9_OUTPUT_DIR:-${BASE}/output/fig9_${SLURM_JOB_ID}}"

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

echo "start $(date -u +%FT%TZ) node $(hostname) job ${SLURM_JOB_ID}"
echo "extra args: ${FIG9_EXTRA_ARGS:-none}"
echo "git $(cat "${BASE}/src/GIT_COMMIT" 2>/dev/null || echo unknown)"
python -c "import weatherisk, numpy, scipy, statsmodels; print('weatherisk', weatherisk.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__)"

/usr/bin/time -v python scripts/reproduce_fig9.py \
    --data-dir "${DATA_DIR}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --output-dir "${OUT_DIR}" \
    --no-download \
    --dpi 300 \
    ${FIG9_EXTRA_ARGS:-}

echo "end $(date -u +%FT%TZ) output ${OUT_DIR}"

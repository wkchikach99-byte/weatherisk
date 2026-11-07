#!/bin/bash
#SBATCH --job-name=edc_prep_sweep
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/edc_prep_sweep_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/edc_prep_sweep_%j.err
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h
# k_EDC under the STL and year-range readings the Extremes paper leaves open.
set -euo pipefail
BASE="/albedo/work/user/wikchi001/weatherisk_pyonly"
DATA=/albedo/pool/CMIP6/data/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212
module purge 2>/dev/null || true
module load python/3.11.7
source "${BASE}/venv/bin/activate"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${BASE}/src"
echo "host $(hostname) start $(date -Is) job ${SLURM_JOB_ID}"
/usr/bin/time -v python hpc/r_reference/edc_preprocessing_sweep.py "${DATA}" --workers 16
echo "end $(date -Is)"

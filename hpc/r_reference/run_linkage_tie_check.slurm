#!/bin/bash
#SBATCH --job-name=linkage_ties
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/linkage_ties_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_pyonly/logs/linkage_ties_%j.err
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h

set -euo pipefail
BASE="/albedo/work/user/wikchi001/weatherisk_pyonly"
R_OUT="${R_OUT:?set R_OUT to an r_edc_fullgrid output directory}"
module purge 2>/dev/null || true
module load python/3.11.7
source "${BASE}/venv/bin/activate"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${BASE}/src"
/usr/bin/time -v python hpc/r_reference/linkage_tie_check.py "${R_OUT}" --perms 4 --port-perms 2

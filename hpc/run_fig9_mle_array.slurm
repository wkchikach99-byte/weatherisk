#!/bin/bash
#SBATCH --job-name=fig9_mle
#SBATCH --partition=batch
#SBATCH --array=0-95
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fig9_mle_%A_%a.out
#SBATCH --error=logs/fig9_mle_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wikchi001@awi.de

# ============================================================
#  SLURM array job: Parallelise MLE across latitude rows.
#
#  The T63 grid is 96 lat × 192 lon.  Each array task processes
#  one latitude row (up to 192 grid points).
#
#  This is the embarrassingly-parallel step that dominates runtime.
#  Pre-processing and clustering are done in the non-array wrapper.
#
#  Submit:
#    # First run preprocessing:
#    sbatch hpc/run_fig9_preprocess.slurm
#    # Then (after preprocessing completes):
#    sbatch --dependency=afterok:<JOB_ID> hpc/run_fig9_mle_array.slurm
#    # Then clustering:
#    sbatch --dependency=afterok:<JOB_ID> hpc/run_fig9_cluster.slurm
# ============================================================

WORK_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "${WORK_DIR}"
mkdir -p logs

module purge 2>/dev/null || true
module load python/3.11 2>/dev/null || module load python/3.10 2>/dev/null || \
    module load anaconda3 2>/dev/null || true

if [ -f "${WORK_DIR}/.venv/bin/activate" ]; then
    source "${WORK_DIR}/.venv/bin/activate"
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Array task ${SLURM_ARRAY_TASK_ID}: latitude row processing"

python3 -c "
import sys, os, numpy as np
sys.path.insert(0, '.')

lat_row = ${SLURM_ARRAY_TASK_ID}
data_dir = 'output/cmip6_fig9'
frechet = np.load(os.path.join(data_dir, 'frechet.npy'))
valid_idx = np.load(os.path.join(data_dir, 'valid_idx.npy'))
n_lat = int(np.load(os.path.join(data_dir, 'n_lat.npy')))
n_lon = int(np.load(os.path.join(data_dir, 'n_lon.npy')))

from weatherisk.cmip6_pipeline import _grid_coords, _local_mle_one_cmip6, CMIP6Config

cfg = CMIP6Config()
grid_coords = _grid_coords(valid_idx, n_lat, n_lon)

# Find cells in this latitude row
row_mask = (valid_idx // n_lon) == lat_row
row_indices = np.where(row_mask)[0]

if len(row_indices) == 0:
    print(f'  Row {lat_row}: no valid cells, skipping')
    sys.exit(0)

print(f'  Row {lat_row}: {len(row_indices)} cells')

est = np.zeros((len(row_indices), 3))
for k, cidx in enumerate(row_indices):
    est[k] = _local_mle_one_cmip6(frechet, cidx, grid_coords, cfg)

outfile = os.path.join(data_dir, f'est_row_{lat_row:03d}.npy')
np.save(outfile, est)
np.save(os.path.join(data_dir, f'idx_row_{lat_row:03d}.npy'), row_indices)
print(f'  Saved {outfile}')
"

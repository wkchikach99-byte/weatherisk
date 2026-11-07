#!/bin/bash
#SBATCH --job-name=r_edc_fullgrid
#SBATCH --partition=smp
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_r_reference/logs/%x_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_r_reference/logs/%x_%j.out
###############################################################################
# Criterion F1, R side: the reference R code's deterministic EDC chain on the full
# Figure 9 input (AWI-ESM-1-1-LR historical pr, 1850-2005, 96 x 192 cells),
# clustered in R's cell order and in weatherisk's cell order.
#
# Submit from the deploy directory:
#   cd /albedo/work/user/wikchi001/weatherisk_r_reference && sbatch hpc/r_reference/run_r_edc_fullgrid.slurm
#
# Resources (estimates from local arm64 benchmarks of the same R code):
#   madogram  c_extrcoeff_matrix, about 8 us/pair x 1.70e8 pairs, about 25 min
#             locally, allow 2x on albedo
#   stl       5.6 ms/cell x 18432, about 2 min;  GEV 1.9 ms/cell, about 1 min
#   hclust    about 10 s per order;  reading 165 files with ncdump, a few min
#   memory    the 18432 x 18432 double matrix is 2.7 GB; upper.tri, the
#             permuted copy and hclust's copy bring the peak to about 10-15 GB
# Single threaded (R loops), so one CPU.
###############################################################################
set -euo pipefail

BASE=/albedo/work/user/wikchi001/weatherisk_r_reference
DATA=/albedo/pool/CMIP6/data/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212
OUT=${BASE}/output/edc_fullgrid_${SLURM_JOB_ID}

module purge
module load r/4.2.2 netcdf-c/4.8.1-gcc12.1.0
export R_LIBS_USER=${HOME}/R/libs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

cd "${BASE}"
echo "host $(hostname)  start $(date -Is)  job ${SLURM_JOB_ID}"
which Rscript ncdump
Rscript --version
sha256sum r_code/functions.R hpc/r_reference/r_edc_fullgrid.R
echo "input files: $(ls ${DATA}/*.nc | wc -l)"

/usr/bin/time -v Rscript hpc/r_reference/r_edc_fullgrid.R \
  --functions r_code/functions.R \
  --netcdf-dir "${DATA}" \
  --year-start 1850 --year-end 2005 \
  --quantile 0.3 --min-k 2 \
  --orders r,py \
  --out "${OUT}"

echo "end $(date -Is)"
ls -la "${OUT}"

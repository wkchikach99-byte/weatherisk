#!/bin/bash
#SBATCH --job-name=r_lec_on_estimates
#SBATCH --partition=smp
#SBATCH --account=paleodyn.paleodyn
#SBATCH --qos=12h
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/albedo/work/user/wikchi001/weatherisk_r_reference/logs/%x_%j.out
#SBATCH --error=/albedo/work/user/wikchi001/weatherisk_r_reference/logs/%x_%j.out
###############################################################################
# Criterion F2, R side: R's LEC clustering stage on weatherisk's full-data
# smoothed estimates.  NOT SUBMITTED YET: it needs the estimates file.
#
# Make the input from weatherisk's pipeline_results.npz (see the comment in
# r_lec_on_estimates.R and the report), then:
#   cd /albedo/work/user/wikchi001/weatherisk_r_reference
#   sbatch --export=ALL,EST=/path/smoothed.f64le.bin,N=18432 hpc/r_reference/run_r_lec_on_estimates.slurm
#
# Reports (a) R as-is, (b) R with the validation patch
# validation_patch_calc_distance_ellipses_fixed.R, (c) the number of i<j
# pairs in the empty-ellipse branch, in lec_summary.csv.
#
# Resources: calc_distance_ellipses took 31.6 us/pair locally, so 1.70e8
# pairs take about 1.5 h per variant locally; two variants run in sequence,
# allow 2x on albedo, about 6 h.  Memory: two 2.7 GB matrices plus hclust's
# copy and the quantile copy, about 12 GB peak.
###############################################################################
set -euo pipefail
: "${EST:?set EST to the f64le estimates file}"
: "${N:?set N to the number of cells}"

BASE=/albedo/work/user/wikchi001/weatherisk_r_reference
OUT=${BASE}/output/lec_on_estimates_${SLURM_JOB_ID}

module purge
module load r/4.2.2
export R_LIBS_USER=${HOME}/R/libs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

cd "${BASE}"
echo "host $(hostname)  start $(date -Is)  job ${SLURM_JOB_ID}"
sha256sum r_code/functions.R hpc/r_reference/r_lec_on_estimates.R \
  hpc/r_reference/validation_patch_calc_distance_ellipses_fixed.R "${EST}"

/usr/bin/time -v Rscript hpc/r_reference/r_lec_on_estimates.R \
  --functions r_code/functions.R \
  --estimates-bin "${EST}" --n "${N}" \
  --quantile 0.3 --min-k 2 --variant both \
  --out "${OUT}"

echo "end $(date -Is)"
ls -la "${OUT}"

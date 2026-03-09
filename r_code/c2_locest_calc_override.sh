#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --array=1-51
#SBATCH --qos=12h
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jcontzen@awi.de

module load conda
source activate myEnv
srun Rscript c2_locest_calc.R $SLURM_ARRAY_TASK_ID $1 $2 $3 $4

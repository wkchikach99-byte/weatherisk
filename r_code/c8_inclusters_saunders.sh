#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --array=1-25
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jcontzen@awi.de

module load conda
source activate myEnv
srun Rscript c8_inclusters_saunders.R $1 $2 $3 $SLURM_ARRAY_TASK_ID

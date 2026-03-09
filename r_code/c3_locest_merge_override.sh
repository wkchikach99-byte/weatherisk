#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jcontzen@awi.de
module load conda
source activate myEnv
srun Rscript c3_locest_merge.R $1 $2 $3 $4
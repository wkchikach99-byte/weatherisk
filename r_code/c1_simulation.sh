#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jcontzen@awi.de
module load conda
source activate myEnv
srun Rscript c1_simulation.R $1
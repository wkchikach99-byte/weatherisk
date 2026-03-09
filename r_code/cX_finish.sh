#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --mail-type=END	
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jcontzen@awi.de
#SBATCH --mem=30G
module load conda
source activate myEnv
srun Rscript cX_finish.R $1
#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 24
#SBATCH --time=2159

srun python vlinc.l1.analysis.py



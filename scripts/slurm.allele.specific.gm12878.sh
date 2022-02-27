#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 1
#SBATCH --time=400

srun allele.specific.analysis.gm12878.sh $1 $2

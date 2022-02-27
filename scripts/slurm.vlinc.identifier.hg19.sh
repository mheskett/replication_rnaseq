#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 1
#SBATCH --time=2159

srun vlinc.identifier.hg19.sh $1 $2 $3 $4 $5

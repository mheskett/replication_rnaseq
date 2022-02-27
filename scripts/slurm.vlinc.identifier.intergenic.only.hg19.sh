#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 1
#SBATCH --time=10:0:0

srun vlinc.identifier.intergenic.only.hg19.sh $1 $2 $3 $4 $5

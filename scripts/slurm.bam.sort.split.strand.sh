#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 8
#SBATCH --time=2159

srun bam.sort.split.strand.sh $1 $2

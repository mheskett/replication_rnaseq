#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 8
#SBATCH --time=2159

srun samtools.process.bam.sh $1 $2

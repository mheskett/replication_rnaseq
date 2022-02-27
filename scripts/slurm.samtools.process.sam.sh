#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 60000
#SBATCH -c 8
#SBATCH --time=1000
#SBATCH -A SpellmanLab

srun samtools.process.sam.sh $1 $2 

#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 8
#SBATCH --time=2159

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/samtools.process.bam.sh $1 $2


#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 1
#SBATCH --time=2159

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/samtools.merge.sh $1


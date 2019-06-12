#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 1
#SBATCH --time=2159

srun /home/groups/Spellmandata/heskett/replication.rnaseq/old.scripts/process.bams.gatk3.no.variant.calling.sh $1 $2


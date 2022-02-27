#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 30000
#SBATCH -c 1
#SBATCH --time=1000
#SBATCH -A SpellmanLab

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/repliseq.analysis.sh $1 $2 $3

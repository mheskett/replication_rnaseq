#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 16000
#SBATCH -c 1
#SBATCH --time=2159
#SBATCH -A SpellmanLab

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/repliseq.analysis.sh $1 $2 $3

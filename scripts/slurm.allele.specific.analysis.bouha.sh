#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 30000
#SBATCH -c 1
#SBATCH --time=6:0:0
#SBATCH -A SpellmanLab

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/allele.specific.analysis.bouha.sh $1 $2

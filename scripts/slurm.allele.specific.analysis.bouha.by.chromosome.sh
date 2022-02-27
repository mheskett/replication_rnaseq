#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 6000
#SBATCH -c 1
#SBATCH --time=5:0:0
#SBATCH -A SpellmanLab

srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/allele.specific.analysis.bouha.by.chromosome.sh $1 $2 $3 $4 $5

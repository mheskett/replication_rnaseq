#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem 30000
#SBATCH -c 1
#SBATCH --time=10:0:0
#SBATCH -A SpellmanLab

newgrp SpellmanLab
srun /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/allele.specific.repliseq.sh $1 $2

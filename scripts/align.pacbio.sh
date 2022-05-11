#!/bin/bash

#SBATCH --time 1000
#SBATCH -c 16
#SBATCH --mem 30000

srun blasr $1 ../../../myron_refs/human_g1k_v37.fasta --sam --out $2 --nproc 16

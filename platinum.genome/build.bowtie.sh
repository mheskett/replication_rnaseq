#!/bin/bash

#SBATCH --mem 20000
#SBATCH --time 1000
#SBATCH -c 6

srun bowtie2-build --threads 6 bcftools.na12878.hap1.fa na12878.hap1
srun bowtie2-build --threads 6 bcftools.na12878.hap2.fa na12878.hap2

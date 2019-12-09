#!/bin/bash


#SBATCH --partition=exacloud
#SBATCH --mem=50000
#SBATCH -c 8
#SBATCH --time=2159
#SBATCH -A SpellmanLab

b=$(basename $1) ## removes /path/to/file
filename=${b%.*} ### removes file extension

out_dir=$2

srun samtools sort -@ 8 -o $outdir$filename.bam -m 5G $1 
srun samtools index $outdir$filename.bam

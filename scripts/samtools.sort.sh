#!/bin/bash


#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 8
#SBATCH --time=2159
#SBATCH -A SpellmanLab

b=$(basename $1) ## removes /path/to/file
filename=${b%.*} ### removes file extension

out_dir=$2

srun samtools sort -@ 8 -o $outdir$filename.sorted.bam -m 5G $1 
srun samtools index $outdir$filename.sorted.bam

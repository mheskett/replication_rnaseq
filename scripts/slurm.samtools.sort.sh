#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 4
#SBATCH --time=2159

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

srun samtools sort -@ 4 -o $out_dir$filename.sorted.bam $input

#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH --time=2159

srun gatk DownsampleSam -I=$1 -O=$2 -P=0.5 -M=downsample.metrics.txt -S=Chained -ACCURACY=0.01

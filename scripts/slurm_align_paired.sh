#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=50000
#SBATCH -c 16
#SBATCH --time=2159

fq1=$1
fq2=$2
out_pref=$3

srun STAR \
  --runThreadN 16 \
  --twopassMode Basic \
  --readFilesIn $fq1 $fq2 \
  --outFileNamePrefix $out_pref \
  --genomeDir /home/groups/Spellmandata/heskett/tools/refdata-cellranger-GRCh38-3.0.0/star/ \
  --outSAMmapqUnique 60 \
  --outSAMtype BAM Unsorted \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI NM MD jM jI \
  --readFilesCommand zcat
#  --limitBAMsortRAM 53130127667

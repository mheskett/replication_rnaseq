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
  --runMode genomeGenerate \
  --genomeDir /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/star.hg19 \
  --genomeFastaFiles /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/hg19.fa \
  --sjdbGTFfile /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.genes.hg19.gtf \
  --sjdbOverhang 100

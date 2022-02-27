#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem=60000
#SBATCH -c 2
#SBATCH --time=2159

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

srun gatk HaplotypeCaller --java-options "-Xmx60G" \
  -I $1 \
  --dbsnp /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.sorted.vcf \
  -O $out_dir/$filename.snps.vcf \
  -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa \
  -stand-call-conf 20

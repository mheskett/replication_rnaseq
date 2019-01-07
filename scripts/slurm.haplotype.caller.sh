#!/bin/bash

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

#SBATCH --partition=exacloud
#SBATCH --mem=30000
#SBATCH -c 2
#SBATCH --time=2159

gatk HaplotypeCaller --java-options "-Xmx30G" \
  -I $1 \
  --dbsnp /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.vcf \
  -O $out_dir/$filename.snps.vcf \
  -R /home/groups/Spellmandata/heskett/tools/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
  -L /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.vcf \
  -ip 100 \
  -stand_call_conf 30

#!/bin/bash

## start from sorted, index and rmduped BAM file.

input=$1
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
out_dir=$2

## get allele specific counts at heterozygous sites
## strand agnostic bc this is DNA-seq

bcftools mpileup -R /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/na12878.hg19.het.bed \
  -f /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta \
  -a DP,AD \ 
  $1 > $out_dir$filename.allele.counts.vcf

gatk VariantsToTable -V $out_dir$filename.allele.counts.vcf -O $out_dir$filename.table -F CHROM -F POS -F REF -F ALT -GF DP -GF AD

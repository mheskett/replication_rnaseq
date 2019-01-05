#!/bin/bash

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
random=$RANDOM

gatk AddOrReplaceReadGroups --java-options "-Xmx30G" -I=$input -O=$out_dir/$filename.rg.sorted.bam \
  --RGID=$random --RGLB=lib1 --RGPL=illumina --RGPU=unit1 --RGSM=$filename -SO=coordinate

gatk MarkDuplicates --java-options "-Xmx30G" -I=$out_dir/$filename.rg.sorted.bam -O=$out_dir/$filename.rg.sorted.markdup.bam \
  --METRICS_FILE=$filename.metrics.txt --ASSUME_SORT_ORDER=coordinate --REMOVE_DUPLICATES=true

# make sure VCF and bam are sorted in the same way with matching contigs etc
/home/exacloud/lustre1/SpellmanLab/heskett/facets/inst/extcode/snp-pileup -q30 -Q30 -r2 /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.vcf \
  $out_dir/$filename.snp.counts.txt $out_dir/$filename.rg.sorted.bam



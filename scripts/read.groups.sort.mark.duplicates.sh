#!/bin/bash

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
random=$RANDOM

### with STAR, used base max map quality score of 60 instead of default 255
### fix misencoded scores to convert from illumina 1.5 to 1.9
### could add base recalibration

gatk AddOrReplaceReadGroups --java-options "-Xmx30G" -I=$input -O=$out_dir/$filename.rg.sorted.bam \
  --RGID=1 --RGLB=lib1 --RGPL=illumina --RGPU=unit1 --RGSM=$filename -SO=coordinate

gatk FixMisencodedBaseQualityReads -I $filename.rg.sorted.bam -O $out_dir/$filename.rg.sorted.fix.bam

gatk MarkDuplicates --java-options "-Xmx30G" -I=$out_dir/$filename.rg.sorted.fix.bam -O=$out_dir/$filename.rg.sorted.markdup.bam \
  --METRICS_FILE=$filename.metrics.txt --ASSUME_SORT_ORDER=coordinate --REMOVE_DUPLICATES=true --CREATE_INDEX=true

gatk SplitNCigarReads --java-options "-Xmx30G" -I $out_dir/$filename.rg.sorted.markdup.bam -O $out_dir/$filename.rg.sorted.markdup.splitn.bam \
  -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa 

# make sure VCF and bam are sorted in the same way with matching contigs etc

#/home/exacloud/lustre1/SpellmanLab/heskett/facets/inst/extcode/snp-pileup -q30 -Q30 -r2 /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.vcf \
#  $out_dir/$filename.snp.counts.txt $out_dir/$filename.rg.sorted.markdup.splitn.bam

rm $out_dir/$filename.rg.sorted.bam
rm $out_dir/filename.rg.sorted.markdup.bam


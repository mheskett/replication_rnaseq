#!/bin/bash

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
random=$RANDOM

#gatk AddOrReplaceReadGroups --java-options "-Xmx30G" -I=$input -O=$out_dir/$filename.rg.sorted.bam \
#  --RGID=$random --RGLB=lib1 --RGPL=illumina --RGPU=unit1 --RGSM=$filename
#VALIDATION_STRINGENCY=silent SORT_ORDER=coordinate

gatk MarkDuplicates --java-options "-Xmx30G" -I=$out_dir/$filename.rg.sorted.bam -O=$out_dir/$filename.rg.sorted.markdup.bam \
  --METRICS_FILE=$filename.metrics.txt --ASSUME_SORT_ORDER=coordinate --REMOVE_DUPLICATES=true


/home/exacloud/lustre1/SpellmanLab/heskett/facets/inst/extcode/snp-pileup -q30 -Q30 -r0 snpfile outputfile sequence file



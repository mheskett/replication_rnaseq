#!/usr/bin/bash

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

##########
gatk AddOrReplaceReadGroups --java-options "-Xmx30G" -I=$input -O=$out_dir/$filename.rg.sorted.bam \
  --RGID=1 --RGLB=lib1 --RGPL=illumina --RGPU=unit1 --RGSM=$filename -SO=coordinate

##########
gatk MarkDuplicates --java-options "-Xmx30G" -I=$out_dir/$filename.rg.sorted.bam -O=$out_dir/$filename.rg.sorted.markdup.bam \
  --METRICS_FILE=$filename.metrics.txt --ASSUME_SORT_ORDER=coordinate --REMOVE_DUPLICATES=true --CREATE_INDEX=true

### for files with correct illumina encoding
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
-T SplitNCigarReads -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa -I $out_dir/$filename.rg.sorted.markdup.bam \
  -o $out_dir/$filename.split.bam \
  -U ALLOW_N_CIGAR_READS

#!/usr/bin/bash

input=$1 ## BAM file from BWA, or whatever.
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## remove reads with less than 20 map quality
## sort
## remove duplicates
samtools view -@ 8 -bhq 20 $input > $out_dir$filename.mq30.bam
samtools sort -@ 8 -m 4G $out_dir$filename.mq30.bam > $out_dir$filename.mq30.sorted.bam

# paired end default
samtools rmdup $out_dir$filename.mq30.sorted.bam $out_dir$filename.samtool.rmdup.bam 

## index bam
samtools index -b $out_dir$filename.samtool.rmdup.bam

rm $out_dir$filename.mq30.bam
rm $out_dir$filename.mq30.sorted.bam

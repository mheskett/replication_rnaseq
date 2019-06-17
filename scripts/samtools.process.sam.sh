#!/usr/bin/bash

input=$1 ## sam file from BWA, or whatever.
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## sam to bam. remove reads with less than 20 map quality. this may require a LOT of memory...
samtools view -@ 8 -bSq 20 $input > $out_dir$filename.bam

samtools sort -@ 8 -m 6G -o $out_dir$filename.sorted.bam $out_dir$filename.bam 

samtools rmdup -S $out_dir$filename.sorted.bam $out_dir$filename.sorted.markdup.bam

## index bam
samtools index -b $out_dir$filename.sorted.markdup.bam

rm $out_dir$filename.bam
rm $out_dir$filename.sorted.bam



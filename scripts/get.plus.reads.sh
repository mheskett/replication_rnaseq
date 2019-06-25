#!/bin/bash

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## get plus strand
samtools view -b -f 131 -F 16 $input > $out_dir$filename.fwd1.bam
samtools index $out_dir$filename.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.fwd2.bam

# now combine, this should contail all plus strand genes
samtools merge -f $out_dir$filename.plus.bam $out_dir$filename.fwd1.bam $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.plus.bam

rm $out_dir$filename.fwd1.bam
rm $out_dir$filename.fwd2.bam

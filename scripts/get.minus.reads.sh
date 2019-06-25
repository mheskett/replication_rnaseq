#!/bin/bash

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## get minus strand
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam

rm $out_dir$filename.rev1.bam
rm $out_dir$filename.rev2.bam

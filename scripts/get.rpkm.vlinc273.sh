#!/bin/bash

input=$1 ## bam file
out_dir=$2
original_bam=$3
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## this considers the original library size as mapped reads, but not deduplicated or quality score filtered
library_size=$(samtools view -c -F 260 $original_bam)
echo $library_size

## get minus strand
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam

counts=$(samtools view -c $out_dir$filename.minus.bam 6:141034658-141219628)

## 184.97kb is vlinc273 length
python -c "print( ${counts} / 184.97 / (${library_size}/1000000) )" > $out_dir$filename.vlinc273.rpkm.txt

#!/bin/bash

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## include header and filter out map quality below 40
samtools view -h -b -q 40 $input > $out_dir$filename.mq40.bed

## get the per base genome coverage, while NOT including inferred coverage of spliced out introns 
bedtools genomecov -bga -split -strand + -ibam $out_dir$filename.mq40.bed  > $out_dir$filename.mq40.plus.cov.bed
bedtools genomecov -bga -split -strand - -ibam $out_dir$filename.mq40.bed  > $out_dir$filename.mq40.minus.cov.bed


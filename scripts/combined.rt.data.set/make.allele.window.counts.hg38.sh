#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

#bedtools makewindows -g hg38.fa.fai -w 250000 -s 250000 | \
#  bedtools subtract -a stdin -b hg38-blacklist.v2.bed | \
#  bedtools subtract -a stdin -b hg38.imprinted.genes.txt > hg38.windows.w250.s250.bed

bedtools makewindows -g hg38.fa.fai -w 250000 -s 250000 | \
  bedtools subtract -a stdin -b hg38-low.mappability.v2.bed > hg38.windows.w250.s250.bed

###
bedtools sort -g hg38.fa.fai -i $1 | 
              bedtools map -a hg38.windows.w250.s250.bed \
              -b stdin -o sum,sum -c 6,7  > $filename.windows.w250.s250.bed

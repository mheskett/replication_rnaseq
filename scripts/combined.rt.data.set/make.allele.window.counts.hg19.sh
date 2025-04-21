#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

bedtools makewindows -g Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai -w 250000 -s 250000 | \
  bedtools subtract -a stdin -b hg19-blacklist.v2.nochr.bed | \
  bedtools subtract -a stdin -b hg19.imprinted.genes.txt > hg19.windows.w250.s250.bed

bedtools sort -g Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai -i $1 | 
              bedtools map -a hg19.windows.w250.s250.bed \
              -b stdin -o sum,sum -c 6,7  > $filename.windows.w250.s250.bed

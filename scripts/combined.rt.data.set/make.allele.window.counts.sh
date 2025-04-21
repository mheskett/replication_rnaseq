#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## hg19.ethan.bedtools.genome.lexicographically2.w250.bed deprecated dont use ethan ref
## human_g1k_v37.w250.s250.bed non overlaping
## human_g1k_v37.w250.s50.bed overlapping then can smooth out (requires CPMs prior to this step)
## for sliding windows
## cant use grep v \. because now you have decimals

bedtools sort -g human_g1k_v37.fasta.fai -i $1 | 
              bedtools map -a human_g1k_v37.w250.s125.bed \
              -b stdin -o sum,sum -c 6,7  > $filename.allele.counts.windows.s125.bed

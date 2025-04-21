#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension


bedtools map -a hg19.ethan.bedtools.genome.lexicographically2.w250.s50.bed \
             -b $1 -o sum,sum -c 6,7 > $filename.allele.counts.windows.bed

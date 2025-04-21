#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

# do this BEFORE library size normalization

bedtools subtract -a $input -b hg19-blacklist.v2.nochr.bed > $filename.rmv.blck.bed

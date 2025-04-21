#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

# do this BEFORE library size normalization

# Should probably do this
#bedtools subtract -a $input -b hg38-blacklist.v2.bed > $filename.rmv.blck.bed

# or do nothing
cat $input > $filename.rmv.blck.bed

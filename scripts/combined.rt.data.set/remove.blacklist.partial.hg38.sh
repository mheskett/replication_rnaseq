#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

# remove partial windows after 250kb windows. do this after normalization and windows
bedtools subtract  -a $input -b hg38-blacklist.v2.bed > $filename.rmv.blck.bed

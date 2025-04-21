#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

# Should probably do this
bedtools subtract -a $input -b hg38.imprinted.genes.txt > $filename.rmv.imprnt.bed

#!/bin/bash

source activate for_bedtools

input=$1 ## bam file
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

# Should probably do this
#bedtools subtract -a $input -b hg38-blacklist.v2.bed | bedtools subtract -a stdin -b hg38.imprinted.genes.txt > $filename.rmv.blck.bed

bedtools subtract -a $input -b hg38-low.mappability.v2.bed > $filename.rmv.blck.bed

# just remove low mappability
#cat $input > $filename.rmv.blck.bed

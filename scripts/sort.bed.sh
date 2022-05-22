#!/bin/bash

#b=$(basename $1) ## removes path to file
b=$1
filename=${b%.*} ### removes file extension


bedtools sort -g /Users/mike/replication_rnaseq/human_g1k_v37.fasta.fai -i $1 > $filename.sorted.bed

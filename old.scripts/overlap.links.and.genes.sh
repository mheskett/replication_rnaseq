#!/bin/bash

input=$1 
out_dir=$2 # include trailing slash
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

bedtools intersect -a $input -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/genes.ucsc.nochr.hg38.bed \
  -wa -wb > $out_dir$filename.annotated.genes.bed
bedtools intersect -a $input -b -b /home/groups/Spellmandata/heskett/replication.rnaseq/data.from.vlinc.paper/kaprakov.vlinc1541.nochr.hg38.sorted.final.bed \
  -wa -wb > $out_dir$filename.annotated.links.bed


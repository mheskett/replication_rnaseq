#!/bin/bash

## OK this script will get all the reads within pre determined windows
## for example vlinc RNAs. Needs to be proper bed format, and strand
## specific is desired.
## Since allele specific counting is required
## If strand is not known, get plus and minus reads at that position

input=$1 ## bam file
out_dir=$2 
windows=$3 # where to count reads. proper bed file with strands
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

## get plus strand
samtools view -b -f 131 -F 16 $input > $out_dir$filename.fwd1.bam
samtools index $out_dir$filename.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.fwd2.bam

# now combine, this should contail all plus strand genes
samtools merge -f $out_dir$filename.plus.bam $out_dir$filename.fwd1.bam $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.plus.bam

## get minus strand
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam
 
 

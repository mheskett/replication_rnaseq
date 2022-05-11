#!/bin/bash

# Separate reads in an RNA-seq alignment based on plus and minus strand. 


input=$1 # Deduplicated BAM file
b=$(basename $input) # removes path to file
name=${b%.*} #removes file extension
out_dir=$2


echo "Input file: " $input 
echo "Sample name: " $name

# Get second in pair reads mapping to forward
# get first in pair reads mapping to reverse
# this should return plus stranded genes
# 131 flag is second in pair, -F 16 excludes reverse strand

echo "Get forward reads"

samtools view -b -f 131 -F 16 $input > $out_dir$name.fwd1.bam
samtools index $out_dir$name.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$name.fwd2.bam
samtools index $out_dir$name.fwd2.bam

# now combine, this should contail all plus strand genes

samtools merge -f $out_dir$name.plus.bam $out_dir$name.fwd1.bam $out_dir$name.fwd2.bam
samtools index $out_dir$name.plus.bam

# repeat for minus strand genes

echo "Get reverse reads"
samtools view -b -f 147 $input > $out_dir$name.rev1.bam
samtools index $out_dir$name.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$name.rev2.bam
samtools index $out_dir$name.rev2.bam

samtools merge -f $out_dir$name.minus.bam $out_dir$name.rev1.bam $out_dir$name.rev2.bam
samtools index $out_dir$name.minus.bam 

rm $name.fwd1.bam*
rm $name.fwd2.bam*
rm $name.rev1.bam*
rm $name.rev2.bam*


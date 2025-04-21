#!/bin/bash

plus_input=$1
minus_input=$2
plus_tls=$3
minus_tls=$4
sample_name=$5

#  bedtools must be command line


a=$(basename $plus_input) ## removes /path/to/file
filename_plus=${a%.*} ### removes file extension

b=$(basename $minus_input) ## removes /path/to/file
filename_minus=${b%.*} ### removes file extension

c=$(basename $plus_tls)
filename_tls=${c%.*}

d=$(basename $minus_tls)
filename_tls=${d%.*}

source activate for_bedtools
###### plus
bedtools sort -i $plus_input -g ../human_g1k_v37.sorted.fasta.fai | bedtools map -a $plus_tls \
             -b stdin \
             -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"plus"}' > $filename_plus.tls.counts.plus.bed
#### minus
bedtools sort -i $minus_input -g ../human_g1k_v37.sorted.fasta.fai | bedtools map -a $minus_tls \
             -b stdin \
       	     -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"minus"}' > $filename_minus.tls.counts.minus.bed

### now cat everything together and sort if you want, and then go to python for "add.stats.genes.py"

cat $filename_plus.tls.counts.plus.bed $filename_minus.tls.counts.minus.bed | sort -k1,1 -k2,2n > $sample_name.as.tl.counts.bed


#!/bin/bash

plus_input=$1
minus_input=$2
total_input=$3
#out_dir=$4

#  bedtools must be command line


b=$(basename $plus_input) ## removes /path/to/file
filename_plus=${b%.*} ### removes file extension

c=$(basename $minus_input) ## removes /path/to/file
filename_minus=${c%.*} ### removes file extension

d=$(basename $total_input) ## removes /path/to/file
filename_total=${d%.*} ### removes file extension

## input is the haplotype resolved allele counts bed file
source activate for_bedtools
###### plus
bedtools sort -i $plus_input -g ../hg38.fa.fai | bedtools map -a ../ucsc.refseq.hg38.txn.whole.gene.plus.bed \
             -b stdin \
             -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"plus"}' > $filename_plus.gene.counts.plus.bed
#### minus
bedtools sort -i $minus_input -g ../hg38.fa.fai | bedtools map -a ../ucsc.refseq.hg38.txn.whole.gene.minus.bed \
             -b stdin \
       	     -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"minus"}' > $filename_minus.gene.counts.minus.bed
#### total
bedtools sort -i $total_input -g ../hg38.fa.fai | bedtools map -a ../ucsc.refseq.hg38.txn.whole.gene.bed \
             -b stdin \
             -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"total"}' > $filename_minus.gene.counts.total.bed

### plus genes with minus RNA antisense
bedtools sort -i $minus_input -g ../hg38.fa.fai | bedtools map -a ../ucsc.refseq.hg38.txn.whole.gene.plus.bed \
             -b stdin \
             -o sum,sum -c 6,7 | awk '$7!="."{print $0}'| awk 'OFS="\t" {print $0,"antisense"}' > $filename_minus.gene.counts.plus.minus.antisense.bed

### minus genes with plus RNA antisense
bedtools sort -i $plus_input -g ../hg38.fa.fai | bedtools map -a ../ucsc.refseq.hg38.txn.whole.gene.minus.bed \
             -b stdin \
             -o sum,sum -c 6,7 | awk '$7!="."{print $0}' | awk 'OFS="\t" {print $0,"antisense"}' > $filename_minus.gene.counts.minus.plus.antisense.bed

### now cat everything together and sort if you want, and then go to python for "add.stats.genes.py"

cat $filename_plus.gene.counts.plus.bed $filename_minus.gene.counts.minus.bed $filename_minus.gene.counts.total.bed $filename_minus.gene.counts.plus.minus.antisense.bed $filename_minus.gene.counts.minus.plus.antisense.bed | sort -k1,1 -k2,2n > $filename_total.comprehensive.gene.counts.bed








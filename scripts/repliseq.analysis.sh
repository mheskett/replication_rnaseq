#!/bin/bash

## start from sorted, index and rmduped BAM file.

early=$1
late=$2
e=$(basename $early) ## removes /path/to/file
early_name=${e%.*} ### removes file extension

l=$(basename $late) ## removes /path/to/file
late_name=${l%.*} ### removes file extension

out_dir=$3

## get genome coverage
## this could be used for a simple replication timing early vs late comparison.
## then take the sum across genome windows (non overlapping)
## then just simply take the log2 ratio of early/late

bedtools makewindows -w 50000 -s 50000 -g /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta.fai > $out_dir'human.g1k.v37.50kb.windows.bed'
bedtools coverage -sorted -g /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta.fai -a $out_dir'human.g1k.v37.50kb.windows.bed' -b $early -counts -F 0.51 > $out_dir$early_name.counts.bed 
bedtools coverage -sorted -g /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta.fai -a $out_dir'human.g1k.v37.50kb.windows.bed' -b $late -counts -F 0.51 > $out_dir$late_name.counts.bed

# use python for plotting

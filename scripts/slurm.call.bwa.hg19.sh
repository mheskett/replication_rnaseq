#!/bin/bash

#SBATCH --partition=exacloud
#SBATCH --mem-per-cpu 12G
#SBATCH -c 8
#SBATCH --time=10:0:0

fq1=$1
fq2=$2
name=$3
ref='/home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta'
out_dir=$4
readgroup='@RG\tID:'$name'\tSM:'1'\tPL:'illumina'\tLB:'$name

source activate for_bwa
bwa mem -R $readgroup -t 8 -Ma $ref $fq1 $fq2 > $out_dir$name.sam

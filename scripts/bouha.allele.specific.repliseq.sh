#!/bin/bash

## start from sorted, index and rmduped BAM file.

input=$1
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
out_dir=$2

##
samtools view -@ 8 -bhq 20 $input > $out_dir$filename.mq20.bam
samtools rmdup $out_dir$filename.mq20.bam $out_dir$filename.samtool.rmdup.bam

## index bam
samtools index -b $out_dir$filename.samtool.rmdup.bam


# replace bed file here with any genotype file
# reference filer here is hg19 with no "chr" prefixes
bcftools mpileup -R /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/na12878.hg19.het.bed \
  -f /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta \
  -a DP,AD \
  -q 30 \
  $out_dir$filename.samtool.rmdup.bam > $out_dir$filename.allele.counts.vcf

# get allele specific counts of all het sites into a table
gatk VariantsToTable -V $out_dir$filename.allele.counts.vcf -O $out_dir$filename.table -F CHROM -F POS -F REF -F ALT -GF AD

## turn this table into a bed file so that it can be intersected with the haplotype resolved file
tail -n +2 $out_dir$filename.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' |
  bedtools intersect -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/na12878.hg19.het.bed -wa -wb > $out_dir$filename.allele.counts.bed

## plug this into python script to filter and align the haplotypes to get final file


rm $out_dir$filename.mq20.bam

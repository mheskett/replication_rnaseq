#!/bin/bash

## start from sorted, index and rmduped BAM file.

input=$1
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
out_dir=$2
chrom_vcf_file=$5
chrom=$4
chrom_bed_file=$3

echo 'running'
echo $chrom_file
echo $chrom
echo 'printed'
##
#samtools view -@ 8 -bhq 20 $input > $out_dir$filename.mq20.bam
#samtools rmdup $out_dir$filename.mq20.bam $out_dir$filename.samtool.rmdup.bam

## index bam
#samtools index -b $out_dir$filename.samtool.rmdup.bam


# replace bed file here with any genotype file
# reference filer here is hg19 with no "chr" prefixes
# use a file with format CHR,POS here
bcftools mpileup -R $chrom_vcf_file \
  -f /home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta \
  -a DP,AD \
  -q 20 \
  $1 > $out_dir$filename.$chrom.allele.counts.vcf

# get allele specific counts of all het sites into a table
gatk VariantsToTable -V $out_dir$filename.$chrom.allele.counts.vcf -O $out_dir$filename.$chrom.table -F CHROM -F POS -F REF -F ALT -GF AD

## turn this table into a bed file so that it can be intersected with the haplotype resolved file

## recent update. for bouhassira data I am making this awk statement 2,2+1. it is possible that 
## the output of bcftools mpileup changed? or actually I think this "homemade" VCF from bouha lab is one position off
## that sounds like the right answer..... so for a regular VCF it is 2-1,2
##3 use a file with format CHR,START,STOP here
tail -n +2 $out_dir$filename.$chrom.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' |
  bedtools intersect -a stdin -b $chrom_bed_file -wa -wb > $out_dir$filename.$chrom.allele.counts.bed

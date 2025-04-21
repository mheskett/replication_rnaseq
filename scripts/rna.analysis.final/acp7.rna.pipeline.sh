#!/bin/bash
###
python align.haplotypes.py --bed acp7_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c4_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c4_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp7_c4_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./

## remove blacklist
./remove.blacklist.hg38.sh acp7_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c4_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c4_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp7_c4_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed

## get gene counts
./get.comprehensive.gene.counts.sh acp7_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh acp7_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh acp7_c4_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c4_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
acp7_c4_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

## cat acp7 tls
cat acp7_c1.tls.minus.bed acp7_c1.tls.plus.bed > acp7_c1.tls.bed
cat acp7_c2.tls.minus.bed acp7_c2.tls.plus.bed > acp7_c2.tls.bed
cat acp7_c4.tls.minus.bed acp7_c4.tls.plus.bed > acp7_c4.tls.bed

# get TL allele specific counts

sort -k1,1 -k2,2n acp7.all.clones.plus.tls.bed > acp7.all.clones.plus.tls.sorted.bed
sort -k1,1 -k2,2n acp7.all.clones.minus.tls.bed > acp7.all.clones.minus.tls.sorted.bed

./get.tl.counts.sh acp7_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7.all.clones.plus.tls.sorted.bed \
  acp7.all.clones.minus.tls.sorted.bed \
  acp7_c1

./get.tl.counts.sh acp7_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7.all.clones.plus.tls.sorted.bed \
  acp7.all.clones.minus.tls.sorted.bed \
  acp7_c2

./get.tl.counts.sh acp7_c4_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7_c4_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp7.all.clones.plus.tls.sorted.bed \
  acp7.all.clones.minus.tls.sorted.bed \
  acp7_c4



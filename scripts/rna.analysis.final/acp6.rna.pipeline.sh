#!/bin/bash
###
python align.haplotypes.py --bed acp6_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c5_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c5_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c5_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c6_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c6_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.bed --out_directory ./
python align.haplotypes.py --bed acp6_c6_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.bed --out_directory ./

# remove black 
./remove.blacklist.hg38.sh acp6_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed 
./remove.blacklist.hg38.sh acp6_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c5_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c5_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c5_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c6_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c6_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c6_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.bed

## get gene counts
./get.comprehensive.gene.counts.sh acp6_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c1_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh acp6_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c2_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh acp6_c5_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c5_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c5_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh acp6_c6_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c6_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c6_rnaAligned.out.samtool.rmdup.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed

## cat TL files together
cat acp6_c1.tls.minus.bed acp6_c1.tls.plus.bed > acp6_c1.tls.bed
cat acp6_c2.tls.minus.bed acp6_c2.tls.plus.bed > acp6_c2.tls.bed
cat acp6_c5.tls.minus.bed acp6_c5.tls.plus.bed > acp6_c5.tls.bed
cat acp6_c6.tls.minus.bed acp6_c6.tls.plus.bed > acp6_c6.tls.bed

# get Allele specific TL counts
sort -k1,1 -k2,2n acp6.all.clones.plus.tls.bed > acp6.all.clones.plus.tls.sorted.bed
sort -k1,1 -k2,2n acp6.all.clones.minus.tls.bed  > acp6.all.clones.minus.tls.sorted.bed


./get.tl.counts.sh acp6_c1_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \ 
  acp6_c1_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6.all.clones.plus.tls.sorted.bed \
  acp6.all.clones.minus.tls.sorted.bed \
  acp6_c1

./get.tl.counts.sh acp6_c2_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c2_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6.all.clones.plus.tls.sorted.bed \
  acp6.all.clones.minus.tls.sorted.bed \
  acp6_c2

./get.tl.counts.sh acp6_c5_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c5_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6.all.clones.plus.tls.sorted.bed \
  acp6.all.clones.minus.tls.sorted.bed \
  acp6_c5

./get.tl.counts.sh acp6_c6_rnaAligned.out.samtool.rmdup.plus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6_c6_rnaAligned.out.samtool.rmdup.minus.allele.counts.all.chrom.haplotype.resolved.counts.rmv.blck.bed \
  acp6.all.clones.plus.tls.sorted.bed \
  acp6.all.clones.minus.tls.sorted.bed \
  acp6_c6

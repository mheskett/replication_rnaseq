#!/bin/bash
###
python align.haplotypes.py --bed gm12878_clone4_rnaAligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone4_rnaAligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone4_rnaAligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone5_rnaAligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone5_rnaAligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone5_rnaAligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./


## remove black
./remove.blacklist.hg38.sh gm12878_clone4_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone4_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone4_rnaAligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone5_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone5_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone5_rnaAligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed

## call get gene counts
./get.comprehensive.gene.counts.sh gm12878_clone4_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm12878_clone4_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm12878_clone4_rnaAligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.sh gm12878_clone5_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm12878_clone5_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm12878_clone5_rnaAligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

## cat TLs together
cat gm_clone4.tls.minus.bed gm_clone4.tls.plus.bed > gm_clone4.tls.bed
cat gm_clone5.tls.minus.bed gm_clone5.tls.plus.bed > gm_clone5.tls.bed

## get allele specific tls
sort -k1,1 -k2,2n gm.all.clones.plus.tls.bed > gm.all.clones.plus.tls.sorted.bed
sort -k1,1 -k2,2n gm.all.clones.minus.tls.bed > gm.all.clones.minus.tls.sorted.bed


gm.all.clones.minus.tls.bed
./get.tl.counts.sh gm12878_clone4_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm12878_clone4_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  gm.all.clones.plus.tls.sorted.bed \
  gm.all.clones.minus.tls.sorted.bed \
  gm_clone4

./get.tl.counts.sh gm12878_clone5_rnaAligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \  
  gm12878_clone5_rnaAligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \   
  gm.all.clones.plus.tls.sorted.bed \
  gm.all.clones.minus.tls.sorted.bed \
  gm_clone5



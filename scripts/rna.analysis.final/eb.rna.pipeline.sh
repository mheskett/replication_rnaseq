#!/bin/bash
## align haplotypes
python align.haplotypes.py --bed eb3_2_clone10Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./ 
python align.haplotypes.py --bed eb3_2_clone13Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone15Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone2Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone3Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone4Aligned.out.samtool.rmdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone10Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone13Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone15Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone2Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone3Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone4Aligned.out.samtool.rmdup.plus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone10Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone13Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone15Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone2Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone3Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed eb3_2_clone4Aligned.out.samtool.rmdup.minus.allele.counts.bed --out_directory ./

# remove blacklist hg19
./remove.blacklist.hg19.sh eb3_2_clone10Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone13Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone15Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone2Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone3Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone4Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone10Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone13Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone15Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone2Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone3Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone4Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone10Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone13Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone15Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone2Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone3Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg19.sh eb3_2_clone4Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.bed


# get gene counts

./get.comprehensive.gene.counts.hg19.sh eb3_2_clone10Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone10Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone10Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.hg19.sh eb3_2_clone13Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone13Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone13Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed \

./get.comprehensive.gene.counts.hg19.sh eb3_2_clone15Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone15Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone15Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
  
./get.comprehensive.gene.counts.hg19.sh eb3_2_clone2Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone2Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone2Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.hg19.sh eb3_2_clone3Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone3Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone3Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

./get.comprehensive.gene.counts.hg19.sh eb3_2_clone4Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone4Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone4Aligned.out.samtool.rmdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

### get gene counts from TLs

sort -k1,1 -k2,2n eb32.all.clones.tls.plus.bed > eb32.all.clones.tls.plus.sorted.bed
sort -k1,1 -k2,2n eb32.all.clones.tls.minus.bed > eb32.all.clones.tls.minus.sorted.bed

./get.tl.counts.hg19.sh eb3_2_clone2Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone2Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone2

./get.tl.counts.hg19.sh eb3_2_clone3Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone3Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone3

./get.tl.counts.hg19.sh eb3_2_clone4Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone4Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone4

./get.tl.counts.hg19.sh eb3_2_clone10Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone10Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone10

./get.tl.counts.hg19.sh eb3_2_clone13Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone13Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone13 

./get.tl.counts.hg19.sh eb3_2_clone15Aligned.out.samtool.rmdup.plus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb3_2_clone15Aligned.out.samtool.rmdup.minus.allele.counts.haplotype.resolved.counts.rmv.blck.bed \
  eb32.all.clones.tls.plus.sorted.bed \
  eb32.all.clones.tls.minus.sorted.bed \
  eb32_clone15


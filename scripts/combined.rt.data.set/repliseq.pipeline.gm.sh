#!/bin/bash

# align haplotypes
python align.haplotypes.py --bed gm12878_clone4_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone4_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone5_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed gm12878_clone5_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./

# remove blacklist SNPs before LSM
./remove.blacklist.hg38.sh gm12878_clone4_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone4_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh gm12878_clone5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed

# do LSM after the bad snps are removed
python lsm.normalize.allele.counts.py gm12878_clone4_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py gm12878_clone4_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py gm12878_clone5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py gm12878_clone5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

## overlap with windows. these windows also have the blacklist subregions removed
./make.allele.window.counts.hg38.sh gm12878_clone4_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh gm12878_clone4_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh gm12878_clone5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh gm12878_clone5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed

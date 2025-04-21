#!/bin/bash

## align haplotypes
python align.haplotypes.py --bed acp6_c1_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c1_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c2_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c2_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c5_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c5_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c6_early.hg38.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed acp6_c6_late.hg38.sorted.markdup.allele.counts.bed --out_directory ./

# remove blacklist SNPs before LSM
./remove.blacklist.hg38.sh acp6_c1_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c1_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed 
./remove.blacklist.hg38.sh acp6_c2_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c2_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed 
./remove.blacklist.hg38.sh acp6_c5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed 
./remove.blacklist.hg38.sh acp6_c6_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.hg38.sh acp6_c6_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.bed

# do LSM after the bad snps are removed
python lsm.normalize.allele.counts.py acp6_c1_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py acp6_c1_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed 
python lsm.normalize.allele.counts.py acp6_c2_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py acp6_c2_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed 
python lsm.normalize.allele.counts.py acp6_c5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py acp6_c5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed 
python lsm.normalize.allele.counts.py acp6_c6_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py acp6_c6_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

## overlap with windows. these windows also have the blacklist subregions removed
./make.allele.window.counts.hg38.sh acp6_c1_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c1_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c2_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c2_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c5_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c5_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c6_early.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg38.sh acp6_c6_late.hg38.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed

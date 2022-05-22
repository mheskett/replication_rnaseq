#!/bin/bash

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha2.rna.500kb.bed
bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.2e.haplotype.resolved.counts.sorted.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.2l.haplotype.resolved.counts.sorted.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed
#####
bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha10.rna.500kb.bed
bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.10e.haplotype.resolved.counts.sorted.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.10l.haplotype.resolved.counts.sorted.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed

###
bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rna.500kb.bed
   bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/4e.combined.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/4l.combined.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.500kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rna.500kb.bed
   bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/5e.combined.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/5l.combined.samtool.rmdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.500kb.bed

#!/bin/bash
## first do
## bedtools sort -i gm12878.rep1.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed -g ../annotation_files/human_g1k_genome.fa.fai etecadfasdf

bedtools map -a gm12878.rep1.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.sorted.bed \
  -b 4l.combined.overlap.na12878.hg19.haplotype.resolved.counts.sorted.bed -o sum,sum -c 6,7 | awk '{FS="\t";}$15!="."{print $0}' | awk '{FS="\t";}$16!="."{print $0}' > 4l.combined.gm12878.rep1.windows.bed

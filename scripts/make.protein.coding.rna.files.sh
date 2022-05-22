bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha2.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha2.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha2.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha2.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha2.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.3Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha3.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.3Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha3.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha3.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha3.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha3.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.4Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha4.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.4Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha4.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha4.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha4.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha4.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha10.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha10.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha10.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha10.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha10.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.13Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha13.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.13Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha13.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha13.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha13.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha13.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.15Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha15.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.15Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha15.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha15.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/bouha15.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > bouha15.protein.coding.all.counts.bed

#######################
bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/bouha2.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.4.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > gm12878.4.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/gm12878.5.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.5.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > gm12878.5.protein.coding.all.counts.bed

bedtools map -a ucsc.genes.plus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts//gm12878.rep1.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.protein.coding.plus.rna.bed
   bedtools map -a ucsc.genes.minus.hg19.bed -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep1.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed \
   -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.protein.coding.minus.rna.bed
cat /Users/mike/replication_rnaseq/all.final.data/gm12878.5.protein.coding.plus.rna.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.protein.coding.minus.rna.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > gm12878.rep1.protein.coding.all.counts.bed

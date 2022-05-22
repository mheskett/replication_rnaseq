sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed | 
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rep1.rep2.vlincs.plus.bed
   
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.4x.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rep1.rep2.vlincs.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rep1.rep2.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rep1.rep2.vlincs.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.4.rep1.rep2.vlincs.all.bed
#####################
 sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed | 
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.plus.bed
   
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.5x.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.all.bed
################### rep 1 counts using rep1rep2 vlincs
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep1.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep1.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed
################### rep 2 counts using rep1rep2 vlincs
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep2.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep2.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.all.bed

### parent cell line. rep1 RNA with rep1 vlincs #### watch out this is overwriting the file from above........
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.rep1.hg19Aligned.out.samtool.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep1.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.rep1.hg19Aligned.out.samtool.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/gm12878.rna.allele.counts/gm12878.rep1.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed
### bouha cell lines when called with all 6 combined
#####
####
###
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.10Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.bed
##
##
##
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.13Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.13Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.bed
#########
####
####
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.15Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.15Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.bed
##
###
####
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.3Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.3Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.bed

##
###
###
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.2Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.bed
##
##3
###
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.4Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.plus.bed

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha.all.rmdup.all.chrom.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/bouhassira.rna.allele.counts/bouha.trim.4Aligned.samtool.rmdup.minus.all.chrom.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.minus.bed

cat /Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.plus.bed /Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.minus.bed |
   sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.bed

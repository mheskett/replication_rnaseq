# sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.plus.1000.10000.50000.vlinc.discovery.bed |
#        bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K9me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
#       -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9me3-human.vlincs.plus.bed

# sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.rep2.rmdup.intergenic.minus.1000.10000.50000.vlinc.discovery.bed |
#    bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K9me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
#     -o sum,sum -c 6,7 > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9me3-human.vlincs.minus.bed

# cat /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9me3-human.vlincs.plus.bed /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9me3-human.vlincs.minus.bed |
#    sort -k1,1 -k2,2n | awk '{ if ( $8 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlinc.H3K9me3-human.bed


   ####
   ####
   ####
#    ###

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K9me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9me3-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H2AFZ-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H2AFZ-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K36me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K36me3-human.bed
#######

sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K4me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K4me3-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K27ac-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K27ac-human.bed

###
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K4me1-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K4me1-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K79me2-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K79me2-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H4K20me1-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H4K20me1-human.bed
####
sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K27me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K27me3-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K4me2-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K4me2-human.bed


sort -k1,1 -k2,2n /Users/mike/replication_rnaseq/scripts/gm12878.rep1.vlincs.all.bed |
       bedtools map -a stdin -b /Users/mike/replication_rnaseq/scripts/H3K9ac-human.all.nochr.allele.counts.haplotype.resolved.counts.bed \
      -o sum,sum -c 6,7 | awk '{ if ( $10 != "." ) { print $0; } }' > /Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.H3K9ac-human.bed




# H3K9me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H2AFZ-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K36me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K4me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K27ac-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K4me1-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K79me2-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H4K20me1-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K27me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K4me2-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
# H3K9ac-human.all.nochr.allele.counts.haplotype.resolved.counts.bed
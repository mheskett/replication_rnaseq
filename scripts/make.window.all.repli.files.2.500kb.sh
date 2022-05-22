bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.15e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.15l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.13e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.13l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.3e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.3l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w500kb.s500kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.4e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.4l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed


### now use 250kb windows
bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w250kb.s250kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.15e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.15l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.250kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w250kb.s250kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.13e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.13l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.250kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w250kb.s250kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.3e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.3l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.250kb.bed

bedtools map -a /Users/mike/replication_rnaseq/all.final.data/window.file/human_g1k_v37.w250kb.s250kb.bed -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.4e.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   bedtools map -a stdin -b /Users/mike/replication_rnaseq/all.final.data/repliseq.haplotype.resolved/bouha.4l.sorted.markdup.allele.counts.haplotype.resolved.counts.bed -o sum,sum -c 6,7 \
        > /Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.250kb.bed

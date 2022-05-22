bedtools map -a hg38.nochr.25kb.windows.bed -b gm12878.rep1Aligned.hcgatk3.overlap.platinum.haplotype.resolved.bed -c 6,7 -o sum | grep -Fv "." | less

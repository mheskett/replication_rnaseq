

bcftools mpileup -R ACP6_final_bothparents_phased.threecol.bed \
  -f hg19.ethan.fa \
  -a DP,AD \
  -q 20 RNA210216MT_3_tACP6c5_hg19Aligned.out.sorted.bam > RNA2102test.vcf

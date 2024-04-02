

## files already appear to be sorted (samtools style) and RMDUPed
samtools.process.bam.sh /Volumes/Untitled/ACP6_Repli/acp6_bam_merged/acp6_c1_m_e.bam /Volumes/Untitled/ACP6_Repli/acp6_bam_merged/


python plot.acp.coding.genes.py RNA210216MT_1_tACP6c1.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_1_tACP6c1.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_1_tACP7c1.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_1_tACP7c1.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_2_tACP6c2.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_2_tACP6c2.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_3_tACP6c5.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_3_tACP6c5.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_4_tACP6c6.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_4_tACP6c6.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_6_tACP7c2.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_6_tACP7c2.coding.aei.txt
python plot.acp.coding.genes.py RNA210216MT_7_tACP7c4.paste.het.snps.haplotype.resolved.counts.genes.bed RNA210216MT_7_tACP7c4.coding.aei.txt





bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_1_tACP6c1.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "." > RNA210216MT_1_tACP6c1.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_1_tACP7c1.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_1_tACP7c1.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_2_tACP6c2.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_2_tACP6c2.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_3_tACP6c5.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_3_tACP6c5.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_4_tACP6c6.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_4_tACP6c6.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_6_tACP7c2.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_6_tACP7c2.paste.het.snps.haplotype.resolved.counts.genes.bed
bedtools map -a ucsc.coding.genes.only.including.isoforms.chr.filtered.final.bed -b RNA210216MT_7_tACP7c4.paste.het.snps.haplotype.resolved.counts.bed -o sum,sum -c 6,7 | \
   grep -Fv "."> RNA210216MT_7_tACP7c4.paste.het.snps.haplotype.resolved.counts.genes.bed


python plot.repliseq.final.acps.py --maternal_rt acp6_c1_m.bedgraph --paternal_rt acp6_c1_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp6_c2_m.bedgraph --paternal_rt acp6_c2_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp6_c5_m.bedgraph --paternal_rt acp6_c5_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp6_c6_m.bedgraph --paternal_rt acp6_c6_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp7_c2_m.bedgraph --paternal_rt acp7_c2_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp7_c4_m.bedgraph --paternal_rt acp7_c4_p.bedgraph
python plot.repliseq.final.acps.py --maternal_rt acp7_c5_m.bedgraph --paternal_rt acp7_c5_p.bedgraph



## make proper bed file from table file
tail -n +2 ACP6_final_bothparents_phased.table | awk 'OFS="\t"{split($5,a,"|");print $1,$2-1,$2,a[1],a[2]}' > ACP6_final_bothparents_phased.haplotypes.bed
tail -n +2 ACP7_final_bothparents_phased.table | awk 'OFS="\t"{split($5,a,"|");print $1,$2-1,$2,a[1],a[2]}' > ACP7_final_bothparents_phased.haplotypes.bed
###
./variants.to.table.get.counts.sh RNA210216MT_6_tACP7c2_rm.allele.counts.vcf
./variants.to.table.get.counts.sh RNA210216MT_7_tACP7c4_rm.allele.counts.vcf
./variants.to.table.get.counts.sh RNA210216MT_1_tACP7c1_rm.allele.counts.vcf


##
 tail -n +2 RNA210216MT_1_tACP6c1.allele.counts.table  | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_1_tACP6c1.allele.counts.bed
 tail -n +2 RNA210216MT_2_tACP6c2.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_2_tACP6c2.allele.counts.bed
 tail -n +2 RNA210216MT_3_tACP6c5.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_3_tACP6c5.allele.counts.bed
 tail -n +2 RNA210216MT_4_tACP6c6.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_4_tACP6c6.allele.counts.bed
###
tail -n +2 RNA210216MT_6_tACP7c2_rm.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_6_tACP7c2_rm.allele.counts.bed
tail -n +2 RNA210216MT_7_tACP7c4_rm.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_7_tACP7c4_rm.allele.counts.bed
tail -n +2 RNA210216MT_1_tACP7c1_rm.allele.counts.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > RNA210216MT_1_tACP7c1_rm.allele.counts.bed



 bedtools intersect -a RNA210216MT_1_tACP6c1.allele.counts.bed -b ACP6_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_1_tACP6c1.paste.het.snps.bed
 bedtools intersect -a RNA210216MT_2_tACP6c2.allele.counts.bed -b ACP6_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_2_tACP6c2.paste.het.snps.bed
 bedtools intersect -a RNA210216MT_3_tACP6c5.allele.counts.bed -b ACP6_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_3_tACP6c5.paste.het.snps.bed
 bedtools intersect -a RNA210216MT_4_tACP6c6.allele.counts.bed -b ACP6_final_bothparents_phased.haplotypes.bed -wa -wb >  RNA210216MT_4_tACP6c6.paste.het.snps.bed

###
bedtools intersect -a RNA210216MT_6_tACP7c2_rm.allele.counts.bed -b ACP7_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_6_tACP7c2.paste.het.snps.bed
bedtools intersect -a RNA210216MT_7_tACP7c4_rm.allele.counts.bed -b ACP7_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_7_tACP7c4.paste.het.snps.bed
bedtools intersect -a RNA210216MT_1_tACP7c1_rm.allele.counts.bed -b ACP7_final_bothparents_phased.haplotypes.bed -wa -wb > RNA210216MT_1_tACP7c1.paste.het.snps.bed




python align.haplotypes.py --bed RNA210216MT_1_tACP6c1.paste.het.snps.bed --out_directory ./
python align.haplotypes.py --bed RNA210216MT_2_tACP6c2.paste.het.snps.bed --out_directory ./
python align.haplotypes.py --bed RNA210216MT_3_tACP6c5.paste.het.snps.bed --out_directory ./
python align.haplotypes.py --bed RNA210216MT_4_tACP6c6.paste.het.snps.bed --out_directory ./

python align.haplotypes.py --bed RNA210216MT_6_tACP7c2.paste.het.snps.bed --out_directory ./
python align.haplotypes.py --bed RNA210216MT_7_tACP7c4.paste.het.snps.bed --out_directory ./
python align.haplotypes.py --bed RNA210216MT_1_tACP7c1.paste.het.snps.bed --out_directory ./




RNA210216MT_1_tACP6c1.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_1_tACP7c1.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_2_tACP6c2.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_3_tACP6c5.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_4_tACP6c6.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_6_tACP7c2.paste.het.snps.haplotype.resolved.counts.bed
RNA210216MT_7_tACP7c4.paste.het.snps.haplotype.resolved.counts.bed

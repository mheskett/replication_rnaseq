awk '$7!=$8{print $0}' mm9.refseq.txt | awk 'OFS="\t"{print $3,$5,$6,$13,0,$4}'  | tail -n +2 | sort -u -k4,4 | less | sort -k1,1 -k2,2n > mm9.refseq.cds.only.first.isoform.bed

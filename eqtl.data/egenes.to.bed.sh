#!/bin/bash
file=$(basename "$1" .txt)

head -n 1 $1 |  awk 'OFS="\t"{print $3,$4,$5,$2,$7,$6,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' > $file.header.txt

echo $file
tail -n +2 $1 | awk 'OFS="\t"{print $3,$4,$5,$2,$7,$6,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' | awk 'BEGIN{OFS="\t";}$28<=0.05{print $0}' | \
  awk 'BEGIN{OFS="\t";}NF==30{print $0}' | bedtools sort -i stdin -g /Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai > $file.significant.bed

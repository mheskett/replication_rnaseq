#!/bin/bash
bedtools multiinter -g hg38.windows.w250.s250.bed \
  -i acp6.rt.vert.sorted.4col.merged.bed \
     acp7.rt.vert.sorted.4col.merged.bed \
     gm.rt.vert.sorted.4col.merged.bed \
     eb.rt.vert.hg38.lifted.sorted.4col.merged.bed \
     mouse.rt.vert.sorted.4col.merged.bed \
  -header -names acp6 acp7 gm12878 eb3_2 c57/b6 | awk '$5!="c57/b6"{print $0}' | grep -v chrX > rt.vert.all.samples.txt

tail -n +2 rt.vert.all.samples.txt | bedtools merge -d 250000 -i stdin  > vert.coordinates.merged.txt

## make colons for gprofiler
awk '{print $1":"$2":"$3}' vert.coordinates.merged.txt > vert.coordinates.merged.colons.txt

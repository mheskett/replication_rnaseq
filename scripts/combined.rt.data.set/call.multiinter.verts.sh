#!/bin/bash
bedtools multiinter -g hg38.windows.w250.s250.bed \
  -i acp6.rt.vert.sorted.4col.merged.bed \
     acp7.rt.vert.sorted.4col.merged.bed \
     gm.rt.vert.sorted.4col.merged.bed \
     eb.rt.vert.hg38.lifted.sorted.4col.merged.bed \
  -header -names acp6 acp7 gm12878 eb3_2 | grep -v chrX > rt.vert.all.samples.txt
tail -n +2 rt.vert.all.samples.txt | bedtools merge -d 250000 -i stdin  > vert.coordinates.merged.all.samples.txt
## make colons for gprofiler
awk '{print $1":"$2":"$3}' vert.coordinates.merged.all.samples.txt > vert.coordinates.merged.all.samples.colons.txt


#### just ACPs
bedtools multiinter -g hg38.windows.w250.s250.bed \
  -i acp6.rt.vert.sorted.4col.merged.bed \
     acp7.rt.vert.sorted.4col.merged.bed \
  -header -names acp6 acp7 | grep -v chrX > rt.vert.acps.txt
tail -n +2 rt.vert.acps.txt | bedtools merge -d 250000 -i stdin  > vert.coordinates.merged.acps.txt
awk '{print $1":"$2":"$3}' vert.coordinates.merged.acps.txt > vert.coordinates.merged.acps.colons.txt


#### just LCLs
bedtools multiinter -g hg38.windows.w250.s250.bed \
  -i gm.rt.vert.sorted.4col.merged.bed \
      eb.rt.vert.hg38.lifted.sorted.4col.merged.bed \
  -header -names gm12878 eb3_2 | grep -v chrX > rt.vert.lcls.txt
tail -n +2 rt.vert.lcls.txt | bedtools merge -d 250000 -i stdin  > vert.coordinates.merged.lcls.txt
awk '{print $1":"$2":"$3}' vert.coordinates.merged.lcls.txt > vert.coordinates.merged.lcls.colons.txt










# not including mouse
##     mouse.rt.vert.sorted.4col.merged.bed \

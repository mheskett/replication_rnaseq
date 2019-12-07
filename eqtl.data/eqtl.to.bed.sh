#!/bin/bash
file=$(basename "$1" .txt)

tail -n +2 $1 | awk 'OFS="\t"{split($1,a,"_")}{print a[1],a[2]-1,a[2],a[3],a[4],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $file.bed


#!/bin/bash

# separates haplotype resolved file into paternal and maternal alleles

awk 'OFS="\t"{split($5,a,"|");print $1,$2-1,$2,a[1],a[2]}' NA12878.nochr.table | tail -n +2 > NA12878.nochr.bed

awk 'BEGIN{OFS="\t";}$4!=$5{print $1,$2,$3,$4,$5}' NA12878.nochr.bed > NA12878.nochr.het.bed

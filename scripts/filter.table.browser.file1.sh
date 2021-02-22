#!/bin/bash

sed 's/chr//g' $1 | grep -v "^GL" | grep -v "^M" | grep -v "^Y" | grep -v "^fix" | grep -v "^alt" | grep -v "^hap" | grep -v "^random" | grep -v "^17_ctg" | grep -v "^4_ctg" | \
  grep -v "^9_gl" | grep -v "^7_gl" | grep -v "^6_qbl_hap6" | grep -v "^6_mcf_hap5" | grep -v "^6_mann_hap4" | grep -v "^6_dbb_hap3" | grep -v "^6_cox_hap2" | grep -v "^6_apd_hap1" | grep -v "^4_gl" | grep -v "^1_gl" | grep -v "^17_gl" | grep -v "^Un_" | grep -v "^6_ssto" | grep -v "^19_gl" | bedtools sort -i stdin -g human_g1k_v37.fasta.fai > $1.filtered.bed

#!/bin/bash

sed 's/chr//g' $1 | grep -v M | grep -v Y | grep -v fix | grep -v alt | grep -v hap | grep -v random | grep -v 17_ctg | grep -v 4_ctg \
  | grep -v Un_ | grep -v 6_ssto | grep -v 19_gl | bedtools sort -i stdin -g  /Users/heskett/replication_rnaseq/annotation.files/genome.fa.fai > $2

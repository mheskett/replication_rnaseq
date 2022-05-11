#!/bin/bash

input=$1
filename=${input%.*}

srun -p light --mem=30000 --time=200 gatk VariantsToTable \
  -O $filename.table \
  -F CHROM -F POS -F REF -F ALT \
  -GF GT \
  -V $input

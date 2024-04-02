#!/bin/bash

input=$1
filename=${input%.*}

/Users/heskett/replication_rnaseq/scripts/gatk-4.2.1.0/gatk VariantsToTable -V $input -O $filename.table -F CHROM -F POS -F REF -F ALT -GF AD

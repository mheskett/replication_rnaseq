#!/bin/bash

srun --mem 8000 vcftools --remove-indels --gzvcf na12878.nochr.hg19.vcf.gz --recode-INFO-all --recode --out ./no.indels

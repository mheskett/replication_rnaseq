#!/bin/bash

sbatch slurm.vlinc.identifier.hg19.sh gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ./ 5000 50000 7500
sbatch slurm.vlinc.identifier.hg19.sh gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ./ 5000 50000 10000
sbatch slurm.vlinc.identifier.hg19.sh gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ./ 8000 50000 7500
sbatch slurm.vlinc.identifier.hg19.sh gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ./ 8000 50000 10000

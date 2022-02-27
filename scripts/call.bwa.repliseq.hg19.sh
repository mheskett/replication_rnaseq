#!/bin/bash

sbatch slurm.call.bwa.hg19.sh ../repli-seq/gilbert.data.nov.19/700/4e.R1.combined.fastq.gz ../repli-seq/gilbert.data.nov.19/700/4e.R2.combined.fastq.gz 4e.combined ./
sbatch slurm.call.bwa.hg19.sh ../repli-seq/gilbert.data.nov.19/700/4l.R1.combined.fastq.gz ../repli-seq/gilbert.data.nov.19/700/4l.R2.combined.fastq.gz 4l.combined ./
sbatch slurm.call.bwa.hg19.sh ../repli-seq/gilbert.data.nov.19/700/5e.R1.combined.fastq.gz ../repli-seq/gilbert.data.nov.19/700/5e.R2.combined.fastq.gz 5e.combined ./
sbatch slurm.call.bwa.hg19.sh ../repli-seq/gilbert.data.nov.19/700/5l.R1.combined.fastq.gz ../repli-seq/gilbert.data.nov.19/700/5l.R2.combined.fastq.gz 5l.combined ./

## 4DN stuff is single ended.
sbatch slurm.call.bwa.single.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/4dn.gm12878/gm12878.4dn.rep1.early.fastq ./gm12878.4dn.rep1.early
sbatch slurm.call.bwa.single.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/4dn.gm12878/gm12878.4dn.rep2.early.fastq ./gm12878.4dn.rep2.early
sbatch slurm.call.bwa.single.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/4dn.gm12878/gm12878.4dn.rep1.late.fastq ./gm12878.4dn.rep1.late
sbatch slurm.call.bwa.single.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/4dn.gm12878/gm12878.4dn.rep2.late.fastq ./gm12878.4dn.rep2.late

## stuff from dec 2020

sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.dec.20/bouha.10.e.r1.fastq.gz ../bouhassira.repliseq.dec.20/bouha.10.e.r2.fastq.gz bouha.10e ../bouhassira.repliseq.dec.20/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.dec.20/bouha.10.l.r1.fastq.gz ../bouhassira.repliseq.dec.20/bouha.10.l.r2.fastq.gz bouha.10l ../bouhassira.repliseq.dec.20/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.dec.20/bouha.2.e.r1.fastq.gz ../bouhassira.repliseq.dec.20/bouha.2.e.r2.fastq.gz bouha.2e ../bouhassira.repliseq.dec.20/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.dec.20/bouha.2.l.r1.fastq.gz ../bouhassira.repliseq.dec.20/bouha.2.l.r2.fastq.gz bouha.2l ../bouhassira.repliseq.dec.20/


## call stuff from november 2021

sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt3E_S1_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt3E_S1_R2_001.fastq.gz bouha.3e ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt3L_S2_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt3L_S2_R2_001.fastq.gz bouha.3l ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt4E_S3_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt4E_S3_R2_001.fastq.gz bouha.4e ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt4L_S4_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt4L_S4_R2_001.fastq.gz bouha.4l ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt13E_S5_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt13E_S5_R2_001.fastq.gz bouha.13e ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt13L_S6_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt13L_S6_R2_001.fastq.gz bouha.13l ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt15E_S7_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt15E_S7_R2_001.fastq.gz bouha.15e ../bouhassira.repliseq.nov.21/
sbatch slurm.call.bwa.hg19.sh ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt15L_S8_R1_001.fastq.gz ../bouhassira.repliseq.nov.21/211115_A01058_0187_BHN5WKDRXY/LIB211029MT/LIB211029MT_mt15L_S8_R2_001.fastq.gz bouha.15l ../bouhassira.repliseq.nov.21/

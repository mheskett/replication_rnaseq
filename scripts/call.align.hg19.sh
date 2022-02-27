
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/p175.x.R1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/p175.x.R2.fq.gz p175.x.hg19

# this is DNA from repli-seq...should be aligned with bwa
#sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/fastq/4el700.R1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/fastq/4el700.R2.fastq.gz 4el700
#sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/fastq/5el700.R1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/repli-seq/fastq/5el700.R2.fastq.gz 5el700

sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM12878.4x_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM12878.4x_2.fq.gz gm12878.4x.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM12878.5x_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM12878.5x_2.fq.gz gm12878.5x.hg19

# these were made by soren impey and are crap
#sbatch slurm.align.paired.hg19.sh ../fastq/a.27.s34.R1.fastq.gz ../fastq/b.27.s34.R2.fastq.gz neb27
#sbatch slurm.align.paired.hg19.sh ../fastq/e.23.s32.R1.fastq.gz ../fastq/f.23.s32.R2.fastq.gz neb23
#sbatch slurm.align.paired.hg19.sh ../fastq/g.22.s31.R1.fastq.gz ../fastq/h.22.s31.R2.fastq.gz neb22
#sbatch slurm.align.paired.hg19.sh ../fastq/i.12.s30.R1.fastq.gz ../fastq/j.12.s30.R2.fastq.gz neb12
#sbatch slurm.align.paired.hg19.sh ../fastq/c.25.s33.R1.fastq.gz ../fastq/d.25.s33.R2.fastq.gz neb25

### Cytosolic RNA-seq to compare with nuclear for GM12878
sbatch slurm.align.paired.hg19.sh ../fastq/wgEncodeCshlLongRnaSeqGm12878CytosolLongnonpolyaFastqRd1Rep1.fastq.gz ../fastq/wgEncodeCshlLongRnaSeqGm12878CytosolLongnonpolyaFastqRd2Rep1.fastq.gz gm12878.cytosol.rep1.hg19
sbatch slurm.align.paired.hg19.sh ../fastq/wgEncodeCshlLongRnaSeqGm12878CytosolLongnonpolyaFastqRd1Rep2.fastq.gz ../fastq/wgEncodeCshlLongRnaSeqGm12878CytosolLongnonpolyaFastqRd2Rep2.fastq.gz gm12878.cytosol.rep2.hg19

sbatch slurm.align.paired.hg19.sh ../fastq/wgEncodeCshlLongRnaSeqGm12878NucleusLongnonpolyaFastqRd1Rep1.fastq.gz ../fastq/wgEncodeCshlLongRnaSeqGm12878NucleusLongnonpolyaFastqRd2Rep1.fastq.gz gm12878.rep1.hg19
sbatch slurm.align.paired.hg19.sh ../fastq/wgEncodeCshlLongRnaSeqGm12878NucleusLongnonpolyaFastqRd1Rep2.fastq.gz ../fastq/wgEncodeCshlLongRnaSeqGm12878NucleusLongnonpolyaFastqRd2Rep2.fastq.gz gm12878.rep2.hg19

### hepg2 and k562 RNA-seq cytosol vs nuclear
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562CytosolLongnonpolyaFastqRd1Rep1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562CytosolLongnonpolyaFastqRd2Rep1.fastq.gz k562.cytosol.rep1.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562CytosolLongnonpolyaFastqRd1Rep2.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562CytosolLongnonpolyaFastqRd2Rep2.fastq.gz k562.cytosol.rep2.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562NucleusLongnonpolyaFastqRd1Rep1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562NucleusLongnonpolyaFastqRd2Rep1.fastq.gz k562.nucleus.rep1.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562NucleusLongnonpolyaFastqRd1Rep2.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqK562NucleusLongnonpolyaFastqRd2Rep2.fastq.gz k562.nucleus.rep2.hg19

sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2NucleusLongnonpolyaFastqRd1Rep1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2NucleusLongnonpolyaFastqRd2Rep1.fastq.gz hepg2.nucleus.rep1.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2NucleusLongnonpolyaFastqRd1Rep2.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2NucleusLongnonpolyaFastqRd2Rep2.fastq.gz hepg2.nucleus.rep2.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2CytosolLongnonpolyaFastqRd1Rep1.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2CytosolLongnonpolyaFastqRd2Rep1.fastq.gz hepg2.cytosol.rep1.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2CytosolLongnonpolyaFastqRd1Rep2.fastq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/wgEncodeCshlLongRnaSeqHepg2CytosolLongnonpolyaFastqRd2Rep2.fastq.gz hepg2.cytosol.rep2.hg19


### gm12878 clones
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_41_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_41_2.fq.gz gm12878.41.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_42_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_42_2.fq.gz gm12878.42.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_43_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_43_2.fq.gz gm12878.43.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_51_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_51_2.fq.gz gm12878.51.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_52_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_52_2.fq.gz gm12878.52.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_53_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/GM878_53_2.fq.gz gm12878.53.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_1_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_1_2.fq.gz p175.1.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_2_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_2_2.fq.gz p175.2.hg19
sbatch slurm.align.paired.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_3_1.fq.gz /home/groups/Spellmandata/heskett/replication.rnaseq/fastq/novogene.rnaseq/H202SC19020025/raw_data/P175_3_2.fq.gz p175.3.hg19

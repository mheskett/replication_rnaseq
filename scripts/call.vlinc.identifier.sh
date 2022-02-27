sbatch slurm.vlinc.identifier.sh p175.xAligned.out.samtool.rmdup.bam ./
sbatch slurm.vlinc.identifier.sh gm12878.5xAligned.out.samtool.rmdup.bam ./
sbatch slurm.vlinc.identifier.sh gm12878.4xAligned.out.samtool.rmdup.bam ./
sbatch slurm.vlinc.identifier.sh gm12878.rep1Aligned.out.samtool.rmdup.bam ./
sbatch slurm.vlinc.identifier.sh gm12878.rep2Aligned.out.samtool.rmdup.bam ./


sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/p175.x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.5x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.4x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/hepg2.nucleus.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/hepg2.nucleus.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/k562.nucleus.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/k562.nucleus.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000

## intergenic
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/p175.x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.5x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.4x.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/hepg2.nucleus.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/hepg2.nucleus.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/k562.nucleus.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/encode.rnaseq.bams/k562.nucleus.rep2.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
## bouhassira cell line clones
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.10Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.13Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.15Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.2Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.3Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.trim.4Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
### all gm12887 and bouhassira compbined. increase the memory for these
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/gm12878.rep1.rep2.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.chr22.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/gm12878.rep1.hg19Aligned.out.samtool.rmdup.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.10.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.11.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.12.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.13.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.14.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.15.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.16.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.17.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.18.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.19.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.1.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.20.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.21.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.22.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.2.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.3.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.4.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.5.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.6.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.7.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.8.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.9.bam ../results/ 1000 10000 50000
sbatch slurm.vlinc.identifier.intergenic.only.hg19.sh /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/bouha.all.rmdup.bam.X.bam ../results/ 1000 10000 50000

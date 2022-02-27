srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.4x.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.5x.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.cytosol.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.cytosol.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh gm12878.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh hepg2.cytosol.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh hepg2.cytosol.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh hepg2.nucleus.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh hepg2.nucleus.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh k562.cytosol.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh k562.cytosol.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh k562.nucleus.rep1.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh k562.nucleus.rep2.hg19Aligned.out.bam ./ &
srun --mem 30000 -c 8 --time 1000 samtools.process.bam.sh p175.x.hg19Aligned.out.bam ./ &

## bouhassira cell lines

srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.10Aligned.out.bam ./ &
srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.13Aligned.out.bam ./ &
srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.15Aligned.out.bam ./ &
srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.2Aligned.out.bam ./ &
srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.3Aligned.out.bam ./ &
srun --mem 50000 -c 8 --time 1000 samtools.process.bam.sh bouha.trim.4Aligned.out.bam ./ &

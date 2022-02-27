srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.2Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.3Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.4Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.10Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.13Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh bouha.trim.15Aligned.out.samtool.rmdup.bam ./ &

srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.2Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.3Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.4Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.10Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.13Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh bouha.trim.15Aligned.out.samtool.rmdup.bam ./ &

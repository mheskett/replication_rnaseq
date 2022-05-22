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

## acp
srun --mem 16000 --time 400 get.plus.reads.sh acp1Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh acp3Aligned.out.samtool.rmdup.bam ./ &  
srun --mem 16000 --time 400 get.plus.reads.sh acp5Aligned.out.samtool.rmdup.bam ./ &  
srun --mem 16000 --time 400 get.plus.reads.sh acp7Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.plus.reads.sh acp2Aligned.out.samtool.rmdup.bam ./ &  
srun --mem 16000 --time 400 get.plus.reads.sh acp4Aligned.out.samtool.rmdup.bam ./ &  
srun --mem 16000 --time 400 get.plus.reads.sh acp6Aligned.out.samtool.rmdup.bam ./ &  
srun --mem 16000 --time 400 get.plus.reads.sh acp8Aligned.out.samtool.rmdup.bam ./ &

srun --mem 16000 --time 400 get.minus.reads.sh acp1Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp3Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp5Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp7Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp2Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp4Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp6Aligned.out.samtool.rmdup.bam ./ &
srun --mem 16000 --time 400 get.minus.reads.sh acp8Aligned.out.samtool.rmdup.bam ./ &

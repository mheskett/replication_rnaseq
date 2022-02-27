srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ nd1.minus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ nd1.plus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ nd2.minus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ nd2.plus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ td1.plus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ td2.minus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ td2.plus.r1.fastq.gz &
srun -c 6 --mem 6000 --time 90 trim_galore --cores 6 -o ./ td.minus.r1.fastq.gz &

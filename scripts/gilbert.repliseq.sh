## see nature protocol Genome-wide analysis of replication timing by next-generation sequencing with E/L Repli-seq 2018.


## Map the .fastq files and generate the coverage files by running the following commands:
for file in *R1.fastq*; do
	bowtie2 -x path_to_your_genome --no-mixed --no-discordant --reorder -X 1000 -1 $file -2 ${file%R1.fastq*}R2.fastq* -S ${file%R1.fastq*}.sam 2>> ${file%R1.fastq*}mapping_log.txt
	samtools view -bSq 20 ${file%R1.fastq*}.sam > ${file%R1.fastq*}.bam
	samtools sort -o ${file%R1.fastq*}_srt.bam ${file%R1.fastq*}.bam
	samtools rmdup -S ${file%R1.fastq*}_srt.bam ${file%R1.fastq*}_rmdup.bam
	bedtools bamtobed -i ${file%R1.fastq*}_rmdup.bam | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > ${file%R1.fastq*}.bed
	x=`wc -l ${file%R1.fastq*}.bed | cut -d' ' -f 1`
	bedtools intersect -sorted -c -b ${file%R1.fastq*}.bed -a your_genome_windows.bed | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > ${file%R1.fastq*}.bg
done

#Calculate RT by running the following commands:
for file in *_E_.bg; do
  paste $file ${file%E_.bg}L_.bg | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%E_.bg}T_.bg
done
echo -e "chr\tstart\tstop\t"`ls *T_.bg` | sed 's/\ /\t/g' > merge_RT.txt

## Merge the RT files. If you have only one sample, run the following command:
cat *T_.bg >> merge_RT.txt

## If you have multiple samples, use the following command:
bedtools unionbedg -filler "NA" -i *T_.bg >> merge_RT.txt

### See R script for next steps involving normalization
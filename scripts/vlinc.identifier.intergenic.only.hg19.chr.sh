#!/bin/bash

## This program utilizes samtools and bedtools to create a call set of non-coding 
## "transcribed loci" that are greater than 50kb in length. We consider
## these as putative very long non-coding RNAs. Important considerations are 
## that memory requirements can get very large when commands are chained,
## thus writing intermediate files to storage is necessary. I found that
## the strand-specific feature on bedtools may not work and will not give
## an error message, so separating files based on strand is easier for now. 
## Arguments include the input BAM files, output directory, and the 
## size parameters for merging reads and filtering. 
## For joint calling of multiple samples, one may combine BAM files before
## running this script, or call individually and then use a secondary
## strategy to filter appropriately and make consensus calls on the 
## start and end points of each result.


## bam file
input=$1
out_dir=$2
## removes /path/to/file
b=$(basename $input) 
### removes file extension
filename=${b%.*} 
### input parameters on the command line
first_merge=$3
second_merge=$4
filter_size=$5
## switching to samtools instead of GATK for remove duplicates. 
library_size=$(samtools view -c -F 260 $input)
echo $library_size
## get plus strand
## this is for a Paired end, strand specific library
samtools view -b -f 131 -F 16 $input > $out_dir$filename.fwd1.bam
samtools index $out_dir$filename.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.fwd2.bam

# now combine, getting all plus strand reads
samtools merge -f $out_dir$filename.plus.bam $out_dir$filename.fwd1.bam $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.plus.bam

## get minus strand
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam

## now do coverage separately. It is easier to work with separate files for each strand
bedtools genomecov -bg -split -ibam $out_dir$filename.plus.bam > $out_dir$filename.plus.cov.bed
bedtools genomecov -bg -split -ibam $out_dir$filename.minus.bam > $out_dir$filename.minus.cov.bed

## remove all overlap with known exons and introns for coding genes (whole gene) in a strand specific manner. 
## remove the poor mapping "blacklisted" regions of the genome

bedtools subtract -a $out_dir$filename.plus.cov.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/encode.blacklist.final.chr.bed \
  | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.plus.hg19.bed \
  > $out_dir$filename.plus.subtract.whole.gene.bed

bedtools subtract -a $out_dir$filename.minus.cov.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/encode.blacklist.final.chr.bed \
  | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.minus.hg19.bed \
  > $out_dir$filename.minus.subtract.whole.gene.bed

## merge segments that are separated by 500bp (caron et al.) or 1000bp (kapranov et al.) or less while subtracting all the cds again incase you merged over some tiny coding exons
bedtools merge -i $out_dir$filename.plus.subtract.whole.gene.bed -d $first_merge | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.plus.hg19.bed \
  > $out_dir$filename.plus.subtract.merge.$first_merge.bed
bedtools merge -i $out_dir$filename.minus.subtract.whole.gene.bed -d $first_merge | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.minus.hg19.bed \
  > $out_dir$filename.minus.subtract.merge.$first_merge.bed

# now merge segments if theyre less than 10kb separated (Caron et al.)
bedtools merge -i $out_dir$filename.plus.subtract.merge.$first_merge.bed -d $second_merge | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.plus.hg19.bed \
  > $out_dir$filename.plus.$first_merge.$second_merge.bed
bedtools merge -i $out_dir$filename.minus.subtract.merge.$first_merge.bed -d $second_merge | bedtools subtract -a stdin -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.minus.hg19.bed \
   > $out_dir$filename.minus.$first_merge.$second_merge.bed

## Need to subtract again after merging to get rid of regions that were joined 
## across tiny little exons.
## Only keep fragments of size 50 kb or larger as defined by Kapranov et al.
## bedtools coverage -s doesn't work as expected so im avoiding it.
## filters by min size 50kb, then removes anything that has overlap with any 
## coding exon or intron (whole gene annotation)
## Get the read coverage from the bam.
## Overlap LINE elements that have 90% intersection with the lncRNA
## Calculate RPKM and print to file

## this step takes alot of memory. For large files may need to split into more
## temporary files on storage.
echo "before big mem"
awk -v var="$filter_size" -v var2="$filename" '{OFS="\t"} $3-$2>=var{print $1, $2, $3, i++"_plus_"var2}' $out_dir$filename.plus.$first_merge.$second_merge.bed > $out_dir$filename.tmp1.plus.bed
	bedtools intersect -f 0.25 -v -a $out_dir$filename.tmp1.plus.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.hg19.bed > $out_dir$filename.tmp2.plus.bed
	bedtools coverage -c -F 0.5 -a $out_dir$filename.tmp2.plus.bed -b $out_dir$filename.plus.bam > $out_dir$filename.tmp3.plus.bed
	bedtools coverage -a $out_dir$filename.tmp3.plus.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.L1s.hg19.chr.bed -F 0.9 |
	awk -v var3="$library_size" 'OFS="\t"{print $1, $2, $3, $4, $5 / ( ($3-$2)/1000 ), "+", $9}' > $out_dir$filename.intergenic.plus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed

awk -v var="$filter_size" -v var2="$filename" -v var3="$library_size" '{OFS="\t"} $3-$2>=var{print $1, $2, $3, i++"_minus_"var2}' $out_dir$filename.minus.$first_merge.$second_merge.bed > $out_dir$filename.tmp1.minus.bed
	bedtools intersect -f 0.25 -v -a $out_dir$filename.tmp1.minus.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.known.genes.cds.chr.hg19.bed > $out_dir$filename.tmp2.minus.bed
	bedtools coverage -c -F 0.5 -a $out_dir$filename.tmp2.minus.bed -b $out_dir$filename.minus.bam > $out_dir$filename.tmp3.minus.bed
	bedtools coverage -a $out_dir$filename.tmp3.minus.bed -b /home/exacloud/gscratch/ThayerLab/heskett/acp_analysis/annotation.files/ucsc.L1s.hg19.chr.bed -F 0.9 |
	awk -v var3="$library_size" 'OFS="\t"{print $1, $2, $3, $4, $5 / ( ($3-$2)/1000 ), "-", $9}' > $out_dir$filename.intergenic.minus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed

## This will create a file specifically for the UCSC genome browser. 
awk 'OFS="\t"{print "chr"$1, $2, $3, $4, 0, $6}' $out_dir$filename.intergenic.plus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI > $out_dir$filename.intergenic.$first_merge.$second_merge.$filter_size.vlinc.discovery.plus.browser.bed
awk 'OFS="\t"{print "chr"$1, $2, $3, $4, 0, $6}' $out_dir$filename.intergenic.minus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI > $out_dir$filename.intergenic.$first_merge.$second_merge.$filter_size.vlinc.discovery.minus.browser.bed

####
####
####
####
## remove the intermediate files if they are not needed for additional visualization
# rm $out_dir$filename.plus.cov.bed
# rm $out_dir$filename.minus.cov.bed

rm $out_dir$filename.tmp1.plus.bed
rm $out_dir$filename.tmp2.plus.bed
rm $out_dir$filename.tmp3.plus.bed
rm $out_dir$filename.tmp1.minus.bed
rm $out_dir$filename.tmp2.minus.bed
rm $out_dir$filename.tmp3.minus .bed

rm $out_dir$filename.plus.subtract.whole.gene.bed
rm $out_dir$filename.minus.subtract.whole.gene.bed

rm $out_dir$filename.plus.subtract.merge.$first_merge.bed
rm $out_dir$filename.minus.subtract.merge.$first_merge.bed

rm $out_dir$filename.plus.$first_merge.$second_merge.bed
rm $out_dir$filename.minus.$first_merge.$second_merge.bed

rm $out_dir$filename.fwd1.bam*
rm $out_dir$filename.fwd2.bam*
rm $out_dir$filename.rev1.bam*
rm $out_dir$filename.rev2.bam*

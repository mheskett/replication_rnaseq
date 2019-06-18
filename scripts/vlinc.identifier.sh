#!/bin/bash

input=$1 ## bam file
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
## currently using file that's already had duplicates removed, print reads, and splitNcigar reads 
## basically the GATK prepwork already done
## switching to samtools instead of GATK for remove duplicates. print reads not necessary unless GATK being used downstream

echo $out_dir$filename
## include header and filter out map quality below 30
## this should be done by preprocessing script
#samtools view -h -b -q 30 $input > $out_dir$filename.mq30.bam

## get the per base genome coverage, while NOT including inferred coverage of spliced out introns 
## this tool is similar to bedtools depth
## genomecov strandedness doesn't seem to work. will investigate this further but for now switch to old method
#bedtools genomecov -bg -split -strand + -ibam $out_dir$filename.mq40.bam  > $out_dir$filename.mq40.plus.cov.bed
#bedtools genomecov -bg -split -strand - -ibam $out_dir$filename.mq40.bam  > $out_dir$filename.mq40.minus.cov.bed

## get plus strand
samtools view -b -f 131 -F 16 $input > $out_dir$filename.fwd1.bam
samtools index $out_dir$filename.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.fwd2.bam

# now combine, this should contail all plus strand genes
samtools merge -f $out_dir$filename.plus.bam $out_dir$filename.fwd1.bam $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.plus.bam

## get minus strand
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam

## now do coverage separately
bedtools genomecov -bga -split -ibam $out_dir$filename.plus.bam > $out_dir$filename.plus.cov.bed
bedtools genomecov -bga -split -ibam $out_dir$filename.minus.bam > $out_dir$filename.minus.cov.bed

## remove all overlap with known exons and introns for coding genes (whole gene) in a strand specific manner. remove poor mapping regions
bedtools subtract -a $out_dir$filename.plus.cov.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/encode.blacklist.nochr.hg38.bed \
  | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/whole.gene.plus.nochr.ucsc.hg38.bed > $out_dir$filename.plus.subtract.whole.gene.bed

bedtools subtract -a $out_dir$filename.minus.cov.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/encode.blacklist.nochr.hg38.bed \
  | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/whole.gene.minus.nochr.ucsc.hg38.bed > $out_dir$filename.minus.subtract.whole.gene.bed

## merge segments that are separated by 500bp (caron) or 1000bp (kapranov) or less
bedtools merge -i $out_dir$filename.plus.subtract.whole.gene.bed -d 1000 > $out_dir$filename.plus.subtract.merge.1kb.bed
bedtools merge -i $out_dir$filename.minus.subtract.whole.gene.bed -d 1000 > $out_dir$filename.minus.subtract.merge.1kb.bed

## now only keep fragments of size 50kb or larger. could consider lowering this.
awk '$3-$2>25000{print $0}' $out_dir$filename.plus.subtract.merge.1kb.bed > $out_dir$filename.plus.subtract.merge.1kb.filter.bed
awk '$3-$2>25000{print $0}' $out_dir$filename.minus.subtract.merge.1kb.bed > $out_dir$filename.minus.subtract.merge.1kb.filter.bed

## now merge anything separated by less than 7kb. could raise this.
## this is the raw data for the results
bedtools merge -i $out_dir$filename.plus.subtract.merge.1kb.filter.bed -d 7000 > $out_dir$filename.plus.vlinc.discovery.bed
bedtools merge -i $out_dir$filename.minus.subtract.merge.1kb.filter.bed -d 7000 > $out_dir$filename.minus.vlinc.discovery.bed

## now create a	new file specifically for genome browser that will be colored by strand 
## this is optional for adding a header
#echo track name=\"vlinc_finder\" description=\"vlincs colored by strand\" visibility=2 colorByStrand=\"255,0,0 0,0,255\" > $out_dir$filename.vlinc.discovery.browser.bed

awk 'OFS="\t"{print "chr"$1,$2,$3,"contig_"NR,0,"+"}' $out_dir$filename.plus.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI > $out_dir$filename.vlinc.discovery.browser.bed
awk 'OFS="\t"{print "chr"$1,$2,$3,"contig_"NR,0,"-"}' $out_dir$filename.minus.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI >> $out_dir$filename.vlinc.discovery.browser.bed



##### remove temp files (will combine steps above for final version making this not necessary)
##### too many temp files will make this slow
##### 

rm $out_dir$filename.mq30.bam
#rm $out_dir$filename.mq30.plus.cov.bed
#rm $out_dir$filename.mq30.minus.cov.bed

rm $out_dir$filename.plus.subtract.whole.gene.bed
rm $out_dir$filename.minus.subtract.whole.gene.bed

rm $out_dir$filename.plus.subtract.merge.1kb.bed
rm $out_dir$filename.minus.subtract.merge.1kb.bed

rm $out_dir$filename.plus.subtract.merge.1kb.filter.bed
rm $out_dir$filename.minus.subtract.merge.1kb.filter.bed

rm $out_dir$filename.fwd1.bam
rm $out_dir$filename.fwd2.bam
rm $out_dir$filename.rev1.bam
rm $out_dir$filename.rev2.bam

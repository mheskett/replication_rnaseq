#!/usr/bin/bash

input=$1 ## bam file from star
out_dir=$2
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

##########
gatk AddOrReplaceReadGroups --java-options "-Xmx30G" -I=$input -O=$out_dir/$filename.rg.sorted.bam \
  --RGID=1 --RGLB=lib1 --RGPL=illumina --RGPU=unit1 --RGSM=$filename -SO=coordinate

##########
gatk MarkDuplicates --java-options "-Xmx30G" -I=$out_dir/$filename.rg.sorted.bam -O=$out_dir/$filename.rg.sorted.markdup.bam \
  --METRICS_FILE=$filename.metrics.txt --ASSUME_SORT_ORDER=coordinate --REMOVE_DUPLICATES=true --CREATE_INDEX=true
#########
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
-T SplitNCigarReads -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa -I $out_dir/$filename.rg.sorted.markdup.bam \
  -o $out_dir/$filename.split.bam \
  -rf ReassignOneMappingQuality \
  -RMQF 255 \
  -RMQT 60 \
  -U ALLOW_N_CIGAR_READs \
  --fix_misencoded_quality_scores

######### gatk3

java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
  -T HaplotypeCaller -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa \
  -I $out_dir/$filename.split.bam --genotyping_mode DISCOVERY \
  -L /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.sorted.noM.vcf \
  --output_mode EMIT_ALL_CONFIDENT_SITES \
  -o $out_dir/$filename.hcgatk3.vcf \
  -stand_call_conf 10.0 \
  -stand_emit_conf 20.0 \
  -dontUseSoftClippedBases

#########

gatk VariantsToTable -V $out_dir/$filename.hcgatk3.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -O $out_dir/$filename.table

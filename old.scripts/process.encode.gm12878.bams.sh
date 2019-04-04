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

######### This is for old files with illumina 1.5 base quality encoding
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
-T SplitNCigarReads -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa -I $out_dir/$filename.rg.sorted.markdup.bam \
  -o $out_dir/$filename.split.bam \
  -rf ReassignOneMappingQuality \
  -RMQF 255 \
  -RMQT 60 \
  -U ALLOW_N_CIGAR_READs \
  --fix_misencoded_quality_scores


######### gatk3 haplotype caller
## caveat is that it wont output sites that are perfectly monoallelic reference which could be many
# should actually use the NA12878 platinum vcf file here if applicable
#   -L /home/groups/Spellmandata/heskett/refs/dbsnp.146.hg38.nochr.sorted.noM.vcf 
# not going to use above, but could make a proper interval list. can make bed file

############# BP resolution mode
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
  -T HaplotypeCaller -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa \
  -I $out_dir/$filename.split.bam --genotyping_mode DISCOVERY \
  -L /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.het.bed \
  -ERC BP_RESOLUTION -o $out_dir/$filename.bp.res.vcf \
  -stand_call_conf 10.0 -stand_emit_conf 20.0 -ip 100 -dontUseSoftClippedBases

#### BP resolution version
### should somehow remove INDELs that don't get called properly?

gatk VariantsToTable -V $out_dir/$filename.bp.res.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -O $out_dir/$filename.table

tail -n +2 $out_dir/$filename.table  | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5,$6}' | grep -Fv \. | grep -v NA > $out_dir/$filename.bed

bedtools intersect -wa -wb -a $out_dir/$filename.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.bed > $out_dir/$filename.overlap.platinum.bed

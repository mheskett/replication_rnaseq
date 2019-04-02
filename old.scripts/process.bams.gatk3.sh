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

### for files with correct illumina encoding
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
-T SplitNCigarReads -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa -I $out_dir/$filename.rg.sorted.markdup.bam \
  -o $out_dir/$filename.split.bam \
  -U ALLOW_N_CIGAR_READS

##### alternatively use samtools pileup
## double check parameters
## only works on het sites
## need script to parse VCF

#bcftools mpileup -f /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa \
#  -R /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.het.bed \
#  --ff SECONDARY -o $out_dir/$filename.pileup.vcf $out_dir/$filename.split.bam

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

#########
#gatk VariantsToTable -V $out_dir/$filename.bp.res.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -O $out_dir/$filename.table

#tail -n +2 $out_dir/$filename.table | awk 'OFS="\t"{split($6,a,",");print $1,$2-1,$2,$3,$4,a[1],a[2]}' | grep -Fv \. | grep -v NA > $out_dir/$filename.bed

#bedtools intersect -wa -wb -a $out_dir/$filename.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.bed > $out_dir/$filename.overlap.platinum.bed

######## python script to arrange the haplotypes

#python /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/haplotyping.py --bed $out_dir/$filename.overlap.platinum.bed --out_directory $out_dir




#### BP resolution version

#gatk VariantsToTable -V $out_dir/$filename.bp.res.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -O $out_dir/$filename.table

tail -n +2 $out_dir/$filename.table | awk '$4!="<NON_REF>"{print $0}' | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5,$6}' | grep -Fv \. | grep -v NA > $out_dir/$filename.bed

bedtools intersect -wa -wb -a $out_dir/$filename.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.bed > $out_dir/$filename.overlap.platinum.bed

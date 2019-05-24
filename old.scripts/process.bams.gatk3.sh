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

############# BP resolution mode for GM12878...
### BP res mode outputs a record for every site, regardless of whether there is a "call"
### if there is a call, gatk will list the allele
### if there is not a call, bu still reads supporting an alt allele, there will be a <NON_REF>
### in gm12878 <NON_REF> alleles are sequencing errors or SNPs

### Use -L dbsnp file from human genome for other
java -Xmx30g -jar /home/groups/Spellmandata/heskett/tools/gatk3.5/GenomeAnalysisTK.jar \
  -T HaplotypeCaller -R /home/groups/Spellmandata/heskett/refs/hg38.10x.nochr.fa \
  -I $out_dir/$filename.split.bam --genotyping_mode DISCOVERY \
  -L /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.het.bed \
  -ERC BP_RESOLUTION -o $out_dir/$filename.bp.res.vcf \
  -stand_call_conf 10.0 -stand_emit_conf 20.0 -ip 100 -dontUseSoftClippedBases

#########
#### BP resolution version

gatk VariantsToTable -V $out_dir/$filename.bp.res.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -O $out_dir/$filename.table

tail -n +2 $out_dir/$filename.table | awk 'OFS="\t"{split($4,a,",");print $1,$2-1,$2,$3,$4,$6}' | awk '$6!="0,0"{print $0}' | awk '$6!="0,0,0"{print $0}' \
  | grep -Fv \. | grep -v NA > $out_dir/$filename.bed

bedtools intersect -wa -wb -a $out_dir/$filename.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/platinum.genome/NA12878.nochr.het.bed > $out_dir/$filename.overlap.platinum.bed

#### now final analyses, output plots, files for circos, etc
#####

python /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/haplotyping.bp.res.py --bed $out_dir/$filename.overlap.platinum.bed --out_directory $out_dir
python /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/allele.rnaseq.py --df $out_dir/$filename.overlap.platinum.haplotype.resolved.bp.res.bed --out_directory $out_dir
/home/groups/Spellmandata/heskett/packages/bin/Rscript /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/test_dnacopy.r $out_dir/$filename.overlap.platinum.haplotype.resolved.bp.res.dnacopy.txt



# remove intermediate files
rm $out_dir/$filename.table
rm $out_dir/$filename.rg.sorted.markdup.bam
rm $out_dir/$filename.rg.sorted.bam
rm $out_dir/$filename.bed

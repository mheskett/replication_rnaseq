#!/usr/bin/bash
#Michael Heskett
## UASAGE:  ./variant_caller.sh normalBAM(full_path) tumorBAM(full_path) /path/to/outputdirectory/
#libs_path='/home/exacloud/lustre1/SpellmanLab/heskett/miniconda3/bin'

normalBAM=$1
n=$(basename $normalBAM) ## removes /path/to/file
normalfname=${n%.*} ### removes file extension
echo $normalfname

tumorBAM=$2
t=$(basename $tumorBAM)
tumorfname=${t%.*}
echo $tumorfname

DATE=`date +%Y-%m-%d`
random=$RANDOM

################
## Make subdirectories within a bigger directory
mkdir $3$tumorfname
output_dir=$3$tumorfname/
echo $output_dir
################


### ALL REFERENCES HAVE 'CHR' PREFIX
#ref='/home/exacloud/lustre1/SpellmanLab/heskett/teratoma/refs/genome.fa'
#indel='/home/exacloud/lustre1/SpellmanLab/heskett/asia_references/1000G_phase1.indels.hg19.sites.vcf'
#dbsnp='/home/exacloud/lustre1/SpellmanLab/heskett/asia_references/dbsnp_138.hg19.vcf'
#cosmic='/home/exacloud/lustre1/SpellmanLab/heskett/asia_references/CosmicCodingMuts.prefix.sorted.vcf'

# references lack chr prefix
ref="/home/groups/Spellmandata/heskett/myron_refs/human_g1k_v37.fasta"
###indel="/home/groups/Spellmandata/heskett/myron_refs/1000g_indels.vcf"
dbsnp="/home/exacloud/tempwork/SpellmanLab/heskett/DBSNP_150_uw.vcf"
cosmic="/home/groups/Spellmandata/heskett/refs/hg19/CosmicCodingMuts.vcf"
#######

snvVCF=$output_dir$tumorfname'.snv.vcf'
callstats=$output_dir$tumorfname'.snv.out'

/usr/lib/jvm/java-1.7.0/jre/bin/java -jar -Xmx16g /home/groups/Spellmandata/atlas/mutation_pipeline/tools/muTect/muTect-1.1.5.jar \
                -T MuTect \
                -I:normal $normalBAM \
                -I:tumor $tumorBAM  \
                -R $ref \
                --dbsnp $dbsnp \
                -L /home/groups/Spellmandata/heskett/teratoma/gatk_scripts/idt_v1.0_exome_hg19_nochr_padding100.bed \
                --cosmic $cosmic \
                --out $callstats \
                --coverage_file $output_dir$tumorfname.coverage.wig.txt \
                -vcf $snvVCF
################
## STEP 8: Filter and annotate variants
## Remove SNPs, rejects, and uncovered sites.
## Annotate filtered variants with ANNOVAR and remove synonymous SNVs
################

#grep -v '\sDBSNP\s' $callstats | grep -vwE "REJECT" | grep -vwE "UNCOVERED" > $callstats.NOVEL.COVERED.KEEP
#awk -F'\t' '{print($1,$2,$2,$4,$5);}' $callstats.NOVEL.COVERED.KEEP > $callstats.ann
#perl /home/exacloud/lustre1/SpellmanLab/heskett/prostate_neoantigen/annovar/annotate_variation.pl -geneanno \
#                        -buildver hg19 \
#                        $callstats.ann \
#                        /home/exacloud/lustre1/SpellmanLab/heskett/prostate_neoantigen/annovar/humandb/ \
#                        -outfile $callstats
#grep -vwE "synonymous SNV" $callstats.exonic_variant_function >> $callstats.nonsyn

#grep -v REJECT $snvVCF > $output_dir$tumorfname'.filtered.vcf'

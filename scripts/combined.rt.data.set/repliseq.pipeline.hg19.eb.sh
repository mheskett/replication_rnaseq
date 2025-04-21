## align haplotypes with eb script
python align.haplotypes.py --bed  eb3_2_clone10_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone10_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone13_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone13_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone15_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone15_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone2_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone2_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone3_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone3_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone4_early_rt.sorted.markdup.allele.counts.bed --out_directory ./
python align.haplotypes.py --bed  eb3_2_clone4_late_rt.sorted.markdup.allele.counts.bed  --out_directory ./

## remove blacklisted snps before LSM hg19
./remove.blacklist.sh eb3_2_clone10_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone10_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone13_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone13_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone15_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone15_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone2_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone2_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone3_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone3_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone4_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed
./remove.blacklist.sh eb3_2_clone4_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.bed

# do LSM after the bad snps are removed
python lsm.normalize.allele.counts.py eb3_2_clone10_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone10_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone13_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone13_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone15_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone15_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone2_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone2_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone3_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone3_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone4_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed
python lsm.normalize.allele.counts.py eb3_2_clone4_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.bed

## overlap with windows. these windows also have the blacklist subregions removed
./make.allele.window.counts.hg19.sh eb3_2_clone10_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone10_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone13_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone13_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone15_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone15_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone2_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone2_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone3_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone3_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone4_early_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed
./make.allele.window.counts.hg19.sh eb3_2_clone4_late_rt.sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.bed

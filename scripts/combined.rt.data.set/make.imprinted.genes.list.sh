### hg38
awk '{print "grep -w \"" $1 "\" ucsc.refseq.hg38.txn.whole.gene.bed >> hg38.imprinted.genes.txt"}' imprinted.high.confidence.genes.txt > get.imprinted.loci.sh 
./get.imprinted.loci.sh

### hg19
awk '{print "grep -w \"" $1 "\" ucsc.refseq.hg19.txn.whole.gene.sorted.nochr.bed >> hg19.imprinted.genes.txt"}' imprinted.high.confidence.genes.txt > get.imprinted.loci.hg19.sh
./get.imprinted.loci.hg19.sh


library(DNAcopy)
test = read.table("/Users/heskett/replication_rnaseq/scripts/gm12878.rep1Aligned.hcgatk3.overlap.platinum.haplotype.resolv.dnacopy.txt",header=TRUE)
test_1 = test[test$chrom=="1",]

CNA.object = CNA(cbind(test_1$logR_dnacopy),test_1$chrom,test_1$stop,data.type="logratio",sampleid="test1")
#smoothed.CNA.object = smooth.CNA(CNA.object)

#segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
segment.CNA.object <- segment(CNA.object, verbose=1,alpha=0.0001) # alpha is significance level to accept change points.
## in general i believe higher alpha means more segments--try alpha between 10**-8 and 5*10**-2
## try using weights based on read depth?
#plot(segment.smoothed.CNA.object, plot.type="w")
plot(segment.CNA.object)


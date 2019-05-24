library(DNAcopy)

arguments = commandArgs(trailingOnly=TRUE)
df = read.table(arguments[1],header=TRUE)

CNA.object = CNA(cbind(df$logR_dnacopy), factor(df$chrom,levels=c(1:22,"X")), df$stop,
	data.type="logratio",
	sampleid="sample1")
smoothed.CNA.object = smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1, alpha=0.0001)
## in general i believe higher alpha means more segments--try alpha between 10**-8 and 5*10**-2
## try using weights based on read depth?
#plot(segment.smoothed.CNA.object, plot.type="w")
pdf(paste(arguments[1],".pdf",sep=""),height=4,width=12)
plotSample(segment.smoothed.CNA.object)
dev.off()

pdf(paste(arguments[1],"allchr.pdf",sep=""),height=9,width=20)
plot(segment.smoothed.CNA.object,plot.type='s')
dev.off()
## now need a cutoff for segment mean to indicate monoallelic

write.table(segment.smoothed.CNA.object$output,
	file=paste(tools::file_path_sans_ext(arguments[1]),".segments.txt",sep=''),
	quote=FALSE,
	sep="\t")

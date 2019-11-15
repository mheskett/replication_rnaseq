library(EMT)

df = read.table("/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed",
				sep="\t",
				col.names=c("chrom","start","stop","line_name","score","strand")
				)
relative_frequencies = table(df$line_name) / nrow(df)
print(relative_frequencies["L1PA2"]) ##kk
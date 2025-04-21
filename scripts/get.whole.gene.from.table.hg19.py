import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

"""
example. this file has chr
chr1    66999043        67216822        SGIP1   0       +
chr1    66999043        67216822        SGIP1   0       +
chr1    66999275        67216822        SGIP1   0       +
chr1    66999275        67216822        SGIP1   0       +
chr1    66999275        67216822        SGIP1   0       +
chr1    66999275        67216822        SGIP1   0       +

"""
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

df = pd.read_csv("ucsc.refseq.hg19.whole.gene.bed",sep="\t",
	names=["chrom","start","stop","gene","score","strand"])

df = df[df["chrom"].isin(chromosomes)]

genes = df["gene"].unique()

result = []
count=0
for gene in genes:
	count+=1
	tmp = df[df["gene"]==gene]
	chrom = tmp["chrom"].unique()[0]
	start=tmp["start"].min()
	stop=tmp["stop"].max()
	strand=tmp["strand"].unique()[0]

	# if len(chrom)>1:
	# 	print(chrom)
	# 	print(gene)
	# 	continue

	# if len(strand)>1:
	# 	print(strand)
	# 	print(gene)
	# print(count)
	result+=[[chrom,start,stop,gene,0,strand]]
# print(result)

df_final = pd.DataFrame(result,columns=["chrom","start","stop","gene","score","strand"])
df_final = df_final.sort_values(by=["chrom","start"])
df_final.to_csv("ucsc.refseq.hg19.txn.whole.gene.bed",sep="\t",header=None,index=None)
os.system("sort -k1,1 -k2,2n ucsc.refseq.hg19.txn.whole.gene.bed > ucsc.refseq.hg19.txn.whole.gene.sorted.bed")
os.system("sort -k1,1 -k2,2n ucsc.refseq.hg19.txn.whole.gene.bed | sed 's/chr//g' > ucsc.refseq.hg19.txn.whole.gene.sorted.nochr.bed")
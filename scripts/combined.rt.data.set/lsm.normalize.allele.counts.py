import numpy as np
import pandas as pd
from sys import argv
## file format is chrom start stop hap1 hap2 hap1_counts hap2_counts

df = pd.read_csv(argv[1],sep="\t",dtype={"chrom":str},
	names=["chrom", "start", "stop" ,"hap1", "hap2", "hap1_counts", "hap2_counts"])
total_reads = (df["hap1_counts"]+df["hap2_counts"]).sum()
df["hap1_lsm"]=(df["hap1_counts"]/total_reads)*10**6
df["hap2_lsm"]=(df["hap2_counts"]/total_reads)*10**6

df.loc[:,["chrom", "start", "stop" ,"hap1", "hap2", "hap1_lsm", "hap2_lsm"]]\
	.to_csv(argv[1].removesuffix(".bed")+".cpms.bed",header=None,index=False,sep="\t")

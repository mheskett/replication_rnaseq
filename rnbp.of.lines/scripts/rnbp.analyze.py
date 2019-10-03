import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns

# for the cluster
# eclip = "/home/groups/Spellmandata/heskett/replication.rnaseq/rnbp.of.lines/big.download/all.eclip.nochr.sorted.bed"
# introns =  "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
# vlincs = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
# lines = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

## for local use
eclip_file = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.hepg2.chr6.nochr.sorted.bed"
introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs_file = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"
eclip = pybedtools.BedTool(eclip_file) # change this to bed file of previously determined windows of interest
introns = pybedtools.BedTool(introns_file)
vlincs = pybedtools.BedTool(vlincs_file)
lines = pybedtools.BedTool(lines_file)

#### make dict of vlinc length
f = open(vlincs_file, 'r')
vlinc_lengths_dict = {}
for line in f:
    tmp = line.strip().split('\t')
    vlinc_lengths_dict[tmp[3]] = int(tmp[2])-int(tmp[1])

f.close()
################
df_eclip = pd.read_table(eclip_file,
					names=["chrom", "start", "stop", "name","score", 
					"strand","signalvalue","pvalue","qvalue","peak"],
					 dtype={"chrom":str, "start":int, "stop":int,
					 "name":str, "score":float, "strand":str,
					 "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})
###############
df_eclip = df_eclip[(~df_eclip.name.str.contains("IDR",case=False))] # filter out replication experiments
eclip_library_size = df_eclip.name.value_counts() ## original library size
print("number experiments: ",len(eclip_library_size))
print("number experiments not including replicates: ",df_eclip.name.str.replace(pat="_rep0[1-2]",repl="",regex=True).nunique())
################

### make data frame, sort by bindings per lncRNA
eclip_within_vlincs_raw_df = eclip.intersect(vlincs,
										f=0.9,
										s=True,
										wa=True,wb=True)\
										.to_dataframe(names=["chrom", "start", "stop", "name_eclip", 
											 "score_eclip", "strand_eclip","signalvalue","pvalue","qvalue","peak",
											 "chrom_lnc","start_lnc","stop_lnc","name_lnc","score_lnc","strand_lnc"],
											 dtype={"chrom":str, "start":int, "stop":int, "name_eclip":str, 
											 "score_eclip":int, "strand_eclip":str,"signalvalue":float,"pvalue":float,"qvalue":float,"peak":float,
											 "chrom_lnc":str,"start_lnc":int,"stop_lnc":int,"name_lnc":str,"score_lnc":int,"strand_lnc":str})
#### filter out insignificant peaks
eclip_within_vlincs_raw_df = eclip_within_vlincs_raw_df[(~eclip_within_vlincs_raw_df.name_eclip.str.contains("IDR",case=False)) & (eclip_within_vlincs_raw_df["pvalue"] >= 2) & (eclip_within_vlincs_raw_df["signalvalue"] >= 1)]

#### turn this into a samples/features matrix
eclip_within_vlincs_df = eclip_within_vlincs_raw_df.groupby("name_lnc")["name_eclip"].value_counts().to_frame().rename(columns={"name_eclip":"counts"}).reset_index().pivot(index="name_lnc",columns="name_eclip",values="counts").fillna(value=0)

################## normalize to significant peaks/100 peaks
eclip_within_vlincs_normalized_df = eclip_within_vlincs_df.divide(eclip_library_size/100, axis="columns")

################## filter the matrix
vlinc_min_binding = np.percentile(eclip_within_vlincs_normalized_df.sum(axis="columns"),25) 
rnbp_min_binding = np.percentile(eclip_within_vlincs_normalized_df.sum(axis="index"),25)
eclip_within_vlincs_normalized_df = eclip_within_vlincs_normalized_df.loc[eclip_within_vlincs_normalized_df.sum(axis="columns") >= vlinc_min_binding, :] # filters vlincs ## k562 value 0.02 ## hepg2 value 0.003
eclip_within_vlincs_normalized_df = eclip_within_vlincs_normalized_df.loc[:, eclip_within_vlincs_normalized_df.sum(axis="index") >= rnbp_min_binding] # filters proteins ## k562 value 0.34 ## hepg2 value 0.14

############### check to see if the results correlate with 1) length of ASAR, 2) size of eclip library
final_proteins = eclip_within_vlincs_normalized_df.columns
eclip_within_vlincs_normalized_df.loc[:,"col_sum"] = eclip_within_vlincs_normalized_df.sum(axis="columns")
eclip_within_vlincs_normalized_df.loc[:,"vlinc_length"] = eclip_within_vlincs_normalized_df.apply(lambda x: vlinc_lengths_dict[x.name], axis=1)
eclip_library_size_final = eclip_library_size[final_proteins]
### now do correlations
print("correlation between feature length and feature binding: ",scipy.stats.spearmanr(a=eclip_within_vlincs_normalized_df["vlinc_length"], 
							b=eclip_within_vlincs_normalized_df["col_sum"])) ## correlation = 0.16 p =.001 not big.
print("correlation between rnbp library size and feature binding: ", scipy.stats.spearmanr(a=eclip_within_vlincs_normalized_df.drop(["col_sum","vlinc_length"],axis=1).sum(axis="index"),
							b=eclip_library_size_final )) ## SpearmanrResult(correlation=-0.17381298381123095, pvalue=0.08061741938032327)
####### output
out1 = os.path.basename(eclip_file.replace(".bed",".most.bound.lncrnas.bed"))
out2 = os.path.basename(eclip_file.replace(".bed",".most.bound.rnbp.txt"))

##### lncrnas
#### stupid trouble here...
df_vlinc = pd.read_table(vlincs_file,names=["chrom","start","stop","name","score","strand"])
print("n unique names ", df_vlinc["name"].value_counts())

df_vlinc = df_vlinc.set_index("name")
# print(df_vlinc)
df_vlinc = df_vlinc.loc[eclip_within_vlincs_normalized_df.index.values,:] ## this line  isn't working how you think it will. NVM
## found bug. vlinc database has repeat names
# print(len(eclip_within_vlincs_normalized_df.index))
# print(df_vlinc)
print(eclip_within_vlincs_normalized_df)
# print(eclip_within_vlincs_normalized_df["col_sum"])
print("null values: ",eclip_within_vlincs_normalized_df.isna().sum().sum()) # zero null values....

out = pd.concat([df_vlinc, eclip_within_vlincs_normalized_df["col_sum"]],axis=1)
print(out)
out.sort_values("col_sum").to_csv(out1,sep="\t")
####

eclip_within_vlincs_normalized_df.drop(["col_sum","vlinc_length"],axis=1).sum(axis="index").sort_values(ascending=False).to_csv(out2, sep="\t")
plt.figure(figsize=(10,2))
f = sns.distplot(eclip_within_vlincs_normalized_df.drop(["col_sum","vlinc_length"],axis=1).sum(axis="columns"),rug=True)
f.set(xlim=(0,None))
f.set_title("Peaks per hundred - feature")
plt.savefig(out1.replace(".bed",".png"))
plt.close()
plt.figure(figsize=(10,2))
f = sns.distplot(eclip_within_vlincs_normalized_df.drop(["col_sum","vlinc_length"],axis=1).sum(axis="index"),rug=True)
f.set_title("Peaks per hundred - RNA binding proteinls -alt")
f.set(xlim=(0,None))
plt.savefig(out2.replace(".txt",".png"))

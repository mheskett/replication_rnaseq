import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools

# for the cluster
# eclip = "/home/groups/Spellmandata/heskett/replication.rnaseq/rnbp.of.lines/big.download/all.eclip.nochr.sorted.bed"
# introns =  "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
# vlincs = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
# lines = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

## for local use
eclip_file = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.k562.chr6.nochr.sorted.bed"
introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs_file = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

eclip = pybedtools.BedTool(eclip_file) # change this to bed file of previously determined windows of interest
introns = pybedtools.BedTool(introns_file)
vlincs = pybedtools.BedTool(vlincs_file)
lines = pybedtools.BedTool(lines_file)


# print("Total number L1s within all vlincs ",len(lines.intersect(vlincs, f=0.9)))
# print("Total number L1s within all introns ", len(lines.intersect(introns,f=0.9)))


### eclip database
### want to get the total number of peaks per protein
df_eclip = pd.read_table("/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.k562.chr6.nochr.sorted.bed",
					names=["chrom", "start", "stop", "name","score", "strand","signalvalue","pvalue","qvalue","peak"],
					 dtype={"chrom":str, "start":int, "stop":int,
					 "name":str, "score":float, "strand":str,
					 "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})


df_eclip = df_eclip[(~df_eclip.name.str.contains("IDR",case=False))]
print(df_eclip.name.value_counts())

## could prepare these two files on command line so i dont have to redo it during development 
eclip_within_vlincs = eclip.intersect(vlincs,
										f=0.9,
										s=True,
										wa=True)


eclip_within_introns = eclip.intersect(introns,
										f=0.9,
										s=True,
										wa=True)
print("done first step")
eclip_bind_vlinc_line = lines.intersect(eclip_within_vlincs, F=0.9, wa=True, wb=True).to_dataframe(names=["chrom", "start", "stop", "name", 
																		 "score", "strand","chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip",
																		 "strand_eclip","signalvalue","pvalue","qvalue","peak"],
																		 dtype={"chrom":str, "start":int, "stop":int,
																		 "name":str, "score":float, "strand":str,
																		 "chrom_eclip":str,"start_eclip":int,"stop_eclip":int,"name_eclip":str,"score_eclip":int,
																		 "strand_eclip":str, "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})

eclip_bind_intron_line = lines.intersect(eclip_within_introns, F=0.9, wa=True, wb=True).to_dataframe(names=["chrom", "start", "stop", "name", 
																		 "score", "strand","chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip",
																		 "strand_eclip","signalvalue","pvalue","qvalue","peak"],
																		 dtype={"chrom":str, "start":int, "stop":int,
																		 "name":str, "score":float, "strand":str,
																		 "chrom_eclip":str,"start_eclip":int,"stop_eclip":int,"name_eclip":str,"score_eclip":int,
																		 "strand_eclip":str, "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})

print("done second step")
eclip_bind_vlinc_line.loc[:,"unique_line_name"] = eclip_bind_vlinc_line.apply(lambda row: str(row["chrom"]) + ":" + str(row["start"]) + "-" + str(row["stop"]),axis=1)
eclip_bind_intron_line.loc[:,"unique_line_name"] = eclip_bind_intron_line.apply(lambda row: str(row["chrom"]) + ":" + str(row["start"]) + "-" + str(row["stop"]),axis=1)

print("total number eclip binding events within vlinc lines: ",len(eclip_bind_vlinc_line))
print("total number eclip binding events within intron lines: ",len(eclip_bind_intron_line))

### just print out some stats
print("total experiments")
print(len(list(eclip_bind_intron_line["name_eclip"].unique())))
print("total different proteins not including replicates")
print(len(list(eclip_bind_intron_line["name_eclip"].str.replace(pat="_rep0[1-2]",repl="",regex=True).unique())))


### i can gather that signal value is log2 (fold enrichment), and pvalue is -log10 pvalue of the enrichment
eclip_bind_vlinc_line = eclip_bind_vlinc_line[(~eclip_bind_vlinc_line.name_eclip.str.contains("IDR",case=False)) 
																& (eclip_bind_vlinc_line["pvalue"] >= 2) 
																& (eclip_bind_vlinc_line["signalvalue"] >= 1)]

eclip_bind_intron_line = eclip_bind_intron_line[(~eclip_bind_intron_line.name_eclip.str.contains("IDR",case=False)) 
																& (eclip_bind_intron_line["pvalue"] >= 2) 
																& (eclip_bind_intron_line["signalvalue"] >= 1)]
########

print("L1s bound in vlincs ", eclip_bind_vlinc_line["unique_line_name"].nunique())
print("L1s bound in introns ",eclip_bind_intron_line["unique_line_name"].nunique())

### value counts gives a counts matrix as a series. to frame makes it a multi index frame. rename gets rid of duplicated column name error. reset index drops multi indexing
### pivot makes rows by columns counts matrix.
df_intron = eclip_bind_intron_line.groupby("unique_line_name").name_eclip.value_counts().to_frame().rename(columns={"name_eclip":"counts"}).reset_index().pivot(index="unique_line_name",columns="name_eclip",values="counts").fillna(value=0)
df_vlinc = eclip_bind_vlinc_line.groupby("unique_line_name").name_eclip.value_counts().to_frame().rename(columns={"name_eclip":"counts"}).reset_index().pivot(index="unique_line_name",columns="name_eclip",values="counts").fillna(value=0)

print("total L1 peaks in introns ", df_intron.sum(axis="columns").sum())
print("total L1 peaks in vlincs ", df_vlinc.sum(axis="columns").sum())



# unique_lines = list(eclip_bind_intron_line["unique_line_name"].unique())
# unique_eclip = list(eclip_bind_intron_line["name_eclip"].unique())
# results = np.zeros( (len(unique_lines), len(unique_eclip) ) ) 



# # for i in range(len(unique_lines)):
# # 	df_unique_line = eclip_binding_lines_in_introns[eclip_binding_lines_in_introns["unique_line_name"]==unique_lines[i]]
# # 	print(df_unique_line["name_eclip"].value_counts())

# # fraction_of_intronic_lines_bound = eclip_binding_lines_in_introns.groupby(["name_eclip"]).size() / total_intronic_l1s
# # fraction_of_vlinc_lines_bound = eclip_binding_lines_in_vlincs.groupby(["name_eclip"]).size() / total_vlinc_l1s

# # print(fraction_of_intronic_lines_bound)
# # print(fraction_of_vlinc_lines_bound)
# # print(eclip_binding_lines_in_introns[eclip_binding_lines_in_introns["unique_line_name"]=="6:534111-534343"]["name_eclip"].value_counts())
# print(eclip_binding_lines_in_introns.groupby("unique_line_name")["name_eclip"].value_counts())

# # fraction_of_vlinc_lines_bound.divide(fraction_of_intronic_lines_bound).to_csv("vlinc_vs_intron_l1_rnbp_binding.txt",sep="\t")





## goal is to get binding sites per LINE for each protein for set of intronic lines and vlinc lines
## 
# print(eclip_binding_lines_in_vlincs.groupby)
# eclip_binding_lines_in_introns = eclip_binding_lines_in_introns[eclip_binding_lines_in_introns["count"]>0]
# eclip_binding_lines_in_vlincs = eclip_binding_lines_in_vlincs[eclip_binding_lines_in_vlincs["count"]>0]

# print(eclip_binding_lines_in_introns.groupby(["name"])["count"].mean())

# df = pd.read_table("/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.chr6.nochr.sorted.bed",
# 	sep="\t",names=["chrom","start","stop","name","score","strand","signal_value","pvalue","qvalue","peak"],
# 	index_col=None)

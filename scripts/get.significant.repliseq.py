import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
import seaborn as sns
import scipy.stats
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
import statsmodels.api as sm
"""
Repliseq directions ( do this separately for early and late)

1. take allele counts bed file and sum over 100kb windows that are tiled by 25kb -- 
2. get autosome wide distribution of logR values
3. only keep top 10% most asynchronous
4. Merge the segments, requiring at least 300kb contiguous with allowable gap of 50kb
4. Merge again, allowing for a 100kb gap.
5. Compare overlap between early and late significant regions. get p-values
"""

ratios=[249250621/249250621,
	243199373/249250621,
	198022430/249250621,
	191154276/249250621,
	180915260/249250621,
	171115067/249250621,
	159138663/249250621,
	146364022/249250621,
	141213431/249250621,
	135534747/249250621,
	135006516/249250621,
	133851895/249250621,
	115169878/249250621,
	107349540/249250621,
	102531392/249250621,
	90354753/249250621,
	81195210/249250621,
	78077248/249250621,
	59128983/249250621,
	63025520/249250621,
	48129895/249250621,
	51304566/249250621,
	155270560/249250621]

lengths = [249250621,
	243199373,
	198022430,
	191154276,
	180915260,
	171115067,
	159138663,
	146364022,
	141213431,
	135534747,
	135006516,
	133851895,
	115169878,
	107349540,
	102531392,
	90354753,
	81195210,
	78077248,
	59128983,
	63025520,
	48129895,
	51304566,
	155270560]

centromere = {"1":124535434,
				"2":95326171,
				"3":93504854,
				"4":52660117,
				"5":49405641,
				"6":61830166,
				"7":61054331,
				"8":46838887,
				"9":50367679,
				"X":61632012,
				"Y":13104553,
				"10":42254935,
				"11":54644205,
				"12":37856694,
				"13":19000000,
				"14":19000000,
				"15":20000000,
				"16":38335801,
				"17":25263006,
				"18":18460898,
				"19":27681782,
				"20":29369569,
				"21":14288129,
				"22":16000000}
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/data/cytoband.nochr.hg19.bed",
							sep="\t",
							names =["chrom","start","stop","arm","band"])

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

arm_dict = get_arms(cytoband)

df_windows = pd.read_csv("4e.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df2_windows = pd.read_csv("4l.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df_windows = df_windows[df_windows["hap1_counts"] + df_windows["hap2_counts"] >= 15]
df2_windows = df2_windows[df2_windows["hap1_counts"] + df2_windows["hap2_counts"] >= 15]

df_windows.loc[:,"logR"] = np.log2( (df_windows["hap1_counts"]+1) / (df_windows["hap2_counts"]+1) )

df2_windows.loc[:,"logR"] = np.log2( (df2_windows["hap1_counts"]+1) / (df2_windows["hap2_counts"]+1) )

df_windows["arm"] = df_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df2_windows["arm"] = df2_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

# get thresholds top 5% skewed on all autosomes.
early_thresholds = np.percentile(df_windows[df_windows["chrom"] != "X"]["logR"],[5,95])
late_thresholds = np.percentile(df2_windows[df2_windows["chrom"] != "X"]["logR"],[5,95])
##
asynch_early_hap1 =  df_windows[(df_windows["logR"] >= early_thresholds[1])]
asynch_early_hap2 = df_windows[(df_windows["logR"] <= early_thresholds[0])]

asynch_late_hap1 =  df2_windows[(df2_windows["logR"] >= late_thresholds[1])]
asynch_late_hap2 = df2_windows[(df2_windows["logR"] <= late_thresholds[0])]

##
asynch_early_hap1_bed = pybedtools.BedTool.from_dataframe(asynch_early_hap1)
asynch_early_hap2_bed = pybedtools.BedTool.from_dataframe(asynch_early_hap2)

asynch_late_hap1_bed = pybedtools.BedTool.from_dataframe(asynch_late_hap1)
asynch_late_hap2_bed = pybedtools.BedTool.from_dataframe(asynch_late_hap2)

###
early_merged_hap1 = asynch_early_hap1_bed.merge(d=25000,c=[4,5],o=["sum","sum"]).filter(lambda x: len(x) >= 250000).merge(d=100000,c=[4,5],o=["sum","sum"])\
			.to_dataframe(names=["chrom","start","stop","hap1","hap2"])
early_merged_hap2 = asynch_early_hap2_bed.merge(d=25000,c=[4,5],o=["sum","sum"]).filter(lambda x: len(x) >= 250000).merge(d=100000,c=[4,5],o=["sum","sum"])\
			.to_dataframe(names=["chrom","start","stop","hap1","hap2"])


late_merged_hap1 = asynch_late_hap1_bed.merge(d=25000,c=[4,5],o=["sum","sum"]).filter(lambda x: len(x) >= 250000).merge(d=100000,c=[4,5],o=["sum","sum"])\
			.to_dataframe(names=["chrom","start","stop","hap1","hap2"])
late_merged_hap2 = asynch_late_hap2_bed.merge(d=25000,c=[4,5],o=["sum","sum"]).filter(lambda x: len(x) >= 250000).merge(d=100000,c=[4,5],o=["sum","sum"])\
			.to_dataframe(names=["chrom","start","stop","hap1","hap2"])

print(early_merged_hap1)
print(early_merged_hap2)

print(late_merged_hap1)
print(late_merged_hap2)


for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))

	ax.scatter(df_windows[df_windows["chrom"]==chromosomes[i]]["start"],df_windows[df_windows["chrom"]==chromosomes[i]]["logR"],s=5,lw=0.2)
	ax.scatter(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"],df2_windows[df2_windows["chrom"]==chromosomes[i]]["logR"],s=5,lw=0.2)

	ax.axhline(y=0)
	for index, row in early_merged_hap1[early_merged_hap1["chrom"]==chromosomes[i]].iterrows():
		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="blue", alpha=0.5)
	for index, row in early_merged_hap2[early_merged_hap2["chrom"]==chromosomes[i]].iterrows():
		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="red", alpha=0.5)
	for index, row in late_merged_hap1[late_merged_hap1["chrom"]==chromosomes[i]].iterrows():
		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="red", alpha=0.5)
	for index, row in late_merged_hap2[late_merged_hap2["chrom"]==chromosomes[i]].iterrows():
		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="blue", alpha=0.5)
	plt.show()
	plt.close()


# ax.scatter(df2_windows[df2_windows["chrom"]=="1"]["start"],df2_windows[df2_windows["chrom"]=="1"]["logR"],s=5,lw=0.2)
# ax.axhline(y=0)
# for index, row in late_merged[late_merged["chrom"]=="1"].iterrows():
# 	ax.axvspan(xmin=row["start"],xmax=row["stop"],facecolor="red",alpha=0.5)

plt.show()
plt.close()

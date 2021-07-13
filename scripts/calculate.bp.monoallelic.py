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
import glob
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
import pybedtools

def sum_bases(df):
	length = df["stop"] - df["start"]

	return length.sum()
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
chromosome_length = {"1":249250621,
"2":243199373,
"3":198022430,
"4":191154276,
"5":180915260,
"6":171115067,
"7":159138663,
"8":146364022,
"9":141213431,
"10":135534747,
"11":135006516,
"12":133851895,
"13":115169878,
"14":107349540,
"15":102531392,
"16":90354753,
"17":81195210,
"18":78077248,
"19":59128983,
"20":63025520,
"21":48129895,
"22":51304566,
"X":155270560}

def add_binom_pval(df):
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
							row["hap1_counts"]+row["hap2_counts"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)
	results = mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")
	df["fdr_pval"] = results[1]
	df["fdr_reject"] = results[0]

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
def helper_func(x):
	if x["total_reads"]==0: # try this for filtering
		return 0
	elif x["hap1_counts"] >= x["hap2_counts"]:
		return x["hap1_counts"]  / x["total_reads"] - 0.5
	else:
		return -x["hap2_counts"]  / x["total_reads"] + 0.5
	return

arm_dict = get_arms(cytoband)
## coding only
coding_files=["bouha2.protein.coding.all.counts.bed",
"bouha3.protein.coding.all.counts.bed",
"bouha4.protein.coding.all.counts.bed",
"bouha10.protein.coding.all.counts.bed",
"bouha13.protein.coding.all.counts.bed",
"bouha15.protein.coding.all.counts.bed"]
dfs = []
for i in range(len(coding_files)):
	df = pd.read_csv(coding_files[i],sep="\t",
							names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
	df["skew"] = df.apply(helper_func, axis = 1)
	df["sample"] = coding_files[i][0:15]
	add_binom_pval(df)
	dfs += [df]
coding_df = pd.concat(dfs)
coding_df = coding_df[coding_df["total_reads"]>=15]

print(coding_df)
## vlincs only
rna_files=["/Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed"]
dfs = []
for i in range(len(rna_files)):
	df = pd.read_csv(rna_files[i],sep="\t",
							names= ["chrom","start","stop","name","rpkm","strand","l1_density","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
	df["skew"] = df.apply(helper_func, axis = 1)
	df["sample"] = os.path.basename(rna_files[i])[0:15]
	add_binom_pval(df)
	dfs += [df]
df = pd.concat(dfs)
df = df[df["total_reads"]>=15]
print("number TLs in gm12878",len(df))
f,ax=plt.subplots(1,1,figsize=(2.5,4))
ax.scatter(np.log2(df["total_reads"]),abs(df["skew"]),c="red",s=20,lw=0.2,edgecolor="black")
ax.scatter(np.log2(df[df["chrom"]!="X"]["total_reads"]),abs(df[df["chrom"]!="X"]["skew"]),c="mediumblue",s=20,lw=0.2,edgecolor="black")
ax.set_ylim([0,0.5])
# ax.set_xticks([])
plt.savefig("fig1.scatter.png",
			dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.show()
#################
## fig 1c 
for i in range(len(chromosomes)):
	f,ax=plt.subplots(1,1,figsize=(10,2))
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_ylim([-0.52,0.52])
	ax.set_xticks(np.linspace(start,stop, 16))
	ax.set_yticks([-0.5,-.25,0,.25,.5])
	print(np.linspace(start,stop,4))
## vlincs########################
	vlincs_tmp = df_vlinc[(df_vlinc["chrom"]==chrom) & (df_vlinc["start"]>=start-2000000) & 
					(df_vlinc["stop"]<=stop+2000000)]
	########################## vlincs
	ax.axhline(y=0,linestyle="--",c="black",lw=0.2)
	print(vlincs_tmp)
	for index,row in vlincs_tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
		ax.add_patch(rect)
	plt.savefig("gm12878.rep1.protein.vlinc."+str(chrom)+"."+str(start)+"."+str(stop)+".png",
	dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()





#######
#######
tmp = pybedtools.BedTool.from_dataframe(df[(df["fdr_reject"]==True) & (abs(df["skew"])>=0.1)].loc[:,["chrom","start","stop","name","rpkm","strand"]])
merged_df = tmp.sort().merge().to_dataframe()
merged_df.columns = ["chrom","start","stop"]

tmp2 = pybedtools.BedTool.from_dataframe(coding_df[(coding_df["fdr_reject"]==True) & (abs(coding_df["skew"])>=0.1)].loc[:,["chrom","start","stop","name","score","strand"]])
merged_coding_df = tmp2.sort().merge().to_dataframe()
merged_coding_df.columns = ["chrom","start","stop"]
print("number of ade vlincs",len(merged_df))
print("total bases noncoding monoallelic",sum_bases(merged_df))
print("fraction of a haploid genome noncoding monoallelic", sum_bases(merged_df)/(3*10**9))

print("number of aDE coding regions",len(merged_coding_df))
print("total bases coding monoallelic",sum_bases(merged_coding_df))
print("fraction of a haploid genome coding monoallelic", sum_bases(merged_coding_df)/(3*10**9))
# sns.set(rc={'figure.figsize':(4,4)})
fig, ax = plt.subplots(figsize=(2.5,2))
# sns.kdeplot(abs(coding_df[coding_df["chrom"]!="X"]["skew"]),cut=0,lw=4,c="blue",ax=ax)
# sns.kdeplot(abs(coding_df[coding_df["chrom"]=="X"]["skew"]),cut=0,lw=4,c="red",ax=ax)
sns.kdeplot(abs(df[df["chrom"]!="X"]["skew"]),cut=0,lw=4,c="mediumblue",ax=ax)
sns.kdeplot(abs(df[df["chrom"]=="X"]["skew"]),cut=0,lw=4,c="red",ax=ax)
plt.xticks(np.linspace(0,0.5,6))
plt.xlim([0,0.5])
plt.savefig("skew.histogram.png",
			dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.show()


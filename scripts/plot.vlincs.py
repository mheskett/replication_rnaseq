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

rna_files=["/Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed"]
#"/Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.all.bed"]
dfs = []
for i in range(len(rna_files)):
	df = pd.read_csv(rna_files[i],sep="\t",
							names= ["chrom","start","stop","name","rpkm","strand","l1_density","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
	df["skew"] = df.apply(helper_func, axis = 1)
	df["sample"] = os.path.basename(rna_files[i])[0:15]
	df["informative_reads_per_kb"] = df["total_reads"] / ((df["stop"] - df["start"])  / 1000)
	add_binom_pval(df)
	dfs += [df]
df = pd.concat(dfs)
##################
df=df[df["total_reads"]>=10]
sns.kdeplot(abs(df[df["chrom"]!="X"]["skew"]),linewidth=4,cut=0)
sns.kdeplot(abs(df[df["chrom"]=="X"]["skew"]),linewidth=4,cut=0)
plt.show()

print(df[df["chrom"]=="X"])
# color_vector_rna
color_vector= []
for index,row in df.iterrows():
	if (row["fdr_reject"]==True) & (abs(row["skew"])>0.1) & (row["chrom"]!="X"):
			color_vector +=[ (0,0,1,1) ]
	if (row["fdr_reject"]==True) & (abs(row["skew"])>0.1) & (row["chrom"]=="X"):
			color_vector +=[ (1,0,0,1) ]
	if (row["fdr_reject"]==True) & (abs(row["skew"])<=0.1) & (row["chrom"]!="X"):
			color_vector +=[ (0,0,1,.05) ]
	if (row["fdr_reject"]==True) & (abs(row["skew"])<=0.1) & (row["chrom"]=="X"):
			color_vector +=[ (1,0,0,.05) ]
	if (row["fdr_reject"]==False) & (row["chrom"]!="X"):
			color_vector +=[ (0,0,1,.05) ]
	if (row["fdr_reject"]==False) & (row["chrom"]=="X"):
			color_vector +=[ (0,0,1,.05)]
# color_vector= []
# for index,row in df.iterrows():
# 	if (row["binom_pval"]<=.05) & (abs(row["skew"]>0.1)) & (row["chrom"]!="X"):
# 			color_vector +=[ (0,0,1,1) ]
# 	if (row["binom_pval"]<=.05) & (abs(row["skew"]>0.1)) & (row["chrom"]=="X"):
# 			color_vector +=[ (1,0,0,1) ]
# 	if (row["binom_pval"]<=.05) & (abs(row["skew"]<=0.1)) & (row["chrom"]!="X"):
# 			color_vector +=[ (0,0,1,.05) ]
# 	if (row["binom_pval"]<=.05) & (abs(row["skew"]<=0.1)) & (row["chrom"]=="X"):
# 			color_vector +=[ (1,0,0,.05) ]
# 	if (row["binom_pval"]>.05) & (row["chrom"]!="X"):
# 			color_vector +=[ (0,0,1,.05) ]
# 	if (row["binom_pval"]>.05) & (row["chrom"]=="X"):
# 			color_vector +=[ (0,0,1,.05)]
df["color"] = color_vector
## example plots
## filtering
f,ax=plt.subplots(figsize=(6,3))
plt.scatter(np.log2(df[df["chrom"]!="X"]["total_reads"]), abs(df[df["chrom"]!="X"]["skew"]),s=15,lw=0.05,color=df[df["chrom"]!="X"]["color"],edgecolor="black")
plt.scatter(np.log2(df[df["chrom"]=="X"]["total_reads"]), abs(df[df["chrom"]=="X"]["skew"]),s=15,lw=0.05,color=df[df["chrom"]=="X"]["color"],edgecolor="black")
plt.xlim([3,16])
plt.xticks(range(3,17))
plt.savefig(os.path.basename("gm12878.parent.rpkm.skew")+".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
###

for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	# tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	# ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_ylim([-.52,.52])
	ax.set_yticks(np.arange(-0.5,.6,.1))
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	tmp  = df[df["chrom"]==chromosomes[i]]
	for index,row in tmp[(tmp["chrom"]==chromosomes[i])].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"],fill=False,hatch="/",edgecolor=row["color"])
		ax.add_patch(rect)
	plt.savefig(os.path.basename(rna_files[0])+str(chromosomes[i])+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()


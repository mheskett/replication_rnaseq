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
## vlincs only
rna_files=["/Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.rep2.vlincs.all.bed"]
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
#######

#######
unique_genes = list(df["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
df_significant_rows = df[(df["binom_pval"]<=0.001) | (df["binom_pval"] <=0.001)]
df_nonsignificant_rows = df[(df["binom_pval"]>=0.001) | (df["binom_pval"] >=0.001)]
### switchers algorithm
for i in range(len(unique_genes)):
	samples = df_significant_rows[df_significant_rows["name"]==unique_genes[i]]
	# samples = samples[(samples["binom_pval_plus"]<=0.05) | (samples["binom_pval_minus"] <=0.05)]
	if len(samples)<=1:
		continue
	# samples = samples.reset_index(drop=True)
	# print(samples)
	hap1_skew,hap2_skew= False,False
	for index,row in samples.iterrows():
		if (row["skew"]>=0.15):
			hap1_skew = True
		if (row["skew"]<=-0.15):
			hap2_skew = True

	if hap1_skew and hap2_skew:
		switchers += [samples]
	elif hap1_skew ^ hap2_skew:
		nonswitchers += [samples]
if len(switchers)>0:
	switchers = pd.concat(switchers)
nonswitchers = pd.concat(nonswitchers)
######
## get the genes list
genes = list(switchers["name"].drop_duplicates())
genes = [x.split(",")[0] for x in genes ]
print("\n".join(genes))
###

# nonswitching gene list
genes = list(nonswitchers["name"].drop_duplicates())
genes = [x.split(",")[0] for x in genes ]
print("\n".join(genes))
print(len(genes))
#plotting shit

color_dict = {"gm12878.4.rna.5":"plum","gm12878.5.rna.5":"olivedrab","bouha13.protein":"r","bouha15.protein":"c",
"bouha10.protein":"y","bouha3.protein.":"g",
"bouha2.protein.":"b","bouha4.protein.":"m"}
color_dict_vlincs = {"bouha.trim.13Al":"r","bouha.trim.15Al":"c","bouha.trim.10Al":"y","bouha.trim.3Ali":"g",
"bouha.trim.2Ali":"b","bouha.trim.4Ali":"m","gm12878.4.rep1.":"plum","gm12878.5.rep1.":"olivedrab"}
switchers["color"]= [color_dict_vlincs[x] for x in switchers["sample"]]
nonswitchers["color"]= [color_dict_vlincs[x] for x in nonswitchers["sample"]]
legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha13',markerfacecolor='r', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha15',markerfacecolor='c', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha3',markerfacecolor='g', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha4',markerfacecolor='m', markersize=10)]

##### plot vlinc switchers
if len(switchers)>0:
	for i in range(len(chromosomes)):
		f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
		plt.suptitle(chromosomes[i])
		# tmp = switchers[switchers["chrom"]==chromosomes[i]]
		# ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
		ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
		ax.set_xlim([0, chromosome_length[chromosomes[i]]])
		ax.set_ylim([-.5,.5])
		ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))

		for index,row in switchers[(switchers["chrom"]==chromosomes[i])].iterrows():
			rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                     facecolor=row["color"], alpha=0.8,fill=False,hatch="/",edgecolor=row["color"])
			ax.add_patch(rect)
		plt.savefig("gm12878.vlinc.switchers.region."+str(chromosomes[i])+ ".png",
			dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
		plt.close()
	switchers.to_csv("gm12878.vlinc.switchers.bed",sep="\t",index=False,header=False)
###### nonswitchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	# tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	# ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_ylim([-.5,.5])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for index,row in nonswitchers[(nonswitchers["chrom"]==chromosomes[i])].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"], alpha=0.8,fill=False,hatch="/",edgecolor=row["color"])
		ax.add_patch(rect)
	plt.savefig("gm12878.vlinc.nonswitchers.region."+str(chromosomes[i])+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
nonswitchers.to_csv("gm12878.vlinc.nonswitchers.bed",sep="\t",index=False,header=False)
############## plot regions with gene labels on top and vlincs
regions=[

			# ["1",187416000,187459201],
			# ["3",174933729,174972546],
			# ["6",40074887,40115741],
			# ["6",73394774,73431597],
			# ["6",77914535,77951564],
			# ["6",120454710,120497591],
			["6",130585933,130623577], ### for paper
			["6",130265000,131160000], ### closer zoo mfor paper
			# ["6",141108314,141151646],
			# ["6",141184172,141225657],
			# ["8",2586332,2625241],
			# ["8",22568321,22607943],
			# ["9",22765327,22807785],
			# ["9",23643340,23683461],
			# ["9",24038517,24077283],
			# ["9",29802888,29843986],
			# ["15",47102705,47144857],
			# ["15",54386886,54427639],
			# ["15",54629690,54670049],
			["15",92002424,92044668], ## for paper
			["15",91900000,92450000], ## closer zoom for paper
			# ["15",96464802,96501683],
			# ["15",97096204,97134188],
]
###############

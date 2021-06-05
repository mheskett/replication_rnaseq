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

rna_files=["bouha2.protein.coding.all.counts.bed",
"bouha3.protein.coding.all.counts.bed",
"bouha4.protein.coding.all.counts.bed",
"bouha10.protein.coding.all.counts.bed",
"bouha13.protein.coding.all.counts.bed",
"bouha15.protein.coding.all.counts.bed"]
dfs = []
for i in range(len(rna_files)):
	df = pd.read_csv(rna_files[i],sep="\t",
							names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
	df["skew"] = df.apply(helper_func, axis = 1)
	df["sample"] = rna_files[i][0:15]
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
"bouha.trim.2Ali":"b","bouha.trim.4Ali":"m","gm12878.4x.hg19":"plum","gm12878.5x.hg19":"olivedrab"}
switchers["color"]= [color_dict[x] for x in switchers["sample"]]
nonswitchers["color"]= [color_dict[x] for x in nonswitchers["sample"]]
legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha13',markerfacecolor='r', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha15',markerfacecolor='c', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha3',markerfacecolor='g', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha4',markerfacecolor='m', markersize=10)]

##### plot coding switchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = switchers[switchers["chrom"]==chromosomes[i]]
	ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	plt.savefig("bouha.coding.switchers.region."+str(chromosomes[i])+ ".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
switchers.to_csv("bouha.coding.switchers.bed",sep="\t",index=False,header=False)
###### nonswitchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	plt.savefig("bouha.coding.nonswitchers.region."+str(chromosomes[i])+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
nonswitchers.to_csv("bouha.coding.nonswitchers.bed",sep="\t",index=False,header=False)
##############

###############


all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed"]
filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
	df_repli = pd.read_csv(all_files_repli[i],sep="\t",
						names= ["chrom","start","stop","hap1_early","hap1_early","hap1_late","hap2_late"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str,
						"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})

	tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
	tmp = tmp.astype(int)
	df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
	df_repli = df_repli.set_index(["chrom","start","stop"])
	df_repli = df_repli[df_repli.sum(axis="columns")!=0]
	df_repli = df_repli.reset_index()
	df_repli["sample"] = filenames_repli[i]
	repli_li.append(df_repli)
repli_df = pd.concat(repli_li)
repli_df.loc[:,"logr_hap1"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
repli_df.loc[:,"logr_hap2"] = repli_df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
repli_df.loc[:,"logr_diff"] = abs(repli_df["logr_hap1"] - repli_df["logr_hap2"]) ## 
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
sig_repli = np.percentile(a = repli_df[repli_df["chrom"]!="X"]["logr_diff"], q = 95)
repli_df["sig_repli"]=["True" if x > sig_repli else "False" for x in repli_df["logr_diff"]]
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector

# coding switchers with asynchrony....
for i in range(len(chromosomes)):
	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = switchers[switchers["chrom"]==chromosomes[i]]
	ax[0].scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j+1].set_ylim([0,1.6])
		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]

		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df[tmp_df["arm"]=="p"])
		smoothed_logrdiff_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
			return_sorted=False, frac = frac1 )

		smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
			return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_logrdiff_q,c="black",zorder=1,lw=0.8)
		ax[j+1].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j+1].add_patch(rect)
	plt.show()
### coding non switchers with asynchrony
for i in range(len(chromosomes)):
	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	ax[0].scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j+1].set_ylim([0,1.6])

		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]

		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df[tmp_df["arm"]=="p"])
		smoothed_logrdiff_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
			return_sorted=False, frac = frac1 )

		smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
			return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_logrdiff_q,c="black",zorder=1,lw=0.8)
		ax[j+1].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j+1].add_patch(rect)
	plt.show()

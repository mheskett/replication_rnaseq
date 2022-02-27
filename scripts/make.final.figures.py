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
import glob
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm

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
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],
							row["hap1_reads"]+row["hap2_reads"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)

	df["fdr_pval"]=mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")[1]
	df["fdr_reject"] =  mt.multipletests(pvals=df["binom_pval"], 
									alpha=0.01,
									method="fdr_bh")[0]
	return

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
#########################
arm_dict = get_arms(cytoband)
df = pd.read_csv("/Users/mike/replication_rnaseq/all.final.data/gm12878.4.rna.50kb.bed",sep="\t",
						names= ["chrom","start","stop","hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str})
df_repli = pd.read_csv("/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.250kb.bed",sep="\t",
						names= ["chrom","start","stop","hap1_early","hap2_early","hap1_late","hap2_late"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str,
						"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})
df_vlinc = pd.read_csv("/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.4x.hg19Aligned.outgm12878.4x.hg19Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed",sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
					"hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
					"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
##################
tmp = df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]].replace(".",0)
tmp = tmp.astype(int)
df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]] = tmp
df = df.set_index(["chrom","start","stop"])
df = df[df.sum(axis="columns")!=0]
df = df.reset_index()
###################
tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
tmp = tmp.astype(int)
df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
df_repli = df_repli.set_index(["chrom","start","stop"])
df_repli = df_repli[df_repli.sum(axis="columns")!=0]
df_repli = df_repli.reset_index()
################
## filter the zeros
df["total_plus"] = df["hap1_counts_plus"] + df["hap2_counts_plus"]
df["total_minus"] = df["hap1_counts_minus"] + df["hap2_counts_minus"]

def helper_func_plus(x):
	if x["total_plus"]==0:
		return 0
	elif x["hap1_counts_plus"] >= x["hap2_counts_plus"]:
		return x["hap1_counts_plus"]  / x["total_plus"] - 0.5
	else:
		return -x["hap2_counts_plus"]  / x["total_plus"] + 0.5
	return

def helper_func_minus(x):
	if x["total_minus"]==0: # try this for filtering
		return 0
	elif x["hap1_counts_minus"] >= x["hap2_counts_minus"]:
		return x["hap1_counts_minus"]  / x["total_minus"] - 0.5
	else:
		return -x["hap2_counts_minus"]  / x["total_minus"] + 0.5
	return


df["skew_plus"] = df.apply(helper_func_plus, axis = 1)
df["skew_minus"] = df.apply(helper_func_minus, axis = 1)

print(arm_dict.items())
df["arm"] = df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_vlinc["arm"] = df_vlinc.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_repli["arm"] = df_repli.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

print(df_repli[(df_repli["chrom"]=="1") & (df_repli["arm"]=="p")])
###########
df_repli.loc[:,"logr_hap1"] = df_repli.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
df_repli.loc[:,"logr_hap2"] = df_repli.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
df_repli.loc[:,"logr_diff"] = abs(df_repli["logr_hap1"] - df_repli["logr_hap2"])
###########

sig_repli = np.percentile(a = df_repli[df_repli["chrom"]!="X"]["logr_diff"], q = 90)
df_repli["sig_repli"]=["True" if x > sig_repli else "False" for x in df_repli["logr_diff"]]
for i in range(len(chromosomes)):
	f,ax = plt.subplots(2,1,figsize=(12,2))
	vlincs_tmp = df_vlinc[df_vlinc["chrom"]==chromosomes[i]]
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_ylim([-4,4])
	ax[0].axhline(y=0,linestyle="--",c="black",lw=0.2)
	
	if len(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]) <= 10:
		frac1 = 1
	else:
		frac1 = 10 / len(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")])
	smoothed_hap1_p = sm.nonparametric.lowess(endog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["logr_hap1"], exog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["start"], 
		return_sorted=False, frac = frac1 )
	smoothed_hap2_p = sm.nonparametric.lowess(endog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["logr_hap2"], exog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["start"],
		return_sorted=False, frac = frac1 )
	smoothed_hap1_q = sm.nonparametric.lowess(endog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["logr_hap1"], exog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"], 
		return_sorted=False, frac = 10/len(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"].index) )
	smoothed_hap2_q = sm.nonparametric.lowess(endog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["logr_hap2"], exog=df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"],
		return_sorted=False, frac = 10/len(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"].index) )
	if len(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]) >= 10:
		ax[0].plot(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["start"],smoothed_hap1_p,c="red",zorder=1,lw=0.8)
		ax[0].plot(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="p")]["start"],smoothed_hap2_p,c="blue",zorder=1,lw=0.8)
	ax[0].plot(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"],smoothed_hap1_q,c="red",zorder=1,lw=0.8)
	ax[0].plot(df_repli[(df_repli["chrom"]==chromosomes[i]) & (df_repli["arm"]=="q")]["start"],smoothed_hap2_q,c="blue",zorder=1,lw=0.8) #hap2 blue
	ax[0].set_xticks(range(0, chromosome_length[chromosomes[i]], 20000000))
	####################### repli
	color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in df_repli.iterrows() ] # red if hap1 early, blue if hap2 early
	df_repli["repli_color"] = color_vector
	for index,row in df_repli[(df_repli["sig_repli"]=="True") & (df_repli["chrom"]==chromosomes[i])].iterrows():
		rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
			facecolor=row["repli_color"],alpha=0.6,fill="True")
		ax[0].add_patch(rect)

	########################## vlincs
	ax[1].set_ylim([-0.5,0.5])
	ax[1].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[1].axhline(y=0,linestyle="--",c="black",lw=0.2)
	for index,row in vlincs_tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor="gray", alpha=0.9,fill=True) ## plot vlincs as rectangles
		ax[1].add_patch(rect)
	#### RNA
	filtered_plus = df[df["total_plus"]>=15]
	filtered_minus = df[df["total_minus"]>=15]
	ax[1].scatter(filtered_plus[filtered_plus["chrom"]==chromosomes[i]]["start"],
				filtered_plus[filtered_plus["chrom"]==chromosomes[i]]["skew_plus"],c="Blue",alpha=0.8,lw=0.2,zorder=2,edgecolor="black",s=5)
	ax[1].scatter(filtered_minus[filtered_minus["chrom"]==chromosomes[i]]["start"],
				filtered_minus[filtered_minus["chrom"]==chromosomes[i]]["skew_minus"],c="Blue",alpha=0.8,lw=0.2,zorder=2,edgecolor="black",s=5)
	ax[1].set_xticks(range(0,chromosome_length[chromosomes[i]],5000000))
	plt.show()




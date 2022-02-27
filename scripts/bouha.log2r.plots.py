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


# ## list of imprinted genes
# df_imprinted = pd.read_table("/Users/heskett/replication_rnaseq/scripts/imprinted.genes.fixed.bed",
# 	sep="\t", names=["chrom","start","stop","gene_name"],dtype={"chrom":str},
# 	header=None,index_col=None)
# print(df_imprinted)

## list of chromosome names
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

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

arm_dict = get_arms(cytoband)

df_windows = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/bouha.2e.100kb.repliseq.haplotype.counts.bed",
						sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df2_windows = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/bouha.2l.100kb.repliseq.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})


df_windows["total_reads"] = df_windows["hap1_counts"] + df_windows["hap2_counts"]
plt.hist(df_windows["hap1_counts"] / df_windows["total_reads"],bins=20)
plt.show()

df_windows["log2r_hap1"] = np.log2((df_windows["hap1_counts"]+1) / (df2_windows["hap1_counts"]+1))
df_windows["log2r_hap2"] = np.log2((df_windows["hap2_counts"]+1) / (df2_windows["hap2_counts"]+1))
df_windows["arm"] = df_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))
	# getting error for the pericentric chromosomes b/c they dont have P arms
	if len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]) <= 10:
		frac1 = 1
	else:
		frac1 = 10 / len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")])

	smoothed_hap1_p = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["log2r_hap1"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"], 
		return_sorted=False, frac = frac1 )
	smoothed_hap2_p = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["log2r_hap2"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"],
		return_sorted=False, frac = frac1 )

	smoothed_hap1_q = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["log2r_hap1"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"], 
		return_sorted = False, frac = 10/len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"].index))
	smoothed_hap2_q = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["log2r_hap2"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"],
	 return_sorted = False, frac = 10/len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"].index))
	
	if len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]) >= 10:
		ax.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"],smoothed_hap1_p,c="red",zorder=1,label="early log(hap1/hap2)",lw=1)
		ax.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"],smoothed_hap2_p,c="green",zorder=1,label="late log(hap1/hap2",lw=1)
	ax.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"],smoothed_hap1_q,c="red",zorder=1,label="early log(hap1/hap2)",lw=1)
	ax.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"],smoothed_hap2_q,c="green",zorder=1,label="late log(hap1/hap2",lw=1)

	plt.show()
	plt.close()
print(df_windows)

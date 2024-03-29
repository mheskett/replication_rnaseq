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
# ####### reads are rna seq reads, counts are repliseq reads
# df = pd.read_table("4e.combined.gm12878.4x.windows.bed",
# 	sep="\t",
# 	header=None,
# 	index_col=None,
# 	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew","hap1_counts","hap2_counts"])
# df2 = pd.read_table("4l.combined.gm12878.4x.windows.bed",
# 	sep="\t",
# 	header=None,
# 	index_col=None,
# 	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew","hap1_counts","hap2_counts"])
# #####
# df = df[(df["total_reads"]>=15) & (abs(df2["skew"]) >= 0.15)]
# df2 = df2[(df2["total_reads"]>=15 ) & (abs(df2["skew"]) >= 0.15)]
# df.loc[:,"logR"] = np.log2( (df["hap1_counts"]+1) / (df["hap2_counts"]+1) )
# df2.loc[:,"logR"] = np.log2( (df2["hap1_counts"]+1) / (df2["hap2_counts"]+1) )
# ###
# df.loc[:,"repli_skew"] = df.apply(lambda x: (x["hap1_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) - 0.5) if (x["hap1_counts"] >= x["hap2_counts"]) else 
# 															(-x["hap2_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) + 0.5), axis = 1)
# df2.loc[:,"repli_skew"] = df2.apply(lambda x: (x["hap1_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) - 0.5) if (x["hap1_counts"] >= x["hap2_counts"]) else 
# 															(-x["hap2_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) + 0.5), axis = 1)

# df = df[df["name_lncrna"].isin(df2["name_lncrna"])].reset_index()
# df2 = df2.reset_index()
# df_plotting = pd.concat( [df.loc[:,["chrom","start","stop","name_lncrna","skew"]].reset_index(drop=True), df["repli_skew"].reset_index(drop=True), df2["repli_skew"].reset_index(drop=True)],axis=1)
# df_plotting.columns = ["chrom","start","stop","name_lncrna","rna_skew","early_skew","late_skew"]
# df_plotting.loc[:,"repli_difference"] = abs(df_plotting["early_skew"] - df_plotting["late_skew"])

# f,ax = plt.subplots()
# ax.scatter(df[df["chrom"]=="X"]["start"],df[df["chrom"]=="X"]["logR"],s=30,c="blue",label = "Early Fraction log(Hap1/Hap2",lw=0.2,edgecolor="black")
# ax.scatter(df2[df2["chrom"]=="X"]["start"],df2[df2["chrom"]=="X"]["logR"],  zorder=2,s=30, c="yellow",label="Late Fraction log(Hap1/Hap2)",lw=0.2,edgecolor="black")
# ax.legend(loc="upper right")
# ax.set_ylim([-3,3])

# ax2 = ax.twinx()
# ax2.scatter(df[df["chrom"]=="X"]["start"], df[df["chrom"]=="X"]["skew"],s=30,zorder=2,c="red",label="RNA Expression Skew",lw=0.2,edgecolor="black")
# ax2.set_ylim([-.5,.5])

# for pos in df[df["chrom"]=="X"]["start"]:
# 	ax.axvline(x=pos,linestyle="--",lw=0.5,c="black")

# ax.axhline(y=0,linestyle="--",c="black",zorder=1)
# ax2.legend(loc="upper left")
# plt.show()
# plt.close()

# print(df[df["chrom"]=="X"]["skew"])
# print(df[df["chrom"]=="X"]["logR"])
# print(df2[df2["chrom"]=="X"]["logR"])

# ## plots difference versus skew
# ## color dots by hap1 vs hap2
# # shape dots by 

# f,ax = plt.subplots()
# ax.scatter(df_plotting[df_plotting["chrom"]=="X"]["rna_skew"], abs(df_plotting[df_plotting["chrom"]=="X"]["early_skew"] - df_plotting[df_plotting["chrom"]=="X"]["late_skew"]),
# 			lw=0.2,edgecolor="black")
# plt.show()
# plt.close()

# rna_skew_mean = np.mean(df_plotting[df_plotting["chrom"] != "X"]["rna_skew"])
# rna_skew_std = np.std(df_plotting[df_plotting["chrom"] != "X"]["rna_skew"])
# el_difference_mean = np.mean(abs(df_plotting[df_plotting["chrom"] != "X"]["early_skew"] - df_plotting[df_plotting["chrom"]!="X"]["late_skew"]))
# el_difference_var = np.var(abs(df_plotting[df_plotting["chrom"] != "X"]["early_skew"] - df_plotting[df_plotting["chrom"]!="X"]["late_skew"]))
# alpha = (((1-el_difference_mean)/el_difference_var) - el_difference_mean**-1) * el_difference_mean**2
# beta=alpha * ( el_difference_mean**-1 - 1)

# rna_skew_dist = scipy.stats.norm(loc=rna_skew_mean, scale=rna_skew_std)
# rna_skew_critical_values = rna_skew_dist.ppf([0.05,0.95]) # given percentiles get value
# el_difference_dist = scipy.stats.beta(a = alpha,b=beta )
# el_difference_critical_value = el_difference_dist.ppf([0.95])

# print(df_plotting[df_plotting["chrom"]!="X"].sort_values(by="repli_difference",
# 	ascending=False))
# colors = []
# for index,row in df_plotting[df_plotting["chrom"]=="15"].iterrows():
# 	# label if hap1 is early or late
# 	if (row["rna_skew"] >= rna_skew_critical_values[1]) and (row["early_skew"] >= 0) and (row["early_skew"] > row["late_skew"]) and (row["repli_difference"] >= el_difference_critical_value):
# 		colors += ["red"]
# 	elif (row["rna_skew"] >= rna_skew_critical_values[1]) and (row["late_skew"] >= 0) and (row["late_skew"] > row["early_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
# 		colors += ["green"]
# 	# label if hap2 is early or late
# 	elif (row["rna_skew"] <= rna_skew_critical_values[0]) and (row["early_skew"] <= 0) and (row["early_skew"] < row["late_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
# 		colors += ["red"] ## red for early
# 	elif (row["rna_skew"] <= rna_skew_critical_values[0]) and (row["late_skew"] <= 0) and (row["late_skew"] < row["early_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
# 		colors += ["green"] # green for late
# 	else:
# 		colors+=["black"]

# f,ax=plt.subplots()
# ax.scatter(df_plotting[df_plotting["chrom"]=="15"]["rna_skew"], abs(df_plotting[df_plotting["chrom"]=="15"]["early_skew"] - df_plotting[df_plotting["chrom"]=="15"]["late_skew"]),
# 			lw=0.2,edgecolor="black",c=colors)
# ax.axvline(x = rna_skew_critical_values[0],linestyle="--",c="black")
# ax.axvline(x = rna_skew_critical_values[1],linestyle="--",c="black")
# ax.axhline(y = el_difference_critical_value,linestyle="--",c="black")

# plt.show()
# plt.close()


####

## plot 50kb windows early vs late, and lncRNA expression

df_windows = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/bouha.2e.100kb.repliseq.haplotype.counts.bed",
						sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
print(df_windows)
df2_windows = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/bouha.2l.100kb.repliseq.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})

df_tiling_expression = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bed.tiling.all.bed",sep="\t",
						names=["chrom","start","stop","hap1_counts","hap2_counts","strand","pval","qval","reject","total_reads","skew"])
print(df_tiling_expression)
df = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.skewed.bed",sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
print(df)
df_windows = df_windows[df_windows["hap1_counts"] + df_windows["hap2_counts"] >= 15]
df2_windows = df2_windows[df2_windows["hap1_counts"] + df2_windows["hap2_counts"] >= 15]

df_windows.loc[:,"logR"] = np.log2( (df_windows["hap1_counts"]+1) / (df_windows["hap2_counts"]+1) )

df2_windows.loc[:,"logR"] = np.log2( (df2_windows["hap1_counts"]+1) / (df2_windows["hap2_counts"]+1) )

df_windows["arm"] = df_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df2_windows["arm"] = df2_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)
plt.rc('figure', titlesize=15)

# ### get asynchronous distributions
# early_thresholds = np.percentile(df_windows[df_windows["chrom"]!="X"]["logR"],[5,95])
# late_thresholds = np.percentile(df2_windows[df2_windows["chrom"]!="X"]["logR"],[5,95])


# asynch_early =  df_windows[(df_windows["chrom"] == chromosomes[i]) & ((df_windows["logR"] >= early_thresholds[1]) | (df_windows["logR"] <= early_thresholds[0]))]
# asynch_late =  df2_windows[(df2_windows["chrom"] == chromosomes[i]) & ((df2_windows["logR"] >= late_thresholds[1]) | (df2_windows["logR"] <= late_thresholds[0]))]

# print(asynch_early)

for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))
	# getting error for the pericentric chromosomes b/c they dont have P arms
	if len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]) <= 10:
		frac1 = 1
	else:
		frac1 = 10 / len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")])

	if len(df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="p")]) <= 10: 
		frac2 = 1
	else:
		frac2 = 10 / len(df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="p")])

	smoothed_early_p = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["logR"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"], 
		return_sorted=False, frac = frac1 )
	smoothed_late_p = sm.nonparametric.lowess(endog=df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="p")]["logR"], exog=df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="p")]["start"],
		return_sorted=False, frac = frac2 )

	smoothed_early_q = sm.nonparametric.lowess(endog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["logR"], exog=df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"], 
		return_sorted = False, frac = 10/len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"].index))
	smoothed_late_q = sm.nonparametric.lowess(endog=df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="q")]["logR"], exog=df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="q")]["start"],
	 return_sorted = False, frac = 10/len(df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="q")]["start"].index))

	## raw data repliseq
	# ax.scatter(df_windows[df_windows["chrom"]==chromosomes[i]]["start"], df_windows[df_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="pink",label="early hap1/hap2",alpha=0.6)
	# ax.scatter(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"], df2_windows[df2_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="green",label="late hap1/hap2",alpha=0.6)
	## smoothed repliseq
	## tiling windows rna-seq
	df_tiling_expression.loc[:,"strand_color"] = df_tiling_expression.apply(lambda x: "pink" if x["strand"]=="+" else "blue", axis = 1)
	ax.scatter(df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["start"],
			df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["skew"],
			label="all RNA Expression Skew",lw=0.2,edgecolor="black",
			c=df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["strand_color"],zorder=3,s=8,alpha=0.4)
	ax.set_ylim([-.52,.52])
	ax.set_yticks([])
	ax.set_xlim([0, chromosome_length[chromosomes[i] ] ] )
	# ax.set_xticks([])
	ax.axhline(y=0,linestyle="--",c="black")
	# smoothed repliseq middle
	ax3 = ax.twinx()
	if len(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]) >= 10:
		ax3.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="p")]["start"],smoothed_early_p,c="red",zorder=1,label="early log(hap1/hap2)",lw=0.6)
		ax3.plot(df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="p")]["start"],smoothed_late_p,c="green",zorder=1,label="late log(hap1/hap2",lw=0.6)
	ax3.plot(df_windows[(df_windows["chrom"]==chromosomes[i]) & (df_windows["arm"]=="q")]["start"],smoothed_early_q,c="red",zorder=1,label="early log(hap1/hap2)",lw=0.6)
	ax3.plot(df2_windows[(df2_windows["chrom"]==chromosomes[i]) & (df2_windows["arm"]=="q")]["start"],smoothed_late_q,c="green",zorder=1,label="late log(hap1/hap2",lw=0.6)


	significant_repliseq = df_windows[(df_windows["chrom"] == chromosomes[i]) & ((df_windows["logR"] >= 0.91) | (df_windows["logR"] <= -0.8))]

	ax3.axvline()
	if chromosomes[i]=="X":
		ax3.set_ylim([-2.5,2.5])
	else:
		ax3.set_ylim([-2,2])
	ax3.set_yticks([])
	ax3.set_xlim([0, chromosome_length[chromosomes[i] ] ] )

	ax3.axhline(y=0, linestyle="--", c="black")
	ax3.margins(x=0,y=0)

	### lncrna goes on top
	ax2 = ax.twinx()
	ax2.scatter(df[df["chrom"]==chromosomes[i]]["start"], df[df["chrom"]==chromosomes[i]]["skew"],c="orange",
		label="lncRNA Expression Skew",lw=0.2,edgecolor="black",zorder=2,s=15)
	ax2.set_ylim([-.5,.5])
	ax2.set_xlim([0, chromosome_length[chromosomes[i] ] ] )

	ax2.set_yticks([])
	# ax2.set_xticks([])

	ax.margins(x=0,y=0)
	ax2.margins(x=0,y=0)
	# ax.ticklabel_format(style='plain')
	# ax2.ticklabel_format(style='plain')
	# ax3.ticklabel_format(style='plain')
	# ax2.set_xticks([])
	### add known imprinted genes as green dots
	# ax4 = ax.twinx()
	# print(df_imprinted[df_imprinted["chrom"]==chromosomes[i]])
	# for index, row in df_imprinted[df_imprinted["chrom"]==chromosomes[i]].iterrows():
	# 	ax4.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="green", alpha=0.8)

	# plt.suptitle("chromosome " + chromosomes[i])	# plt.show()
	# plt.show()

	plt.show()
	plt.savefig("bouha.2.repli."+chromosomes[i]+".png", dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()


"""
now just make the classic log2r plots
"""



# ### plots 
# f,ax = plt.subplots()
# ax.scatter(df[df["chrom"]=="6"]["skew"], abs(df[df["chrom"]=="6"]["logR"] - df2[df2["chrom"]=="6"]["logR"]),
# 		lw=0.2,edgecolor="black" )
# plt.show()
# plt.close()
# ########
# f,ax = plt.subplots()
# ax.scatter(df[df["chrom"]=="6"]["start"],df[df["chrom"]=="6"]["logR"],s=30,c="blue",label = "Early Fraction log(Hap1/Hap2",lw=0.2,edgecolor="black")
# ax.scatter(df2[df2["chrom"]=="6"]["start"],df2[df2["chrom"]=="6"]["logR"], zorder=2,s=30,c="yellow",label="Late Fraction log(Hap1/Hap2)",lw=0.2,edgecolor="black")
# ax.legend(loc="upper right")
# ax.set_ylim([-3,3])

# ax2 = ax.twinx()
# ax2.scatter(df[df["chrom"]=="6"]["start"], df[df["chrom"]=="6"]["skew"],s=30, zorder=2,c="red",label="RNA Expression Skew",lw=0.2,edgecolor="black")
# ax2.set_ylim([-.5,.5])

# for pos in df[df["chrom"]=="6"]["start"]:
# 	ax.axvline(x=pos,linestyle="--",lw=0.5,c="black",zorder=1)



# ax.axhline(y=0,linestyle="--",c="black")
# ax2.legend(loc="upper left")
# plt.show()
# plt.close()
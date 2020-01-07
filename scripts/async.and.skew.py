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

## reads are rna seq reads, counts are repliseq reads
df = pd.read_table("4e.combined.gm12878.4x.windows.bed",
	sep="\t",
	header=None,
	index_col=None,
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew","hap1_counts","hap2_counts"])
df2 = pd.read_table("4l.combined.gm12878.4x.windows.bed",
	sep="\t",
	header=None,
	index_col=None,
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew","hap1_counts","hap2_counts"])
##
df.loc[:,"logR"] = np.log2( (df["hap1_counts"]+1) / (df["hap2_counts"]+1) )
df2.loc[:,"logR"] = np.log2( (df2["hap1_counts"]+1) / (df2["hap2_counts"]+1) )
###
df.loc[:,"repli_skew"] = df.apply(lambda x: (x["hap1_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) - 0.5) if (x["hap1_counts"] >= x["hap2_counts"]) else 
															(-x["hap2_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) + 0.5), axis = 1)
df2.loc[:,"repli_skew"] = df2.apply(lambda x: (x["hap1_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) - 0.5) if (x["hap1_counts"] >= x["hap2_counts"]) else 
															(-x["hap2_counts"]  / (x["hap1_counts"] + x["hap2_counts"]) + 0.5), axis = 1)

df = df[df["name_lncrna"].isin(df2["name_lncrna"])].reset_index()
df2 = df2.reset_index()
df_plotting = pd.concat( [df.loc[:,["chrom","start","stop","name_lncrna","skew"]].reset_index(drop=True), df["repli_skew"].reset_index(drop=True), df2["repli_skew"].reset_index(drop=True)],axis=1)
df_plotting.columns = ["chrom","start","stop","name_lncrna","rna_skew","early_skew","late_skew"]
df_plotting.loc[:,"repli_difference"] = abs(df_plotting["early_skew"] - df_plotting["late_skew"])

f,ax = plt.subplots()
ax.scatter(df[df["chrom"]=="X"]["start"],df[df["chrom"]=="X"]["logR"],s=30,c="blue",label = "Early Fraction log(Hap1/Hap2",lw=0.2,edgecolor="black")
ax.scatter(df2[df2["chrom"]=="X"]["start"],df2[df2["chrom"]=="X"]["logR"],  zorder=2,s=30,c="yellow",label="Late Fraction log(Hap1/Hap2)",lw=0.2,edgecolor="black")
ax.legend(loc="upper right")
ax.set_ylim([-3,3])

ax2 = ax.twinx()
ax2.scatter(df[df["chrom"]=="X"]["start"], df[df["chrom"]=="X"]["skew"],s=30,zorder=2,c="red",label="RNA Expression Skew",lw=0.2,edgecolor="black")
ax2.set_ylim([-.5,.5])

for pos in df[df["chrom"]=="X"]["start"]:
	ax.axvline(x=pos,linestyle="--",lw=0.5,c="black")

ax.axhline(y=0,linestyle="--",c="black",zorder=1)
ax2.legend(loc="upper left")
plt.show()
plt.close()

# print(df[df["chrom"]=="X"]["skew"])
# print(df[df["chrom"]=="X"]["logR"])
# print(df2[df2["chrom"]=="X"]["logR"])

### plots difference versus skew
### color dots by hap1 vs hap2
## shape dots by 

# f,ax = plt.subplots()
# ax.scatter(df_plotting[df_plotting["chrom"]=="X"]["rna_skew"], abs(df_plotting[df_plotting["chrom"]=="X"]["early_skew"] - df_plotting[df_plotting["chrom"]=="X"]["late_skew"]),
# 			lw=0.2,edgecolor="black")
# plt.show()
# plt.close()

rna_skew_mean = np.mean(df_plotting[df_plotting["chrom"] != "X"]["rna_skew"])
rna_skew_std = np.std(df_plotting[df_plotting["chrom"] != "X"]["rna_skew"])
el_difference_mean = np.mean(abs(df_plotting[df_plotting["chrom"] != "X"]["early_skew"] - df_plotting[df_plotting["chrom"]!="X"]["late_skew"]))
el_difference_var = np.var(abs(df_plotting[df_plotting["chrom"] != "X"]["early_skew"] - df_plotting[df_plotting["chrom"]!="X"]["late_skew"]))
alpha = (((1-el_difference_mean)/el_difference_var) - el_difference_mean**-1) * el_difference_mean**2
beta=alpha* ( el_difference_mean**-1 - 1)

rna_skew_dist = scipy.stats.norm(loc=rna_skew_mean, scale=rna_skew_std)
rna_skew_critical_values = rna_skew_dist.ppf([0.05,0.95]) # given percentiles get value
el_difference_dist = scipy.stats.beta(a = alpha,b=beta )
el_difference_critical_value = el_difference_dist.ppf([0.95])

print(df_plotting[df_plotting["chrom"]!="X"].sort_values(by="repli_difference",ascending=False))
colors = []
for index,row in df_plotting[df_plotting["chrom"]=="15"].iterrows():
	# label if hap1 is early or late
	if (row["rna_skew"] >= rna_skew_critical_values[1]) and (row["early_skew"] >= 0) and (row["early_skew"] > row["late_skew"]) and (row["repli_difference"] >= el_difference_critical_value):
		colors += ["red"]
	elif (row["rna_skew"] >= rna_skew_critical_values[1]) and (row["late_skew"] >= 0) and (row["late_skew"] > row["early_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
		colors += ["green"]
	# label if hap2 is early or late
	elif (row["rna_skew"] <= rna_skew_critical_values[0]) and (row["early_skew"] <= 0) and (row["early_skew"] < row["late_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
		colors += ["red"] ## red for early
	elif (row["rna_skew"] <= rna_skew_critical_values[0]) and (row["late_skew"] <= 0) and (row["late_skew"] < row["early_skew"]) and (row["repli_difference"]>= el_difference_critical_value):
		colors += ["green"] # green for late
	else:
		colors+=["black"]

f,ax=plt.subplots()
ax.scatter(df_plotting[df_plotting["chrom"]=="15"]["rna_skew"], abs(df_plotting[df_plotting["chrom"]=="15"]["early_skew"] - df_plotting[df_plotting["chrom"]=="15"]["late_skew"]),
			lw=0.2,edgecolor="black",c=colors)
ax.axvline(x = rna_skew_critical_values[0],linestyle="--",c="black")
ax.axvline(x = rna_skew_critical_values[1],linestyle="--",c="black")
ax.axhline(y = el_difference_critical_value,linestyle="--",c="black")

plt.show()
plt.close()


####

## plot 50kb windows early vs late, and lncRNA expression

df_windows = pd.read_csv("4e.50kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df2_windows = pd.read_csv("4l.50kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})


df_windows.loc[:,"logR"] = np.log2( (df_windows["hap1_counts"]+1) / (df_windows["hap2_counts"]+1) )

df2_windows.loc[:,"logR"] = np.log2( (df2_windows["hap1_counts"]+1) / (df2_windows["hap2_counts"]+1) )



chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]

for i in range(len(chromosomes)):
	f,ax = plt.subplots()
	smoothed_early = sm.nonparametric.lowess(endog=df_windows[df_windows["chrom"]==chromosomes[i]]["logR"], exog=df_windows[df_windows["chrom"]==chromosomes[i]]["start"], 
		return_sorted=False, frac = 20/len(df_windows[df_windows["chrom"]==chromosomes[i]]["start"].index))
	smoothed_late = sm.nonparametric.lowess(endog=df2_windows[df2_windows["chrom"]==chromosomes[i]]["logR"], exog=df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"],
	 return_sorted=False, frac = 20/len(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"].index))

	ax.scatter(df_windows[df_windows["chrom"]==chromosomes[i]]["start"], df_windows[df_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="blue",label="early hap1/hap2",alpha=0.6)
	ax.scatter(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"], df2_windows[df2_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="orange",label="late hap1/hap2",alpha=0.6)
	ax.plot(df_windows[df_windows["chrom"]==chromosomes[i]]["start"],smoothed_early)
	ax.plot(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"],smoothed_late)
	ax.set_ylim([-2,2])

	ax.axhline(y=0, linestyle="--", c="black")
	print(df[df["chrom"]==chromosomes[i]])
	ax2 = ax.twinx()
	ax2.scatter(df[df["chrom"]==chromosomes[i]]["start"], df[df["chrom"]==chromosomes[i]]["skew"],s=20,c="red",
		label="RNA Expression Skew",lw=0.2,edgecolor="black")
	ax2.set_ylim([-.5,.5])

	ax.legend()
	plt.show()
	plt.close()




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
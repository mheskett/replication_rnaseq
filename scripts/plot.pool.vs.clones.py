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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import statsmodels.stats.multitest as mt
from matplotlib_venn import venn3

# should do bedtools subtract -A

# pool = pybedtools.BedTool("gm12878.rep1.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed")
clone4 = pybedtools.BedTool("gm12878.4x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed")
clone5 = pybedtools.BedTool("gm12878.5x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed")

clone4_exclude = pybedtools.BedTool("/Users/heskett/replication_rnaseq/annotation.files/gm12878.4x.exclude.regions.bed")
clone5_exclude = pybedtools.BedTool("/Users/heskett/replication_rnaseq/annotation.files/gm12878.5x.exclude.regions.bed")


df_clone4 = clone4.subtract(clone4_exclude, A=True).to_dataframe(names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone4"])
df_clone5 = clone5.subtract(clone5_exclude, A=True).to_dataframe(names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone5"])

df_pool = pd.read_csv("gm12878.rep1.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed", 
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_pool"],
	sep="\t")

###
set1 = set(df_pool[df_pool["chrom"]!="X"]["name_lncrna"])
set2 = set(df_clone4[df_clone4["chrom"]!="X"]["name_lncrna"])
set3 = set(df_clone5[df_clone5["chrom"]!="X"]["name_lncrna"])

venn3([set1, set2, set3], ('GM12878 Parent Cell Line', 'GM12878 Single Cell Clone #4', 'GM12878 Single Cell Clone #5'))
plt.show()
plt.close()
###


# df_clone4 = pd.read_csv("gm12878.4x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed", 
# 	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone4"],
# 	sep="\t")
# df_clone5 = pd.read_csv("gm12878.5x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed",
# 	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone5"],
# 	sep="\t")

df_pool_skewed = df_pool[(abs(df_pool["skew_pool"]) >= 0.25) & (df_pool["chrom"]!="X") ]
df_clone4_skewed = df_clone4[(abs(df_clone4["skew_clone4"])>=0.25) & (df_clone4["chrom"]!="X") ]
df_clone5_skewed = df_clone5[(abs(df_clone5["skew_clone5"]) >= 0.25) & (df_clone5["chrom"]!="X") ]

print(len(df_pool_skewed.index),"skewed in pool")
print(len(df_clone4_skewed.index),"skewed in 4")
print(len(df_clone5_skewed.index),"skewed in 5")

skewed_list = list(df_pool_skewed["name_lncrna"].values) + list(df_clone4_skewed["name_lncrna"].values) \
					+ list(df_clone5_skewed["name_lncrna"].values)
# print(df_clone4[ (df_clone4["chrom"] == "11")  & (df_clone4["start"] >= 0) & (df_clone4["stop"] <= 2*10**7) ] ) 
a = df_pool[df_pool["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_pool"]].set_index("name_lncrna")
b = df_clone4[df_clone4["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_clone4"]].set_index("name_lncrna")
c = df_clone5[df_clone5["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_clone5"]].set_index("name_lncrna")

print(len(list(set(df_pool_skewed["name_lncrna"]) &  set(df_clone4_skewed["name_lncrna"]))), "skewed pool and 4")
print(len(list(set(df_pool_skewed["name_lncrna"]) &  set(df_clone5_skewed["name_lncrna"]))), "skewed pool and 5")
print(len(list(set(df_clone4_skewed["name_lncrna"]) &  set(df_clone5_skewed["name_lncrna"]))), "skewed 4 and 5")

#####
set1 = set(df_pool_skewed["name_lncrna"])
set2 = set(df_clone4_skewed["name_lncrna"])
set3 = set(df_clone5_skewed["name_lncrna"])

venn3([set1, set2, set3], ('GM12878 Parent Cell Line', 'GM12878 Single Cell Clone #4', 'GM12878 Single Cell Clone #5'))
plt.show()
plt.close()

###
all_skewed = pd.concat([a,b,c],axis=1).dropna()

print(len(all_skewed.index))
all_skewed.loc[:,"range"] = all_skewed.apply(lambda x: abs(max(x) - min(x)),axis=1)
all_skewed = all_skewed.sort_values("skew_pool",ascending=False)
## try only plotting the "changers"
all_skewed_changers = all_skewed[all_skewed["range"] >= 0.75]
all_skewed_changers = all_skewed_changers.sort_values("range",ascending=False)
print(all_skewed_changers)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
f,ax = plt.subplots(1)
ax.scatter(range(len(all_skewed_changers.index)),all_skewed_changers["skew_pool"],c="Red",s=40,label="Pool")
ax.scatter(range(len(all_skewed_changers.index)),all_skewed_changers["skew_clone4"],c="Blue",s=40,label="Clone4")
ax.scatter(range(len(all_skewed_changers.index)),all_skewed_changers["skew_clone5"],c="Green",s=40,label="Clone5")


ax.set_yticks([-1,-.75,-.5,-.25,0,.25,.5,.75,1])
# for i in range(len(all_skewed.index)):
# 	if all_skewed.at[all_skewed.index[i],"range"] >= 0.6:
# 		ax.axvline(x=i,color="red",lw=0.6)
# 	else:
# 		ax.axvline(x=i,color="black",lw=0.9)

## changers only
for i in range(len(all_skewed_changers.index)):
	if abs(all_skewed_changers.at[all_skewed_changers.index[i],"skew_pool"] ) <= 0.25:
		ax.axvline(x=i,color="red",lw=1)
	else:
		ax.axvline(x=i,color="black",lw=1)


plt.legend(prop={'size': 16})
plt.show()
plt.close()

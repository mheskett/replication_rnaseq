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

df_pool = pd.read_csv("gm12878.rep1.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed", 
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_pool"],
	sep="\t")

df_clone4 = pd.read_csv("gm12878.4x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed", 
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone4"],
	sep="\t")
df_clone5 = pd.read_csv("gm12878.5x.hg19Aligned.outgm12878.rep1.hg19Aligned.out.samtool.rmdup.1000.10000.50000.vlinc.discovery.all.bed",
	names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads","binom_pval","fdr_pval","fdr_reject","total_reads","skew_clone5"],
	sep="\t")

skewed_list = list(df_pool[df_pool["skew_pool"]>=0.25]["name_lncrna"].values) + list(df_clone4[df_clone4["skew_clone4"]>=0.25]["name_lncrna"].values) + list(df_clone5[df_clone5["skew_clone5"]>=0.25]["name_lncrna"].values)

a = df_pool[df_pool["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_pool"]].set_index("name_lncrna")
b = df_clone4[df_clone4["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_clone4"]].set_index("name_lncrna")
c = df_clone5[df_clone5["name_lncrna"].isin(skewed_list)].loc[:,["name_lncrna","skew_clone5"]].set_index("name_lncrna")
all_skewed = pd.concat([a,b,c],axis=1).dropna()
f,ax = plt.subplots(1)
ax.scatter(range(len(all_skewed.index)),all_skewed["skew_pool"],c="Red",s=20,label="Pool")
ax.scatter(range(len(all_skewed.index)),all_skewed["skew_clone4"],c="Blue",s=20,label="Clone4")
ax.scatter(range(len(all_skewed.index)),all_skewed["skew_clone5"],c="Green",s=20,label="Clone5")
for i in range(len(all_skewed.index)):
	ax.axvline(x=i,color="black",lw=0.2)
plt.legend()
plt.show()

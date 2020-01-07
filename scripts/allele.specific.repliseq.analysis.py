import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import scipy.stats

import statsmodels.api as sm

# df_late = pd.read_csv("4l.combined.overlap.na12878.hg19.haplotype.resolved.counts.bed",
# 						sep="\t",header=None,index_col=None,
# 						names=["chrom","start","stop","hap1","hap2","hap1_counts","hap2_counts"],
# 						dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"hap1_counts":int,"hap2_counts":int})
# df_early = pd.read_csv("4e.combined.overlap.na12878.hg19.haplotype.resolved.counts.bed",sep="\t",header=None,index_col=None,
# 						names=["chrom","start","stop","hap1","hap2","hap1_counts","hap2_counts"])


# df_late.loc[:,"binom_pval"] = df_late.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
# 							row["hap1_counts"]+row["hap2_counts"],
# 							p=0.5,
# 							alternative="two-sided"),axis=1)

df = pd.read_csv("4e.50kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df2 = pd.read_csv("4l.50kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})



df2.loc[:,"logR"] = np.log2( (df2["hap1_counts"]+1) / (df2["hap2_counts"]+1) )

df.loc[:,"logR"] = np.log2( (df["hap1_counts"]+1) / (df["hap2_counts"]+1) )


chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]

for i in range(len(chromosomes)):
	f,ax = plt.subplots()
	smoothed_early = sm.nonparametric.lowess(endog=df[df["chrom"]==chromosomes[i]]["logR"], exog=df[df["chrom"]==chromosomes[i]]["start"], 
		return_sorted=False, frac = 20/len(df[df["chrom"]==chromosomes[i]]["start"].index))
	smoothed_late = sm.nonparametric.lowess(endog=df2[df2["chrom"]==chromosomes[i]]["logR"], exog=df2[df2["chrom"]==chromosomes[i]]["start"],
	 return_sorted=False, frac = 20/len(df2[df2["chrom"]==chromosomes[i]]["start"].index))

	ax.scatter(df[df["chrom"]==chromosomes[i]]["start"], df[df["chrom"]==chromosomes[i]]["logR"],s=5,c="blue",label="early hap1/hap2",alpha=0.3)
	ax.scatter(df2[df2["chrom"]==chromosomes[i]]["start"], df2[df2["chrom"]==chromosomes[i]]["logR"],s=5,c="orange",label="late hap1/hap2",alpha=0.3)
	ax.plot(df[df["chrom"]==chromosomes[i]]["start"],smoothed_early)
	ax.plot(df2[df2["chrom"]==chromosomes[i]]["start"],smoothed_late)
	
	ax.axhline(y=0, linestyle="--", c="black")
	ax.legend()
	plt.show()
	plt.close()
df.loc[:,["chrom","start","logR"]].to_csv("4e.50kb.dnacopy.bed",index=None,header=True,sep="\t")
df2.loc[:,["chrom","start","logR"]].to_csv("4l.50kb.dnacopy.bed",index=None,header=True,sep="\t")

# print(df)

# df.loc[:,"percent_hap1"] = df["hap1_counts"] / (df["hap1_counts"] + df["hap2_counts"])
# df.loc[:,"binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
# 							row["hap1_counts"]+row["hap2_counts"],
# 							p=0.5,
# 							alternative="two-sided"), # v slow for some reason 
# 							axis=1)

# df.loc[:,"log_binom"] = -np.log2(df["binom_pval"])
# print(df)

# plt.scatter(df[df["chrom"]=="1"]["start"], df[df["chrom"]=="1"]["percent_hap1"])
# plt.axhline(y=0.5)
# plt.show()
# plt.close()
# plt.scatter(df[df["chrom"]=="1"]["start"], df[df["chrom"]=="1"]["log_binom"])

# plt.show()
# plt.close()

# df_late_chr6 = df_late[df_late["chrom"]=="X"]
# df_late_chr6 = df_late_chr6[df_late_chr6["hap1_counts"]+df_late_chr6["hap2_counts"]>=20]
# df_late_chr6.loc[:,"percent_hap1"] = df_late_chr6["hap1_counts"] / (df_late_chr6["hap1_counts"] + df_late_chr6["hap2_counts"])
# print(df_late_chr6)
# plt.scatter(df_late_chr6["start"],df_late_chr6["percent_hap1"],s=2)
# plt.axhline(y=0.5)
# plt.show()
# plt.plot()


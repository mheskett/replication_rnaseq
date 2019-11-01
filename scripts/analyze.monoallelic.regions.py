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

df = pd.read_csv("gm12878.rep1.hg19Aligned.out.tiling.skewed.bed",
														sep="\t",header=None,
														names=["chrom","start","stop","hap1_reads","hap2_reads",
														"strand","pval","fdrpval","reject","total_reads","skew"],
														dtype={"chrom":str,"start":int,"stop":int,"hap1_reads":int,"hap2_reads":int,
														"strand":str,"pval":float,"fdrpval":float,"reject":str,"total_reads":int,"skew":float})
df_autosomes = df[df["chrom"]!="X"]
df_autosomes_skewed = df_autosomes[(df_autosomes["pval"]<=10**-6) & df_autosomes["skew"]>=0.4]
df_X = df[df["chrom"]=="X"]
plt.hist( [abs(df_autosomes["skew"]),abs(df_X["skew"]) ], 
							bins = 30,
							label=["autosomes","X"])
#plt.hist( abs(df_X["skew"]) , bins = 30 , alpha=0.5, label="X")
sum_autosomes_skewed = (df_autosomes_skewed["stop"] - df_autosomes_skewed["start"]).sum()
plt.legend(loc="upper right" ,
		prop={'size': 15})
plt.xticks(size = 15)
plt.yticks(size = 15)
print(sum_autosomes_skewed)
plt.show()
import os
import re
import csv
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.stats.multitest as mt
from sys import argv

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



df = pd.read_csv(argv[1],sep="\t",
	names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
	dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
df=df[df["total_reads"]>=20]

add_binom_pval(df)

df.to_csv(argv[2],sep="\t",index=False,header=True)

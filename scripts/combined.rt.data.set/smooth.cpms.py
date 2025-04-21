import csv
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import statsmodels.api as sm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from scipy.stats import ttest_ind
import glob
import re
def add_binom_pval(df):
    """
    example:
    >>> scipy.stats.binomtest(5,10,p=0.5)
    BinomTestResult(k=5, n=10, alternative='two-sided', statistic=0.5, pvalue=1.0)
    """
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binomtest(row["paternal"],
                                row["paternal"]+row["maternal"],
                            p=0.5,
                            alternative="two-sided").pvalue, # v slow for some reason 
                            axis=1)
    results = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.05,
                                method="fdr_bh")
    df["fdr_pval"] = results[1]
    df["fdr_reject"] = results[0]

def smooth_vector(x,y):
    y_smooth = []
    if len(x) <= 4:
        frac = 1
    elif len(x) >4:
        frac= 4 / len(x)
    if len(x) > 0:
        y_smooth = sm.nonparametric.lowess(endog=y, 
                exog=x, 
                return_sorted=False, frac = frac )
    return y_smooth

samples = glob.glob("./acp6_113024/acp6*haplotype*counts*rmv.blck*.cpms.*s125.bed") + \
			glob.glob("./acp7_120324/acp7*haplotype*counts*rmv.blck*.cpms.*s125.bed") + \
			glob.glob("eb*rt*counts*rmv*blck*cpm*windows.s125.bed") + \
			glob.glob("gm*rt*counts*rmv.blck*cpm*windows.s125.bed")			

### this script reads a hap resolved CPM file sample line by line
### for each line, throw away any blank data marked by a "."
### make sure theres not a gap (or new chromosome)before the current window
### gather 4 (or N) data points, input to smooth_lowess
### and return the output. save the whole output to a new file

for file in samples:
	f = open(file,"r").readlines()
	dat=[]
	for line in f:
		tmp = re.split(r'[;,\t\s]\s*',line.rstrip())
		if (tmp[3]!="." and tmp[4]!="."):
			dat+=[[tmp[0], int(tmp[1]), int(tmp[2]), float(tmp[3]), float(tmp[4])]] # change data types here

	df = pd.DataFrame(dat,columns=["chrom","start","stop","paternal","maternal"])
	# print(dat)
	# print(df)
	# print(dat)
	maternal_presmoothed=[]
	paternal_presmoothed=[]
	maternal_smoothed=[]
	paternal_smoothed=[]
	for i in range(len(dat)):

		chrom = dat[i][0]
		start = dat[i][1]
		stop = dat[i][2]
		paternal = dat[i][3]
		maternal = dat[i][4]

		if i==0: # first iteration
			maternal_presmoothed += [maternal]
			paternal_presmoothed += [paternal]
			continue

		if i == len(dat)-1: # last iteration
			maternal_presmoothed += [maternal]
			paternal_presmoothed += [paternal]
			maternal_smoothed += list(smooth_vector(range(len(maternal_presmoothed)), maternal_presmoothed))
			paternal_smoothed += list(smooth_vector(range(len(paternal_presmoothed)), paternal_presmoothed))
			continue

		if (chrom == dat[i-1][0]) and ((start - int(dat[i-1][1])) == 125000): # make sure chrom is good
			maternal_presmoothed += [maternal]
			paternal_presmoothed += [paternal]

		else: # if chrom and start not good, smooth the previous and then add this
			maternal_smoothed += list(smooth_vector(range(len(maternal_presmoothed)), maternal_presmoothed))
			maternal_presmoothed = [maternal]
			paternal_smoothed += list(smooth_vector(range(len(paternal_presmoothed)), paternal_presmoothed))
			paternal_presmoothed = [paternal]

	df["paternal_smoothed"] = paternal_smoothed
	df["maternal_smoothed"] = maternal_smoothed
	df.loc[:,["chrom","start","stop","paternal_smoothed","maternal_smoothed"]]\
			.to_csv(file.removesuffix(".bed")+".smoothed.bed",sep="\t",header=None,index=None)

	if len(maternal_smoothed) != len(df):
		print(" ERROR ")
		exit()
	if len(paternal_smoothed) != len(df):
		print(" ERROR ")
		exit()	
	# print("mat len",len(maternal_smoothed))
	# print("pat len",len(paternal_smoothed))
	# print("len df",len(df))
	# exit()

	# fig,ax=	plt.subplots(nrows=2,ncols=1,figsize=(12,2))
	# plt.suptitle(file)
	# ax[0].plot(df[df["chrom"]=="10"]["start"],df[df["chrom"]=="10"]["paternal"])
	# ax[0].plot(df[df["chrom"]=="10"]["start"],df[df["chrom"]=="10"]["maternal"])
	# ax[1].plot(df[df["chrom"]=="10"]["start"],df[df["chrom"]=="10"]["paternal_smoothed"])
	# ax[1].plot(df[df["chrom"]=="10"]["start"],df[df["chrom"]=="10"]["maternal_smoothed"])
	# plt.show()

	# plt.subplots(figsize=(12,2))
	# plt.bar(x=df[df["chrom"]=="10"]["start"],height=df[df["chrom"]=="10"]["paternal_smoothed"],width=250000)
	# plt.bar(x=df[df["chrom"]=="10"]["start"],height=df[df["chrom"]=="10"]["maternal_smoothed"],width=250000)

	# plt.show()


	### testing
	# add_binom_pval(df)
	# print(df)


	# plt.bar(df[df["chrom"]=="1"]["start"],height=-np.log10(df["fdr_pval"]))



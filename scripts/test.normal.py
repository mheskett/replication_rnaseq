import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
from scipy.stats import norm

import seaborn as sns
from sys import argv
import glob
import statsmodels.api as sm
from sklearn.cluster import KMeans
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from sklearn.linear_model import LinearRegression


def add_normal_pval(df):
	# still wrong.....
	## num reads deviation from 50/50, given X total reads
	### create a normal dist based on 
    num_reads_deviation = abs(df["hap1_counts"] - df["total_reads"]/2)

    return 
def intersect_tables(df1,df2):
    ### return all df1 rows that intersect df2 by >0bs
    ### run time is n_squared ish.....dont use this for big DFs. just for small ones
    ### three possibilities of overlap: 1) start is between start/stop. 2) stop is between start/stop. 3) start <= start AND stop >= stop
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    return result

def add_binom_pval(df):
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
                            row["hap1_counts"]+row["hap2_counts"],
                            p=0.5,
                            alternative="two-sided"), # v slow for some reason 
                            axis=1)
    switcherss = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = switcherss[1]
    df["fdr_reject"] = switcherss[0]

def get_arms(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes)):
        # should be (p end, q end)
        arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
        cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
    return arm_dict
def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["hap1_counts"] >= x["hap2_counts"]:
        return x["hap1_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["hap2_counts"]  / x["total_reads"] + 0.5
    return
def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

def smooth_repli(df):
    ## returns the Y- values after smoothing
    p=[]
    if len(df[df["arm"]=="p"]) <= 10:
        frac_p = 1
    elif len(df[df["arm"]=="p"]) >10:
        frac_p= 6 / len(df[df["arm"]=="p"])

    if len(df[df["arm"]=="p"]) > 0:
        p = sm.nonparametric.lowess(endog=df[df["arm"]=="p"]["logr"], 
                exog=df[df["arm"]=="p"]["start"], 
                return_sorted=False, frac = frac_p )
    ###
    q=[]
    print(len(df[df["arm"]=="q"]) )
    if len(df[df["arm"]=="q"]) <= 10:
        frac_q = 1
    elif len(df[df["arm"]=="q"]) > 10:
        frac_q = 6 / len(df[df["arm"]=="q"])
    if len(df[df["arm"]=="p"]) > 0:
        q = sm.nonparametric.lowess(endog=df[df["arm"]=="q"]["logr"], 
            exog=df[df["arm"]=="q"]["start"], 
            return_sorted=False, frac = frac_q) 
    return p,q

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
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)

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

vlinc_files=["/Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.bed"]
dfs = []
for i in range(len(vlinc_files)):
    df = pd.read_csv(vlinc_files[i],sep="\t",
                            names= ["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
    df["skew"] = df.apply(helper_func, axis = 1)
    df["sample"] = os.path.basename(vlinc_files[i])[0:8]
    add_binom_pval(df)
    dfs += [df]
df = pd.concat(dfs)
add_normal_pval(df)
# exclude X
df=df[df["chrom"]!="X"]
########## variance by read count size groups

print(df)
group_size = []
prev_group = 0
current = 0


for i in range(1,5000):
	current = i
	sample_size = len(df[df["total_reads"].between(prev_group,current)])
	if sample_size < 100:
		continue
	if sample_size >=100:
		group_size+=[current]
		prev_group= current
print(group_size)

reads_vector = group_size
print("STD for all is : ", abs(df[df["total_reads"]>=10]["skew"]).std())
df["deviant_reads"] = abs(df["hap1_counts"] - df["total_reads"]/2)
variance_by_reads=[]
for i in range(len(reads_vector)):
	if i==0:
		variance_by_reads += [abs(df[df["total_reads"].between(0,reads_vector[i])]["deviant_reads"]).std()]
	else:
		variance_by_reads += [abs(df[df["total_reads"].between(reads_vector[i-1],reads_vector[i])]["deviant_reads"]).std()]
f,ax=plt.subplots()
ax.scatter(reads_vector,
	variance_by_reads)
plt.show()
plt.close()
##############
df["color"] = ["red" if x<=0.01 else "blue" for x in df["binom_pval"]]
x = np.array(df["total_reads"]).reshape((-1, 1))
y = np.array(abs(df["hap1_counts"] - df["total_reads"]/2))
model = LinearRegression()
model.fit(x, y)
y_pred = model.predict(np.array(range(0,5000)).reshape(-1,1))
### get significant thresholds?
min_deviation=[]
for i in range(0,5000):
	min_deviation += [model.predict(np.array([i]).reshape(1,-1))*2.5]
plt.scatter(df["total_reads"],abs(df["hap1_counts"] - df["total_reads"]/2),c=df["color"],s=15,lw=0.02,edgecolor="black")
# plt.plot(range(0,10000),[1/2*x for x in list(range(0,10000))],linestyle="--")
plt.plot(range(0,5000),y_pred,linestyle="--")
plt.plot(range(0,5000),min_deviation,linestyle="--",c="red")
plt.show()





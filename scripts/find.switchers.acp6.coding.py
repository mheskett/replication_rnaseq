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
import seaborn as sns
from scipy.stats import norm

from sys import argv
import glob
import statsmodels.api as sm
from sklearn.cluster import KMeans
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
chromosome_length = {"chr1":249250621,
"chr2":243199373,
"chr3":198022430,
"chr4":191154276,
"chr5":180915260,
"chr6":171115067,
"chr7":159138663,
"chr8":146364022,
"chr9":141213431,
"chr10":135534747,
"chr11":135006516,
"chr12":133851895,
"chr13":115169878,
"chr14":107349540,
"chr15":102531392,
"chr16":90354753,
"chr17":81195210,
"chr18":78077248,
"chr19":59128983,
"chr20":63025520,
"chr21":48129895,
"chr22":51304566,
"chrX":155270560}

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
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["paternal_counts"],
                            row["paternal_counts"]+row["maternal_counts"],
                            p=0.5,
                            alternative="two-sided"), # v slow for some reason 
                            axis=1)
    switcherss = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = switcherss[1]
    df["fdr_reject"] = switcherss[0]
    return

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
    elif x["paternal_counts"] >= x["maternal_counts"]:
        return x["paternal_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["maternal_counts"]  / x["total_reads"] + 0.5
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

arm_dict = get_arms(cytoband)

########
########
########
########
########
coding_files = glob.glob("RNA*ACP7*samtool.rmdup.genes.bed")
print(coding_files)

coding_dfs = []
for i in range(len(coding_files)):
    coding_df = pd.read_csv(coding_files[i],sep="\t",
                            names= ["chrom","start","stop","name","score","strand",
                                    "cds_start","cds_stop","paternal_counts",
                                    "maternal_counts","strand_reads","total_reads","skew","sample","informative_reads_per_kb",
                                    "binom_pval","fdr_pval","fdr_reject"],
                            dtype = {"chrom":str,"start":int,"stop":int,"paternal_counts":int,"maternal_counts":int,"name":str})
    # coding_df["total_reads"] = coding_df["paternal_counts"] + coding_df["maternal_counts"]
    # coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
    # coding_df["sample"] = coding_files[i][0:9]
    # add_binom_pval(coding_df)
    coding_dfs += [coding_df]
df_coding = pd.concat(coding_dfs)
df_coding = df_coding[df_coding["total_reads"]>=10]
df_coding = df_coding[df_coding["chrom"]!="chrX"]
# df_coding["reads_per_kb"] = df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 )
# df_coding  = df_coding[df_coding["reads_per_kb"]>=1]
# df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["paternal_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
#     .reshape(1,-1))*2.5 else False,
#     axis=1)

print(df_coding)

#####
unique_genes = list(df_coding["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
### switchers algorithm
### should i do sense only or antisense only or total reads only?
for i in range(len(unique_genes)):
    samples = df_coding[(df_coding["name"]==unique_genes[i]) &
     (~df_coding["strand_reads"].isin(["antisense","total"])) & 
     # (df_coding["strand_reads"].isin(["total"])) & 

     (df_coding["fdr_reject"]==True)]
    # samples = samples[(samples["binom_pval_plus"]<=0.05) | (samples["binom_pval_minus"] <=0.05)]
    if len(samples)<=1:
        continue
    # samples = samples.reset_index(drop=True)
    # print(samples)
    hap1_skew,hap2_skew= False,False
    for index,row in samples.iterrows():
        if (row["skew"]>=0.1):
            hap1_skew = True
        if (row["skew"]<=-0.1):
            hap2_skew = True

    if hap1_skew and hap2_skew:
        switchers += [samples]
    elif hap1_skew ^ hap2_skew:
        nonswitchers += [samples]

print(switchers)
print(len(switchers))
# print(nonswitchers)


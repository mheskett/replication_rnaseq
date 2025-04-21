import csv
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import matplotlib.patheffects as path_effects
import scipy.stats
import statsmodels.api as sm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.stats.multitest as mt
import glob

def get_switchers(df):
    # unique_genes = np.unique(df["name"])
    unique_genes = np.unique(df["name"])
    switchers = [] 
    df_significant_rows = df[(df["fdr_pval"]<=0.05) & (abs(df["aei"])>=0.2)]
    ### switchers algorithm
    for i in range(len(unique_genes)):
        samples = df_significant_rows[df_significant_rows["name"]==unique_genes[i]]
        if len(samples)<=1:
            continue

        hap1_skew, hap2_skew = False, False
        for index, row in samples.iterrows():
            if (row["aei"] >= 0.2):
                hap1_skew = True
            if (row["aei"] <= -0.2):
                hap2_skew = True

        if hap1_skew and hap2_skew:
            switchers += [samples]
    print(switchers)

    switchers = pd.concat(switchers)

    return switchers

def add_binom_pval(df):
    """
    example:
    >>> scipy.stats.binomtest(5,10,p=0.5)
    BinomTestResult(k=5, n=10, alternative='two-sided', statistic=0.5, pvalue=1.0)
    """
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binomtest(row["paternal_counts"],
                                row["paternal_counts"]+row["maternal_counts"],
                            p=0.5,
                            alternative="two-sided").pvalue, # v slow for some reason 
                            axis=1)
    results = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.05,
                                method="fdr_bh")
    df["fdr_pval"] = results[1]
    df["fdr_reject"] = results[0]


def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["paternal_counts"] >= x["maternal_counts"]:
        return x["paternal_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["maternal_counts"]  / x["total_reads"] + 0.5
    return

def sum_region_length(df):
    diffs = df["stop"] - df["start"]
    return diffs.sum()


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


chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]


arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.hg38.txt",sep="\t",
                            names =["chrom","start","stop","arm","band"])

arm_dict = get_arms(cytoband)

acp7_samples = glob.glob("./acp7*as.tl.counts.bed")

dfs=[]
dfs_with_x=[]
for x in acp7_samples:
    samp = os.path.basename(x).split(".")[0]
    df = pd.read_csv(x,sep="\t",
                            names= ["chrom","start","stop","name","reads/kb","strand","percent_l1",
                            "paternal_counts","maternal_counts","strand_reads"],
                            dtype = {"chrom":str,"strand_reads":str,"start":int,"stop":int})
    df = df[df["paternal_counts"]!="."]
    df = df[df["maternal_counts"]!="."]
    df["paternal_counts"] = df["paternal_counts"].astype(int)
    df["maternal_counts"] = df["maternal_counts"].astype(int)
    df["total_reads"] = df["paternal_counts"] + df["maternal_counts"]
    df["aei"] = df.apply(helper_func, axis = 1)
    df["sample"] = samp
    # get rid of total
    df = df[df["strand_reads"]!= "total"]
    df = df[df["strand_reads"]!= "antisense"]

    df["informative_reads_per_kb"] = df["total_reads"] / ((df["stop"] - df["start"])  / 1000)
    df = df[df["total_reads"] >= 10]
    df_with_x = df.copy()
    df = df[df["chrom"]!="chrX"]
    # filter out chrX to lower P-values on autosomes
    add_binom_pval(df)
    add_binom_pval(df_with_x)
    dfs += [df]
    dfs_with_x +=[df_with_x]

df_acp7 = pd.concat(dfs)
df_acp7_with_x = pd.concat(dfs_with_x)
df_acp7=df_acp7.sort_values(["chrom","start","sample"])
df_acp7_with_x=df_acp7_with_x.sort_values(["chrom","start","sample"])

# df_acp7["name_strand"] = df_acp7["name"]+"_"+df_acp7["strand_reads"]

df_acp7["70_percent_aei"] = abs(df_acp7["aei"]) >= 0.20
df_acp7["80_percent_aei"] = abs(df_acp7["aei"]) >= 0.30
df_acp7["90_percent_aei"] = abs(df_acp7["aei"]) >= 0.40
df_acp7["95_percent_aei"] = abs(df_acp7["aei"]) >= 0.45


# sns.kdeplot(df_acp7["informative_reads_per_kb"],clip=(0,20))
# plt.show()
# plt.close()
# sns.kdeplot(df_acp7["total_reads"],clip=(0,1000))
# plt.show()

# switchers=get_switchers(df_acp7)
# pd.Series(switchers["name"].unique()).to_csv("acp7.tls.switchers.txt",sep="\t",index=False,header=False)

df_acp7.to_csv("acp7.as.tl.counts.all.txt",sep="\t",index=None)
df_acp7.to_csv("acp7.as.tl.counts.all.bed",sep="\t",index=None,header=None)





# color_dict_acp7 = {
# "acp7_c1_rnaAligned":"red", 
# "acp7_c2_rnaAligned":"green",  
# "acp7_c4_rnaAligned":"blue"}

# #### make scatter plot dots including faded dots
# tmp = df_acp7[df_acp7["name_strand"].isin(switchers["name"])]
# tmp["color"]= [color_dict_acp7[x] for x in tmp["sample"]]
# tmp["alpha"] = tmp.apply(lambda row: 1 if (row["fdr_reject"]==True) and (abs(row["aei"])>=0.2) else 0.1, axis=1)
# tmp["unique_pos"] = [row["chrom"]+":"+row["name"] for index,row in tmp.iterrows()]
# f,ax=plt.subplots(1,1,figsize=(3,1))
# ax.scatter(tmp["unique_pos"],tmp["aei"],c=tmp["color"],s=15,edgecolor="black",lw=0.1,zorder=3,alpha=tmp["alpha"])
# for index,row in tmp.drop_duplicates(["unique_pos"]).iterrows():
#     ax.axvline(x=row["unique_pos"],linestyle="--",lw=0.4,c="black")
# plt.xticks(rotation = 280,fontsize=3)
# ax.margins(x=.015,y=0)
# ax.set_ylim([-0.53,.53])
# ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
# ax.set_yticks([-0.5,-.25,0,.25,.5])
# # f.subplots_adjust(bottom=2)
# # f.subplots_adjust(left=0.09, bottom=0.2, right=0.1, top=0.0)
# # f.subplots_adjust(right=0.7) 
# plt.savefig("acp7.rna.switchers.png",
#         dpi=300,transparent=True,bbox_inches="tight",pad_inches=0)
plt.close()
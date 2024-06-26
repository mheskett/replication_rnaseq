import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
from scipy import stats
import scipy.stats
import matplotlib.pyplot as plt
import pybedtools
import glob
import pickle
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
import pybedtools
import statistics
import math

def intersect_tables(df1,df2):
    ### return all df1 rows that intersect df2 by >0bs
    ### run time is n_squared ish.....dont use this for big DFs. just for small ones
    ### three possibilities of overlap: 1) start is between start/stop. 2) stop is between start/stop. 3) start <= start AND stop >= stop
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    #b = pybedtools.BedTool.slop(b,b=250000,g="human_g1k_v37.fasta.fai")
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    return result
def sum_bases(df):
    length = df["stop"] - df["start"]

    return length.sum()
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres                
# cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
#                             names =["chrom","start","stop","arm","band"])
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
def helper_func_rna(x):
    if x["total_reads_rna"]==0: # try this for filtering
        return 0
    elif x["hap1_counts_rna"] >= x["hap2_counts_rna"]:
        return x["hap1_counts_rna"]  / x["total_reads_rna"] - 0.5
    else:
        return -x["hap2_counts_rna"]  / x["total_reads_rna"] + 0.5
    return
rna_files=["acp6.genes.all.bed.fixed.bed"]
dfs = []
for i in range(len(rna_files)):
    df = pd.read_csv(rna_files[i],sep="\t",
                            names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
    df["skew"] = df.apply(helper_func, axis = 1)
    df["sample"] = os.path.basename(rna_files[i])[0:15]
    add_binom_pval(df)
    dfs += [df]
df = pd.concat(dfs)
df = df[df["total_reads"]>=30]
#############################
####################
print("number TLs in acp",len(df))

df["significant_deviation"] = (df["fdr_pval"]<=0.001)

# df[df["significant_deviation"]].to_csv("gm12878.tls.txt",sep="\t",header=True,index=False)
df_auto = df[df["chrom"]!="X"] 
df_auto["color"] = [(0,0,1,1) if x==True else (0,0,1,0.1) for x in df_auto["significant_deviation"]]
df["color"] = [(1,0,0,1) if x==True else (1,0,0,0.1) for x in df["significant_deviation"]]

f,ax=plt.subplots(1,1,figsize=(5,2.5))
# arrow = mpatches.FancyArrowPatch((tl187_coords[0]+1, tl187_coords[1]+0.1), (tl187_coords[0], tl187_coords[1]), mutation_scale=100)
# ax.add_patch(arrow)
ax.scatter(np.log2(df[df["chrom"]=="X"]["total_reads"]),abs(df[df["chrom"]=="X"]["skew"]),c=df[df["chrom"]=="X"]["color"],s=20,lw=0.2,edgecolor="black")
ax.scatter(np.log2(df_auto["total_reads"]),abs(df_auto["skew"]),c=df_auto["color"],s=20,lw=0.2,edgecolor="black")
ax.set_ylim([0,0.5])
# ax.set_xticks([])
plt.savefig(os.path.basename(rna_files[0])+".scatter.png",
            dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()

### add the unique TL names
df_unique_name = df.drop_duplicates(["name"])
df_unique_name["tl_name"] = ["acp:"+row["chrom"] +"-"+ str(round(row["start"]/10**6,1)) for index,row in df_unique_name.iterrows()]
df_unique_name["tl_name"] = df_unique_name.tl_name.str.cat(
    df_unique_name.groupby(['tl_name']).cumcount().add(1).astype(str),
    sep='_')
test_dict = {row["name"]:row["tl_name"] for (index,row) in df_unique_name.iterrows()}
df["tl_name"] = [test_dict[row["name"]] for index,row in df.iterrows()]

####

print(" num autosomal TLs DAE: ",len(df_auto[df_auto["significant_deviation"]==True]))
print("bases tLE with moraml Z method: ", sum_bases(df_auto[df_auto["significant_deviation"]==True]))
print("number DAE TLs per megabase ", len(df_auto[df_auto["significant_deviation"]==True])/3000)
print("number DAE TLs on X chromosome: ", len(df[(df["chrom"]=="X") & df["significant_deviation"]==True]))
print("number dae tls on 1 chromosome: ", len(df[(df["chrom"]=="1") & df["significant_deviation"]==True]))
print("number TLs on X chromosome: ",len(df[(df["chrom"]=="X")]))
print("fraction of TLs on X chromosome that are DAE: ", len(df[(df["chrom"]=="X") & df["significant_deviation"]==True]) / (len(df[(df["chrom"]=="X")])))
# df[["chrom","start","stop","tl_name","rpkm",
#         "strand","rpkm","l1_density",
#         "hap1_counts","hap2_counts","skew","sample","binom_pval",
#         "fdr_pval","fdr_reject","significant_deviation"]].to_csv("gm12878.tls.txt",sep="\t",header=True,index=False)

# df_auto[df_auto["significant_deviation"]==True]["chrom"].value_counts().plot(kind="bar")
# plt.show()
##################

#######
for i in range(len(chromosomes)):
    f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
    plt.suptitle(chromosomes[i])
    # tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
    # ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
    ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax.set_xlim([0, chromosome_length[chromosomes[i]]])
    ax.set_ylim([-.52,.52])
    ax.set_yticks(np.arange(-0.5,.6,.1))
    ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
    tmp  = df[df["chrom"]==chromosomes[i]]
    for index,row in tmp[(tmp["chrom"]==chromosomes[i])].iterrows():
        rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"],fill=False,hatch="/",edgecolor=row["color"])
        ax.add_patch(rect)
    plt.savefig(os.path.basename(rna_files[0])+"vlinc.only."+str(chromosomes[i])+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

plt.close()

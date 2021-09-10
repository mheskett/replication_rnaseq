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
import scipy.stats
import pickle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
import pybedtools
import statistics

def closest_gene(df_lnc,df_coding):
    a = pybedtools.BedTool.from_dataframe(df_lnc).sort()
    b = pybedtools.BedTool.from_dataframe(df_coding).sort()
    result = a.closest(b,d=True).to_dataframe(names=list(df_lnc.columns) + [x+'1' for x in df_coding.columns]+["dist"],
        dtype={"chrom":str})
    return result

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
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
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



###
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))
arm_dict = get_arms(cytoband)

rna_files=["/Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed"]
dfs = []
for i in range(len(rna_files)):
    df = pd.read_csv(rna_files[i],sep="\t",
                            names= ["chrom","start","stop","name","rpkm","strand","l1_density","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
    df["skew"] = df.apply(helper_func, axis = 1)
    df["sample"] = os.path.basename(rna_files[i])[0:15]
    add_binom_pval(df)
    dfs += [df]
df = pd.concat(dfs)
df = df[df["total_reads"]>=10]
df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
    axis=1)
df_auto = df[df["chrom"]!="X"] 

#####################################'
## do a quick plot of coding & lncrna
####
coding_files=["gm12878.rep1.protein.coding.all.counts.bed"]
coding_dfs = []
for i in range(len(coding_files)):
    coding_df = pd.read_csv(coding_files[i],sep="\t",
                            names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    coding_df["total_reads"] = coding_df["hap1_counts"] + coding_df["hap2_counts"]
    coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
    coding_df["sample"] = coding_files[i][0:7]
    add_binom_pval(coding_df)
    coding_dfs += [coding_df]
df_coding = pd.concat(coding_dfs)

#######
df_coding = df_coding[df_coding["total_reads"]>=10]
df_coding = df_coding[df_coding["chrom"]!="X"]
#### figure that shows coding gene and lncrnas together
####
df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
    .reshape(1,-1))*2.5 else False,
    axis=1)

#####
# plot all lncs and genes
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
    tmp_lnc  = df_auto[(df_auto["chrom"]==chromosomes[i]) & (df_auto["significant_deviation"]==True)]
    
    for index,row in tmp_lnc.iterrows():
        rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor="mediumblue",fill=True,edgecolor="black")
        ax.add_patch(rect)
    tmp_coding = df_coding[(df_coding["chrom"]==chromosomes[i]) & (df_coding["significant_deviation"]==True)]
    for index,row in tmp_coding.iterrows():
        rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor="red",fill=True,hatch="/",edgecolor="red")
        ax.add_patch(rect)
    # plt.savefig(os.path.basename(rna_files[0])+"vlinc.only."+str(chromosomes[i])+".png",
    #     dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    # plt.show()
    plt.close()
## what about this: get closest gene and check cis-trans relationships with lncRNAs
dae_lnc_coding = closest_gene(df_auto[df_auto["significant_deviation"]==True],df_coding[df_coding["significant_deviation"]==True])
dae_lnc_coding["cis_trans"] = dae_lnc_coding.apply(lambda x: "cis" if np.sign(x["skew"])==np.sign(x["skew1"]) else "trans",axis=1)
cis_genes = dae_lnc_coding[dae_lnc_coding["cis_trans"]=="cis"].sort_values(["dist"])
cis_genes_head = cis_genes.head(30)
print(cis_genes_head)
for index,row in cis_genes_head.iterrows():
    f,ax = plt.subplots(1,1,figsize=(4,2),sharex=False)
    chrom=row["chrom"]
    start=row["start"]
    stop=row["stop"]
    print(chrom,start,stop)
    # plt.suptitle()
    # tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
    # ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
    ax.axhline(y=0, linestyle="--", lw=0.4, c="black")
    ax.set_xlim([start-100000, stop+100000])
    ax.set_ylim([-.52,.52])
    ax.set_yticks(np.arange(-0.5,.6,.1))
    ax.set_xticks(np.linspace(start-100000, stop+100000, 4))
    tmp_lnc = df_auto[(df_auto["chrom"]==chrom) & (df_auto["start"].between(start-100000,start+100000)) ]
    for index2,row2 in tmp_lnc.iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.025), width=row2["stop"]-row2["start"], height=0.05,
                     facecolor="red" if row2["strand"]=="+" else "blue",fill=True,edgecolor="black",lw=0.5)
        ax.add_patch(rect)
        ax.text(row2["start"]+10000,row2["skew"]-0.1,str("TL:")+str(chrom),size=10) ## to get the gene names

    tmp_coding = df_coding[(df_coding["chrom"]==chrom) & (df_coding["start"].between(start-1000000,start+1000000)) ]
    for index3,row3 in tmp_coding.iterrows():
        rect=Rectangle((row3["start"], row3["skew"]-.025), width=row3["stop"]-row3["start"], height=0.05,
                     facecolor="red" if row2["strand"]=="+" else "blue",fill=True,edgecolor="black",lw=0.5)
        ax.add_patch(rect)
        ax.text(row3["start"]+10000,row3["skew"]-0.1,row3["name"][0:15],size=10) ## to get the gene names
    # plt.savefig(os.path.basename(rna_files[0])+"vlinc.only."+str(chromosomes[i])+".png",
    #     dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.show()
    plt.close()
        # x = range(start,stop)
        #     y = range(row["start"],row["stop"])
        #     # if range(max(x[0], y[0]), min(x[-1], y[-1])+1):
        #     #   ax.text(row["start"]+10000,row["skew"]-0.1,row["name"][0:15],size=5) ## to get the gene names


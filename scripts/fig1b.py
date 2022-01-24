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
import glob
import pickle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
import pybedtools
import statistics

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
#############################
####################
print("number TLs in gm12878",len(df))
df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
    axis=1)
df["significant_deviation"] = (df["significant_deviation"]==True) & (df["fdr_pval"]<=0.01)
df_auto = df[df["chrom"]!="X"] 
df_auto["color"] = [(0,0,1,1) if x==True else (0,0,1,0.1) for x in df_auto["significant_deviation"]]
df["color"] = [(1,0,0,1) if x==True else (1,0,0,0.1) for x in df["significant_deviation"]]

f,ax=plt.subplots(1,1,figsize=(5,2.5))
ax.scatter(np.log2(df[df["chrom"]=="X"]["total_reads"]),abs(df[df["chrom"]=="X"]["skew"]),c=df[df["chrom"]=="X"]["color"],s=20,lw=0.2,edgecolor="black")
ax.scatter(np.log2(df_auto["total_reads"]),abs(df_auto["skew"]),c=df_auto["color"],s=20,lw=0.2,edgecolor="black")
ax.set_ylim([0,0.5])
# ax.set_xticks([])
plt.savefig("fig1.scatter.png",
            dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
print(" num autosomal TLs with normal Z method: ",len(df_auto[df_auto["significant_deviation"]==True]))
print("bases tLE with moraml Z method: ", sum_bases(df_auto[df_auto["significant_deviation"]==True]))
print("number DAE TLs per megabase ", len(df_auto[df_auto["significant_deviation"]==True])/3000)
print("number DAE TLs on X chromosome: ", len(df[(df["chrom"]=="X") & df["significant_deviation"]==True]))
print("number dae tls on 1 chromosome: ", len(df[(df["chrom"]=="1") & df["significant_deviation"]==True]))

# df_auto[df_auto["significant_deviation"]==True]["chrom"].value_counts().plot(kind="bar")
# plt.show()
##################
exit()
df_auto.to_csv("gm12878.rep1.vlincs.dae.list.bed",sep="\t")

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
    tmp  = df_auto[df_auto["chrom"]==chromosomes[i]]
    for index,row in tmp[(tmp["chrom"]==chromosomes[i])].iterrows():
        rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"],fill=False,hatch="/",edgecolor=row["color"])
        ax.add_patch(rect)
    plt.savefig(os.path.basename(rna_files[0])+"vlinc.only."+str(chromosomes[i])+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

plt.close()
#### do intergenic analysis plot
df_genes = pd.read_csv("/Users/mike/replication_rnaseq/scripts/ucsc.genes.cds.only.filtered.bed",sep="\t",header=None,index_col=None,
    names=["chrom","start","stop","name","score","strand"])
df_genes = df_genes.drop_duplicates(["chrom","start","strand"])

lnc_intersect_coding = intersect_tables(df_auto[df_auto["significant_deviation"]==True],df_genes)

lnc_start_to_coding_anything = [min(abs(row["start"]-row["start1"]),abs(row["start"]-row["stop1"])) 
                                    if row["strand"]=="+" else min(abs(row["stop"]-row["start1"]),abs(row["stop"]-row["stop1"])) for index,row in lnc_intersect_coding.iterrows()]

def closest_gene(df_lnc,df_coding):
    a = pybedtools.BedTool.from_dataframe(df_lnc).sort()
    b = pybedtools.BedTool.from_dataframe(df_coding).sort()
    result = a.closest(b,d=True).to_dataframe(names=list(df_lnc.columns) + [x+'1' for x in df_coding.columns]+["dist"])
    return result

dae_lnc_closest_gene = closest_gene(df_auto[df_auto["significant_deviation"]==True],df_genes)
print("average distance of DAE TLs to coding gene: ",statistics.mean(dae_lnc_closest_gene["dist"]))
print("median distance of DAE TLs to coding gene: ",statistics.median(dae_lnc_closest_gene["dist"]))

# print(dae_lnc_closest_gene)
# plt.figure()
# sns.kdeplot(dae_lnc_closest_gene["dist"],cut=0)

# plt.show()

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
                     facecolor="blue",fill=False,hatch="/",edgecolor="blue")
        ax.add_patch(rect)
    tmp_coding = df_coding[(df_coding["chrom"]==chromosomes[i]) & (df_coding["significant_deviation"]==True)]
    for index,row in tmp_coding.iterrows():
        rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor="red",fill=False,hatch="/",edgecolor="red")
        ax.add_patch(rect)
    # plt.savefig(os.path.basename(rna_files[0])+"vlinc.only."+str(chromosomes[i])+".png",
    #     dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.show()
    plt.close()
## what about this: get closest gene and check cis-trans relationships with lncRNAs
dae_lnc_coding = closest_gene(df_auto[df_auto["significant_deviation"]==True],df_coding[df_coding["significant_deviation"]==True])
dae_lnc_coding["cis_trans"] = dae_lnc_coding.apply(lambda x: "cis" if np.sign(x["skew"])==np.sign(x["skew1"]) else "trans",axis=1)
print(dae_lnc_coding.groupby("cis_trans").count())
print("average distance for lncRNA that have cis dae coding genes nearest to them ",dae_lnc_coding[dae_lnc_coding["cis_trans"]=="cis"]["dist"].mean())
print("average distance for lncRNA that have trans dae coding genes nearest to them ",dae_lnc_coding[dae_lnc_coding["cis_trans"]=="trans"]["dist"].mean())

print("median distance for lncRNA that have cis dae coding genes nearest to them ",dae_lnc_coding[dae_lnc_coding["cis_trans"]=="cis"]["dist"].median())
print("median distance for lncRNA that have trans dae coding genes nearest to them ",dae_lnc_coding[dae_lnc_coding["cis_trans"]=="trans"]["dist"].median())

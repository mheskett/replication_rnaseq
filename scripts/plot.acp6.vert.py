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

arm_dict = get_arms(cytoband)
all_files_repli = [ "acp6_c1_m_el.log2rt.as.only.bed", 
                    "acp6_c2_p_el.log2rt.as.only.bed", 
                    "acp6_c6_m_el.log2rt.as.only.bed",
                    "acp6_c1_p_el.log2rt.as.only.bed", 
                    "acp6_c5_m_el.log2rt.as.only.bed", 
                    "acp6_c6_p_el.log2rt.as.only.bed", 
                    "acp6_c2_m_el.log2rt.as.only.bed", 
                    "acp6_c5_p_el.log2rt.as.only.bed"]

filenames_repli=[os.path.basename(x)[0:12] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",
                        names= ["chrom","start","stop","log2rt"],
                        dtype = {"chrom":str,"start":int,"stop":int,"log2rt":float})

    # tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
    # tmp = tmp.astype(int)
    # df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
    df_repli = df_repli.set_index(["chrom","start","stop"])
    # df_repli = df_repli[df_repli.sum(axis="columns")!=0]
    df_repli = df_repli.reset_index()
    df_repli["sample"] = filenames_repli[i]
    repli_li.append(df_repli)

repli_df = pd.concat(repli_li)
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
repli_df["haplotype"] = ["paternal" if "_p_" in x else "maternal" for x in repli_df["sample"]]
zscore = lambda x: (x - x.mean()) / x.std()

# print(repli_df)

### investigate dropna vs non dropna data. its probably like 10-15% of rows
## be aware that youre droping NAs in ANY row here...maybe shouldnt?
df_normalized_paternal = quantile_normalize(repli_df[repli_df["haplotype"]=="paternal"].pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any")).reset_index()
df_normalized_maternal = quantile_normalize(repli_df[repli_df["haplotype"]=="maternal"].pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any")).reset_index()
df_normalized_paternal["std_dev"] = df_normalized_paternal.filter(like="acp6",axis=1).std(axis="columns")
df_normalized_maternal["std_dev"] = df_normalized_maternal.filter(like="acp6",axis=1).std(axis="columns")
df_normalized_paternal["arm"] = df_normalized_paternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_normalized_maternal["arm"] = df_normalized_maternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)


# print(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any"))

df_normalized_both_haps = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any")).reset_index()
df_normalized_both_haps["std_dev"] = df_normalized_both_haps.filter(like="acp6",axis=1).std(axis="columns")
df_normalized_both_haps["arm"] = df_normalized_both_haps.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

df_normalized_both_haps["std_dev_zscore"] = df_normalized_both_haps["std_dev"].transform(zscore)
df_normalized_both_haps["p_values"] = scipy.stats.norm.sf(abs(df_normalized_both_haps["std_dev_zscore"])) #one-sided
# df_normalized_both_haps.index.name=None
df_normalized_both_haps = df_normalized_both_haps.rename_axis(None, axis=1)

# print(df_normalized_both_haps)
# print(df_normalized_paternal)
# print(df_normalized_maternal)
# print(df_normalized_both_haps.sort_values("p_values",ascending=True))


## change appropriately. color is clonal sample, then dotted vs colored is maternal paternal
df_normalized_both_haps.to_csv("acp6.rt.vert.txt",sep="\t",header=True,index=False)
df_normalized_both_haps.to_csv("acp6.rt.vert.bed",sep="\t",header=False,index=False)
######
#####
######

### get the top 30 P-value VERT regions and do plot zooms on them. many will be
### the same regions just different windows. 
samples=repli_df["sample"].unique()
color_dict_repli = {"acp6_c1_m_el":"red", 
"acp6_c1_p_el":"red",  
"acp6_c2_m_el":"cyan",  
"acp6_c2_p_el":"cyan",
"acp6_c5_m_el":"yellow",  
"acp6_c5_p_el":"yellow",  
"acp6_c6_m_el":"green",  
"acp6_c6_p_el":"green"}
# get top 30
# then do a bedtools merge, then reconvert back to DF, reference those positions but use
## the original dataframe for plotting the locations?
### zoomed into regions plots
tops = df_normalized_both_haps[df_normalized_both_haps["std_dev_zscore"]>=2.5]
print(pybedtools.BedTool.from_dataframe(tops.sort_values(["chrom","start"])\
    .loc[:,["chrom","start","stop"]]).merge(d=600000).to_dataframe())
for index, row in tops.iterrows():
    f, ax = plt.subplots(1,1,figsize=(3,4))
    for sample in samples:

        repli_sample = df_normalized_both_haps.loc[:,["chrom","start","stop",sample,"arm","std_dev_zscore"]]
        repli_sample = repli_sample.sort_values("start",ascending=True)
        repli_sample = repli_sample[(repli_sample["chrom"]==row["chrom"]) 
                                    & (repli_sample["start"]>=row["start"]-3000000)
                                    & (repli_sample["stop"]<=row["stop"]+3000000)]
        print( " repli sample ")
        print(repli_sample)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10) 

        ax.axhline(y=0,lw=0.5,c="black")
        ### maybe the data is not sorted before being smoothed??? ...fixed...
        ### try without smoothing?
        ax.plot(repli_sample[repli_sample["arm"]=="p"]["start"],
                    #smooth_vector(repli_sample[repli_sample["arm"]=="p"]["start"], repli_sample[repli_sample["arm"]=="p"][sample]), 
                    repli_sample[repli_sample["arm"]=="p"][sample],
                    c=color_dict_repli[sample],lw=0.8,linestyle="--" if "_p_" in sample else "-")

        ax.plot(repli_sample[repli_sample["arm"]=="q"]["start"],
                #smooth_vector(repli_sample[repli_sample["arm"]=="q"]["start"],repli_sample[repli_sample["arm"]=="q"][sample]), 
                repli_sample[repli_sample["arm"]=="q"][sample],
                c=color_dict_repli[sample],lw=0.8,linestyle="--" if "_p_" in sample else "-" ) ## -- line style is haplotype 2
    
    for index,row in repli_sample[(repli_sample["std_dev_zscore"]>=2.5)].iterrows():
        rect=Rectangle((row["start"],-5),width=row["stop"]-row["start"],height=10,
            facecolor="gray",alpha=0.6,fill="True")
        ax.add_patch(rect)      
    # for index,row in df_chrom_q[(df_chrom_q["logr_diff_abs_sample_zscore"]>=2.5)].iterrows():
    #     rect=Rectangle((row["start"],min_val),width=row["stop"]-row["start"],height=max_val+abs(min_val),
    #         facecolor=row["repli_color"],alpha=0.6,fill="True")


    ax.set_ylim([-4,4])
    ax.set_yticks([-4,-3,-2,-1,0,1,2,3,4])
    ax.set_xlim(row["start"]-2500000,row["stop"]+2500000)
    ax.set_xticks(np.linspace(row["start"]-2500000,row["stop"]+2500000,5))
    plt.suptitle("acp6.rt.vert.zoom"+row["chrom"]+"-"+str(row["start"]))
    plt.savefig("acp6.rt.vert.zoom."+row["chrom"]+"-"+str(row["start"])+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()
######
######
######
#####


####
####
# whole chromosome plots
#####
# #######
for chrom in chromosomes:
    repli_df_chrom = df_normalized_both_haps[(df_normalized_both_haps["chrom"]==chrom)]
    f, ax = plt.subplots(1,1,figsize=(30,3))
    # print("df chrom no sample selection",repli_df_chrom)
    for sample in samples:
        # print("sample:",sample)
        # print("df working with: ",repli_df_chrom)
        repli_sample = repli_df_chrom.loc[:,["chrom","start","stop",sample,"arm","std_dev_zscore"]]
        repli_sample = repli_sample.sort_values("start",ascending=True)
        # repli_sample.to_csv("test.txt",sep="\t")
        # print(repli_sample)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10) 

        # repli_df_chrom = repli_sample[(repli_sample["chrom"]==chrom)]
        ax.axhline(y=0,lw=0.5,c="black")
        ### maybe the data is not sorted before being smoothed???
        ax.plot(repli_sample[repli_sample["arm"]=="p"]["start"],
                    smooth_vector(repli_sample[repli_sample["arm"]=="p"]["start"], 
                        repli_sample[repli_sample["arm"]=="p"][sample]), c=color_dict_repli[sample],lw=0.8,linestyle="--" if "_p_" in sample else "-")
        # print("x coords: ", repli_sample[repli_sample["arm"]=="p"]["start"])
        # print("smoothed: ",smooth_vector(repli_sample[repli_sample["arm"]=="p"]["start"],repli_sample[repli_sample["arm"]=="p"][sample]))
        ax.plot(repli_sample[repli_sample["arm"]=="q"]["start"],
                smooth_vector(repli_sample[repli_sample["arm"]=="q"]["start"],
                    repli_sample[repli_sample["arm"]=="q"][sample]), c=color_dict_repli[sample],lw=0.8,linestyle="--" if "_p_" in sample else "-" ) ## -- line style is haplotype 2
    

    for index,row in repli_df_chrom[(repli_df_chrom["std_dev_zscore"]>=2.5)].iterrows():
        rect=Rectangle((row["start"],-5),width=row["stop"]-row["start"],height=10,
            facecolor="gray",alpha=0.6,fill="True")
        ax.add_patch(rect)      
    # for index,row in df_chrom_q[(df_chrom_q["logr_diff_abs_sample_zscore"]>=2.5)].iterrows():
    #     rect=Rectangle((row["start"],min_val),width=row["stop"]-row["start"],height=max_val+abs(min_val),
    #         facecolor=row["repli_color"],alpha=0.6,fill="True")


    ax.set_ylim([-4,4])
    ax.set_yticks([-4,-3,-2,-1,0,1,2,3,4])
    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length[chrom],20))
    plt.suptitle("acp6.rt.vert"+chrom)

    plt.savefig("acp6.rt.vert."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()
##############





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
############## get the repliseq with logr diff for each sample
all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.500kb.bed"]
filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",
                        names= ["chrom","start","stop","hap1_early","hap2_early","hap1_late","hap2_late"],
                        dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str,
                        "hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})

    tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
    tmp = tmp.astype(int)
    df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
    df_repli = df_repli.set_index(["chrom","start","stop"])
    df_repli = df_repli[df_repli.sum(axis="columns")!=0]
    df_repli = df_repli.reset_index()
    df_repli["sample"] = filenames_repli[i]
    repli_li.append(df_repli)
repli_df = pd.concat(repli_li)
repli_df.loc[:,"logr_hap1"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
repli_df.loc[:,"logr_hap2"] = repli_df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
repli_df.loc[:,"logr_diff_abs"] = abs(repli_df["logr_hap1"] - repli_df["logr_hap2"]) ## 
repli_df.loc[:,"logr_diff_raw"] = repli_df["logr_hap1"] - repli_df["logr_hap2"] # positive if hap1 early, negative if hap2 early
repli_df.loc[:,"logr"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+x["hap2_early"]+1) / (x["hap1_late"]+x["hap2_late"]+1)), axis=1 )
sig_repli = np.percentile(a = repli_df[repli_df["chrom"]!="X"]["logr_diff_abs"], q = 95)
repli_df["sig_repli"]=["True" if x > sig_repli else "False" for x in repli_df["logr_diff_raw"]]
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector
asynchronous_regions = repli_df[repli_df["sig_repli"]=="True"]
#####
zscore = lambda x: (x - x.mean()) / x.std()

# for i in repli_df["samples"].unique():
#     mean=repli_df[(repli_df["chrom"]!="X")&(repli_df["sample"]==i)]["logr_diff_abs"].mean()
#     std=repli_df[(repli_df["chrom"]!="X")&(repli_df["sample"]==i)]["logr_diff_abs"].std()


repli_df = repli_df[repli_df["chrom"]!="X"]
repli_df["logr_diff_abs_sample_zscore"] = repli_df.groupby("sample")["logr_diff_abs"].transform(zscore)
print(repli_df)
###################################

### now get all the TLs an DAE TLs in EB2 and EB10
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
unique_genes = list(df["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
df_significant_rows = df[df["binom_pval"]<=0.001]
df_nonsignificant_rows = df[df["binom_pval"] >=0.001]
### switchers algorithm
for i in range(len(unique_genes)):
    samples = df_significant_rows[df_significant_rows["name"]==unique_genes[i]]
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
switchers = pd.concat(switchers)
color_dict = {"bouha.4.":"r","bouha.15":"c","bouha.10":"orange","bouha.3.":"g",
"bouha.2.":"b","bouha.13":"green"}
switchers["color"]= [color_dict[x] for x in switchers["sample"]]
df["color"]=[color_dict[x] for x in df["sample"]]


##########
print("standard deviation of logr difference")



#########
### organize this code because its actually pretty decent!
tmp1 = intersect_tables(df[(df["fdr_reject"]==True) & (df["chrom"]!="X") & (abs(df["skew"])>=0.1) & (df["sample"].isin(["bouha.2.","bouha.10"]))],asynchronous_regions)
f,ax = plt.subplots(figsize=(4,4))
plt.scatter(abs(tmp1["skew"]),tmp1["logr_diff_abs1"],s=30,edgecolor="black",lw=0.2)
# plt.show()
plt.close()
##########
repli_diff = repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_diff_raw").reset_index()
repli_diff["bouha2.bouha10"] = repli_diff["bouha.2.repli.5"] - repli_diff["bouha.10.repli."]
tmp = intersect_tables(df[(df["fdr_reject"]==True) & (df["chrom"]!="X") & (abs(df["skew"])>=0.1) & (df["sample"].isin(["bouha.2.","bouha.10"]))],
    repli_diff)
#######
f,ax = plt.subplots(figsize=(4,4))
plt.scatter(abs(tmp["skew"]),abs(tmp["bouha2.bouha101"]),s=30,edgecolor="black",lw=0.2)
# plt.show()
plt.close()
head = tmp.sort_values("bouha2.bouha101").head(20).append(tmp.sort_values("bouha2.bouha101").tail(20))

#################
most_asrt = tmp1.sort_values(["logr_diff_abs1"]).head(20)
# Using the variation in log difference To make error bars
standard_deviation_2 = repli_df[(repli_df["sample"]=="bouha.2.repli.5")]["logr_diff_abs"].std()
standard_deviation_10 = repli_df[(repli_df["sample"]=="bouha.10.repli.")]["logr_diff_abs"].std()
for index,row in most_asrt.drop_duplicates(["chrom","start","stop"]).iterrows():
    f,ax = plt.subplots(1,1,figsize=(4,2),sharex=False)
    start=row["start"]
    stop=row["stop"]
    chrom=str(row["chrom"])
    plt.suptitle(chrom)
    for index2,row2 in most_asrt[(most_asrt["chrom"]==chrom) & (most_asrt["start"]==start) & (most_asrt["stop"]==stop) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.05), width=row2["stop"]-row2["start"], height=0.1,
                     facecolor=row2["color"], edgecolor=row2["color"],hatch="/",fill=False) ## plot vlincs as rectangles
            # rectvar = Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.16,
         #                 facecolor="gray", edgecolor="gray",hatch="/",fill=False) ## plot vlincs as rectangles
        ax.add_patch(rect)
    ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax.set_xlim([max(0,start-2000000),stop+2000000])
    ax.set_ylim([-0.6,0.6])
    ax.set_xticks(np.linspace(max(0,start-2000000),stop+2000000, 8))

    ax2=ax.twinx()
    # (repli_df["sample"]=="bouha.2.repli.5") & (repli_df["chrom"]==chrom) 
    bouha2 = repli_df[(repli_df["sample"]=="bouha.2.repli.5") &(repli_df["chrom"]==chrom)]
    bouha10 = repli_df[(repli_df["sample"]=="bouha.10.repli.")& (repli_df["chrom"]==chrom)]# & (repli_df["start"]>=start-3000000) & (repli_df["stop"]<=stop+3000000)]
    # for index,row in bouha2.iterrows():
    # ax2.plot(bouha2["start"],smooth_vector(bouha2["start"],bouha2["logr_diff_abs"]))#smooth_vector(list(bouha2["start"]),list(bouha2["logr_diff_abs"])))
    # ax2.plot(bouha10["start"],smooth_vector(bouha10["start"],bouha10["logr_diff_abs"]))
    ## unsmoothed
    ax2.plot(bouha2["start"],bouha2["logr_diff_abs"])#smooth_vector(list(bouha2["start"]),list(bouha2["logr_diff_abs"])))
    ax2.plot(bouha10["start"],bouha10["logr_diff_abs"])
    #plot thick lines for STD ?
    ax2.fill_between(bouha2["start"],bouha2["logr_diff_abs"]- standard_deviation_2/2, bouha2["logr_diff_abs"] + standard_deviation_2/2,
                 color='blue', alpha=0.1)
    ax2.fill_between(bouha10["start"],bouha10["logr_diff_abs"]- standard_deviation_10/2, bouha10["logr_diff_abs"] + standard_deviation_10/2,
                 color='orange', alpha=0.1)
    ##########0 
    ax2.set_ylim([0,max(max(bouha2[(bouha2["start"]>=start-3000000) & (bouha2["start"]<=stop+3000000)]["logr_diff_abs"]),
        max(bouha10[(bouha10["start"]>=start-3000000) & (bouha10["start"]<=stop+3000000)]["logr_diff_abs"]) )+0.1])
       
    plt.savefig("vlinc.logr.diff."+str(chrom)+"."+str(start)+"."+str(stop)+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()
#########
## now plot 10 examples with lncRNA and as-RT switching
## plot lncRNA rectangles and allele specific RT
#  problem is that smoothing data shifts that compared to on smooth data
# use absolute value of logr diff diff for error bars
standard_deviation = repli_diff["bouha2.bouha10"].std()
print("head")
print(head.drop_duplicates(["chrom","start","stop"]))
for index,row in head.drop_duplicates(["chrom","start","stop"]).iterrows():
    f,ax = plt.subplots(2,1,figsize=(2,4),sharex=False)
    start=row["start"]
    stop=row["stop"]
    chrom=str(row["chrom"])
    plt.suptitle(chrom)
    ax_lnc = ax[0].twinx()
    for index2, row2 in head[(head["chrom"]==chrom) & (head["start"]==start) & (head["stop"]==stop) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.05), width=row2["stop"]-row2["start"], height=0.1,
                     facecolor=row2["color"], edgecolor=row2["color"],hatch="/",fill=False) ## plot vlincs as rectangles
            # rectvar = Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.16,
         #                 facecolor="gray", edgecolor="gray",hatch="/",fill=False) ## plot vlincs as rectangles
        ax_lnc.add_patch(rect)
    ax_lnc.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax_lnc.set_xlim([max(0,start-3000000),stop+3000000])
    ax_lnc.set_ylim([-0.6,0.6])
    ax_lnc.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 6))

    ax[0].set_xlim([max(0,start-3000000),stop+3000000])
    ax[0].set_ylim([-0.6,0.6])
    ax[0].set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 6))
    # (repli_df["sample"]=="bouha.2.repli.5") & (repli_df["chrom"]==chrom) 
    bouha2 = repli_df[(repli_df["sample"]=="bouha.2.repli.5") &(repli_df["chrom"]==chrom)&(repli_df["start"]>=start-4000000) & (repli_df["stop"]<=stop+4000000)]
    bouha10 = repli_df[(repli_df["sample"]=="bouha.10.repli.")& (repli_df["chrom"]==chrom) & (repli_df["start"]>=start-4000000) & (repli_df["stop"]<=stop+4000000)]
    # for index,row in bouha2.iterrows():
    ### shuld smooth this right below
    ax[0].plot(repli_diff[repli_diff["chrom"]==chrom]["start"],
    repli_diff[repli_diff["chrom"]==chrom]["bouha2.bouha10"])
    # plot standard deviation on the logr diff diff plot
    ax[0].fill_between(repli_diff[repli_diff["chrom"]==chrom]["start"],
        repli_diff[repli_diff["chrom"]==chrom]["bouha2.bouha10"]- standard_deviation/2, 
    repli_diff[repli_diff["chrom"]==chrom]["bouha2.bouha10"] + standard_deviation/2,
                 color='blue', alpha=0.1)
    try:
        ax[0].set_ylim([-max(repli_diff[repli_diff["chrom"]==chrom]["bouha2.bouha10"]),
            max(repli_diff[repli_diff["chrom"]==chrom]["bouha2.bouha10"])])
    except:
        ax[0].set_ylim([-2,2])

    ax[1].set_xlim([max(0,start-3000000),stop+3000000])
    ax[1].set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 6))
    ax[1].set_ylim([min(min(bouha2["logr_hap1"]),min(bouha2["logr_hap2"]),min(bouha10["logr_hap1"]),min(bouha10["logr_hap2"])),
        max(max(bouha2["logr_hap1"]),max(bouha2["logr_hap2"]),max(bouha10["logr_hap1"]),max(bouha10["logr_hap2"]))])

    ax[1].plot(bouha2["start"],
    smooth_vector(list(bouha2["start"]),
    list(bouha2["logr_hap1"])),
    c="red",linestyle="--")
    ax[1].plot(bouha2["start"],
    smooth_vector(list(bouha2["start"]),
    list(bouha2["logr_hap2"])),
    c="blue",linestyle="--")
    ax[1].plot(bouha10["start"],
    smooth_vector(list(bouha10["start"]),
    list(bouha10["logr_hap1"])),
    c="red")
    ax[1].plot(bouha10["start"],
    smooth_vector(list(bouha10["start"]),
    list(bouha10["logr_hap2"])),
    c="blue")
    plt.savefig("vlinc.logr.diff.diff."+str(chrom)+"."+str(start)+"."+str(stop)+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()


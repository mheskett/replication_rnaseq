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
from matplotlib.lines import Line2D
from sys import argv
import glob
import statsmodels.api as sm
from sklearn.cluster import KMeans

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
    if len(x) <= 10:
        frac = 1
    elif len(x) >10:
        frac= 6 / len(x)
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

##############
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


repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector

## remove X
repli_df= repli_df[repli_df["chrom"]!="X"]
#############################
#### tasks. plot all repli samples normalized with allele specific and non allele specific. normalized
#### look for sample specific e/l switchers
#### plot individual cell lines non-normalized to look for AS-RT
#### find AS-RT switching
df_normalized_logr_hap1 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap1")).reset_index()
df_normalized_logr_hap2 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap2")).reset_index()
df_normalized_logr_hap1["bouha.2.bouha.10"] = df_normalized_logr_hap1["bouha.2.repli.5"] - df_normalized_logr_hap1["bouha.10.repli."]
df_normalized_logr_hap1["gm12878.4.gm12878.5"] = df_normalized_logr_hap1["gm12878.4.repli"] - df_normalized_logr_hap1["gm12878.5.repli"]
df_normalized_logr_hap2["bouha.2.bouha.10"] = df_normalized_logr_hap2["bouha.2.repli.5"] - df_normalized_logr_hap2["bouha.10.repli."]
df_normalized_logr_hap2["gm12878.4.gm12878.5"] = df_normalized_logr_hap2["gm12878.4.repli"] - df_normalized_logr_hap2["gm12878.5.repli"]
print(df_normalized_logr_hap1)
df_normalized_logr_hap1["arm"] = df_normalized_logr_hap1.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_normalized_logr_hap2["arm"] = df_normalized_logr_hap2.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

df_logr = repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr").reset_index()
df_logr["arm"] = df_logr.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
####
## calculate the sum of the absolute difference between hap1s and hap2s
zscore = lambda x: (x - x.mean()) / x.std()

sum_difference_gm = np.log2(abs(df_normalized_logr_hap1["gm12878.4.gm12878.5"]) + abs(df_normalized_logr_hap2["gm12878.4.gm12878.5"]) )
sum_difference_eb = np.log2(abs(df_normalized_logr_hap1["bouha.2.bouha.10"]) + abs(df_normalized_logr_hap2["bouha.2.bouha.10"]) )
df_normalized_logr_hap1["bouha.epigenetic.difference"] = sum_difference_eb.transform(zscore)
df_normalized_logr_hap1["gm.epigenetic.difference"] = sum_difference_gm.transform(zscore)
###



plt.subplots(figsize=(2,2))
# plt.hist(df_normalized_logr_hap1["bouha.epigenetic.difference"],lw=4,bins=30)
# plt.hist(df_normalized_logr_hap1["gm.epigenetic.difference"],lw=4,bins=30)
sns.kdeplot(df_normalized_logr_hap1["bouha.epigenetic.difference"],cut=0)
sns.kdeplot(df_normalized_logr_hap1["gm.epigenetic.difference"],cut=0)

plt.xlim([-3.5,3.5])
# plt.xticks()
plt.savefig("sum.epigenetic.differences.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
####
color_dict = {"gm12878.4.repli":"plum","gm12878.5.repli":"olivedrab","bouha.10.repli.":"y",
"bouha.2.repli.5":"b"}
color_dict = {"gm12878.4.repli":"red","gm12878.5.repli":"red","bouha.10.repli.":"b",
"bouha.2.repli.5":"b"}
legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10)]
comparisons = ["bouha.2.bouha.10","gm12878.4.gm12878.5"]
####
####

print(df_normalized_logr_hap1)
for i in range(len(chromosomes)):
    f,ax=plt.subplots(figsize=(10,2))
    df_chrom=df_normalized_logr_hap1[df_normalized_logr_hap1["chrom"]==chromosomes[i]]

    ax.plot(df_chrom[df_chrom["arm"]=="p"]["start"],
        smooth_vector(list(df_chrom[df_chrom["arm"]=="p"]["start"]),
        list(abs(df_chrom[df_chrom["arm"]=="p"]["bouha.epigenetic.difference"]))),
        c="red")
    ax.plot(df_chrom[df_chrom["arm"]=="q"]["start"],
        smooth_vector(list(df_chrom[df_chrom["arm"]=="q"]["start"]),
            list(abs(df_chrom[df_chrom["arm"]=="q"]["bouha.epigenetic.difference"]))),
            c="red")
    ax.plot(df_chrom[df_chrom["arm"]=="p"]["start"],
        smooth_vector(list(df_chrom[df_chrom["arm"]=="p"]["start"]),
        list(abs(df_chrom[df_chrom["arm"]=="p"]["gm.epigenetic.difference"]))),
        c="blue")
    ax.plot(df_chrom[df_chrom["arm"]=="q"]["start"],
        smooth_vector(list(df_chrom[df_chrom["arm"]=="q"]["start"]),
            list(abs(df_chrom[df_chrom["arm"]=="q"]["gm.epigenetic.difference"]))),
            c="blue")
    ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
    ax.set_xlim([0,chromosome_length[chromosomes[i]]])
plt.close()
for i in range(len(comparisons)):
    for j in range(len(chromosomes)):
        haplotype1 = df_normalized_logr_hap1[df_normalized_logr_hap1["chrom"]==chromosomes[j]]
        haplotype2 = df_normalized_logr_hap2[df_normalized_logr_hap2["chrom"]==chromosomes[j]]
        f,ax=plt.subplots(1,1, figsize=(10,2) )

        ax.plot(haplotype1[haplotype1["arm"]=="p"]["start"],
            smooth_vector(list(haplotype1[haplotype1["arm"]=="p"]["start"]),
            list(abs(haplotype1[haplotype1["arm"]=="p"][comparisons[i]]))),
            c="red")
        ax.plot(haplotype1[haplotype1["arm"]=="q"]["start"],
            smooth_vector(list(haplotype1[haplotype1["arm"]=="q"]["start"]),
                list(abs(haplotype1[haplotype1["arm"]=="q"][comparisons[i]]))),
                c="red")
        ax.plot(haplotype2[haplotype2["arm"]=="p"]["start"],
            smooth_vector(list(haplotype2[haplotype2["arm"]=="p"]["start"]),
                list(abs(haplotype2[haplotype2["arm"]=="p"][comparisons[i]]))),
            c="blue")
        ax.plot(haplotype2[haplotype2["arm"]=="q"]["start"],
            smooth_vector(list(haplotype2[haplotype2["arm"]=="q"]["start"]),
                list(abs(haplotype2[haplotype2["arm"]=="q"][comparisons[i]]))),
            c="blue")
        ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[j]], 16))
        ax.set_xlim([0,chromosome_length[chromosomes[j]]])
        # for index,row in haplotype1[haplotype1["zscore"]>=2].iterrows():
        #     rect=Rectangle((row["start"], -5), width=row["stop"]-row["start"], height=10,
        #              facecolor=row["repli_color"],fill=True)
        #     ax.add_patch(rect)  
        # for index,row in haplotype2[haplotype2["zscore"]>=2].iterrows():
        #     rect=Rectangle((row["start"], -5), width=row["stop"]-row["start"], height=10,
        #              facecolor=row["repli_color"],fill=True)
        #     ax.add_patch(rect)

        plt.savefig("haplotype.diff."+comparisons[i]+str(chromosomes[j])+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
        plt.close()
plt.close()
############################
# plot two samples hap1 hap2
filenames_repli = filenames_repli[2:]
for j in range(len(chromosomes)):
    f,ax=plt.subplots(1,1, figsize=(10,2) )
    for i in range(len(filenames_repli)):
        haplotype1 = df_normalized_logr_hap1[df_normalized_logr_hap1["chrom"]==chromosomes[j]  ]
        haplotype2 = df_normalized_logr_hap2[df_normalized_logr_hap2["chrom"]==chromosomes[j] ]
        styles=["solid","dotted","dashed","dashdot"]
        # ax.plot(haplotype1[haplotype1["arm"]=="p"]["start"],
        #     smooth_vector(list(haplotype1[haplotype1["arm"]=="p"]["start"]),
        #     list(haplotype1[haplotype1["arm"]=="p"][filenames_repli[i]])),
        #     c="red",linestyle=styles[i])
        # ax.plot(haplotype1[haplotype1["arm"]=="q"]["start"],
        #     smooth_vector(list(haplotype1[haplotype1["arm"]=="q"]["start"]),
        #         list(haplotype1[haplotype1["arm"]=="q"][filenames_repli[i]])),
        #         c="red",linestyle=styles[i])
        ax.plot(haplotype2[haplotype2["arm"]=="p"]["start"],
            smooth_vector(list(haplotype2[haplotype2["arm"]=="p"]["start"]),
                list(haplotype2[haplotype2["arm"]=="p"][filenames_repli[i]])),
            c="blue",linestyle=styles[i])
        ax.plot(haplotype2[haplotype2["arm"]=="q"]["start"],
            smooth_vector(list(haplotype2[haplotype2["arm"]=="q"]["start"]),
                list(haplotype2[haplotype2["arm"]=="q"][filenames_repli[i]])),
            c="blue",linestyle=styles[i])
        ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[j]], 16))
        ax.set_xlim([0,chromosome_length[chromosomes[j]]])       
    plt.show()
    plt.close()



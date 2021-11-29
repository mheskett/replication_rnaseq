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
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed",
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
#############################
#### tasks. plot all repli samples normalized with allele specific and non allele specific. normalized
#### look for sample specific e/l switchers
#### plot individual cell lines non-normalized to look for AS-RT
#### find AS-RT switching
df_normalized_logr = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr")).reset_index()
df_normalized_logr["arm"] = df_normalized_logr.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_logr = repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr").reset_index()
df_logr["arm"] = df_logr.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)


####
color_dict = {"gm12878.4.repli":"plum","gm12878.5.repli":"olivedrab","bouha.10.repli.":"y",
"bouha.2.repli.5":"b"}
color_dict = {"gm12878.4.repli":"red","gm12878.5.repli":"red","bouha.10.repli.":"b",
"bouha.2.repli.5":"b","bouha.3.repli.5":"b","bouha.4.repli.5":"b","bouha.13.repli.":"b","bouha.15.repli.":"b"} ### add full color dict here....l
legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10)]
#### to commpare non allele specific samples against each other before and after normalization
# unnormalized smoothed
## good plot!
for j in range(len(chromosomes)):
    f,ax=plt.subplots(2,1,figsize=(12,4))
    for i in range(len(filenames_repli)):
        tmp = df_logr[(df_logr["chrom"]==chromosomes[j])].loc[:,["chrom","start","stop",filenames_repli[i],"arm"]]
        ax[0].plot(tmp[tmp["arm"]=="p"]["start"],smooth_vector(list(tmp[tmp["arm"]=="p"]["start"]),list(tmp[tmp["arm"]=="p"][filenames_repli[i]])),c=color_dict[filenames_repli[i]])
        ax[0].plot(tmp[tmp["arm"]=="q"]["start"],smooth_vector(list(tmp[tmp["arm"]=="q"]["start"]),list(tmp[tmp["arm"]=="q"][filenames_repli[i]])),c=color_dict[filenames_repli[i]])

    #normalized smoothed
    for i in range(len(filenames_repli)):
        tmp = df_normalized_logr[(df_normalized_logr["chrom"]==chromosomes[j])].loc[:,["chrom","start","stop",filenames_repli[i],"arm"]]
        ax[1].plot(tmp[tmp["arm"]=="p"]["start"],smooth_vector(list(tmp[tmp["arm"]=="p"]["start"]),list(tmp[tmp["arm"]=="p"][filenames_repli[i]])),c=color_dict[filenames_repli[i]])
        ax[1].plot(tmp[tmp["arm"]=="q"]["start"],smooth_vector(list(tmp[tmp["arm"]=="q"]["start"]),list(tmp[tmp["arm"]=="q"][filenames_repli[i]])),c=color_dict[filenames_repli[i]])
    # plt.legend(handles=legend)
    plt.show()
    plt.close()

##################
## try: k-means clustering of logr diff. two clusters total for all autosomes.
###

# try: taking z-scores (log (logr diff +1))
logged_logr_diff_abs = np.log2(repli_df[(repli_df["sample"]==filenames_repli[0]) & (repli_df["logr_diff_abs"]!=0)]["logr_diff_abs"])
zscores = scipy.stats.zscore(repli_df[(repli_df["sample"]==filenames_repli[0]) & (repli_df["logr_diff_abs"]!=0)]["logr_diff_abs"])
plt.hist(logged_logr_diff_abs,bins=20)
plt.show()
## AS repliseq in individual samples
for i in range(len(filenames_repli)):
    tmp= repli_df[(repli_df["sample"]==filenames_repli[i])  & (repli_df["logr_diff_abs"]!=0)]
    ## try log and then zscore
    tmp["zscore"]=list(scipy.stats.zscore(np.log2(tmp[(tmp["chrom"]!="X")]["logr_diff_abs"])))\
    + [3]*len(tmp[tmp["chrom"]=="X"])
    print(tmp[(tmp["zscore"]>=2) * (tmp["chrom"]!="X")])
    ## tried kmeans
    autosome_cluster_labels = KMeans(n_clusters=2).fit_predict(tmp[tmp["chrom"]!="X"]["logr_diff_abs"].values.reshape(-1,1))
    tmp["label"]= list(autosome_cluster_labels) + [2]*len(tmp[tmp["chrom"]=="X"])
    tmp["label"] = tmp["label"].astype(int)
    ###
    col = {0:(1,0,0,0),1:"blue",2:"green"}
    tmp["label_color"] = [col[x] for x in tmp["label"]]
    for j in range(len(chromosomes)):
        tmp_chrom = tmp[tmp["chrom"]==chromosomes[j]]
        f,ax=plt.subplots(1,1, figsize=(12,2) )

        ax.plot(tmp_chrom[tmp_chrom["arm"]=="p"]["start"],
            smooth_vector(list(tmp_chrom[tmp_chrom["arm"]=="p"]["start"]),
            list(tmp_chrom[tmp_chrom["arm"]=="p"]["logr_hap1"])),
            c="red")
        ax.plot(tmp_chrom[tmp_chrom["arm"]=="q"]["start"],
            smooth_vector(list(tmp_chrom[tmp_chrom["arm"]=="q"]["start"]),
                list(tmp_chrom[tmp_chrom["arm"]=="q"]["logr_hap1"])),
                c="red")
        ax.plot(tmp_chrom[tmp_chrom["arm"]=="p"]["start"],
            smooth_vector(list(tmp_chrom[tmp_chrom["arm"]=="p"]["start"]),
                list(tmp_chrom[tmp_chrom["arm"]=="p"]["logr_hap2"])),
            c="blue")
        ax.plot(tmp_chrom[tmp_chrom["arm"]=="q"]["start"],
            smooth_vector(list(tmp_chrom[tmp_chrom["arm"]=="q"]["start"]),
                list(tmp_chrom[tmp_chrom["arm"]=="q"]["logr_hap2"])),
            c="blue")
        # for index,row in tmp_chrom[tmp_chrom["zscore"]>=2].iterrows():
        #     rect=Rectangle((row["start"], -5), width=row["stop"]-row["start"], height=10,
        #              facecolor=row["repli_color"],fill=True)
        #     ax.add_patch(rect)
        ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[j]], 16))
        ax.set_xlim([0, chromosome_length[chromosomes[j]]])
        plt.savefig("as.repli.plot."+str(filenames_repli[i])+str(chromosomes[j])+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
        plt.close()

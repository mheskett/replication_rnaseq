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
import pickle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Shadow
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from sklearn.cluster import KMeans
from scipy import signal
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
def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["hap1_counts"] >= x["hap2_counts"]:
        return x["hap1_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["hap2_counts"]  / x["total_reads"] + 0.5
    return
def intersect_tables(df1,df2):
    ### return all df1 rows that intersect df2 by >0bs
    ### run time is n_squared ish.....dont use this for big DFs. just for small ones
    ### three possibilities of overlap: 1) start is between start/stop. 2) stop is between start/stop. 3) start <= start AND stop >= stop
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    ## add slop to the AS-RT region which is b. a is the lncrnas
    b = pybedtools.BedTool.slop(b,b=250000,g="human_g1k_v37.fasta.fai")
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    return result
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
## load vlincs

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed"]

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
repli_df = repli_df.dropna(how="any",axis="index")
#############################
### we want to find any regions that have high variation of individual haplotypes(outlier of std dev per haplotype?), 
### or high difference between all hap1s and all hap2s (difference of mean?)
df_normalized_logr_hap1 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap1")).reset_index()
df_normalized_logr_hap2 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap2")).reset_index()

## youre accidentally including the chrom coordinates in setd dev calculation here dumbo
df_normalized_logr_hap1["std_dev"] = df_normalized_logr_hap1.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap2["std_dev"] = df_normalized_logr_hap2.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap1["arm"] = df_normalized_logr_hap1.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_normalized_logr_hap2["arm"] = df_normalized_logr_hap2.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

## calculate the sum of the absolute difference between hap1s and hap2s
zscore = lambda x: (x - x.mean()) / x.std()

abs_diff_of_hap_means = abs(df_normalized_logr_hap1.filter(like="bouha",axis=1).mean(axis="columns") - df_normalized_logr_hap2.filter(like="bouha",axis=1).mean(axis="columns"))
df_normalized_logr_hap1["std_dev_zscore"] = df_normalized_logr_hap1["std_dev"].transform(zscore)
df_normalized_logr_hap2["std_dev_zscore"] = df_normalized_logr_hap2["std_dev"].transform(zscore)
df_normalized_logr_hap1["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as below
df_normalized_logr_hap2["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as above

color_dict_repli = {"bouha.10.repli.":"red",
"bouha.2.repli.5":"cyan",
"bouha.3.repli.5":"yellow",
"bouha.4.repli.5":"green",
"bouha.13.repli.":"plum",
"bouha.15.repli.":"olivedrab"}

x1=df_normalized_logr_hap1[df_normalized_logr_hap1["chrom"]=="1"]["bouha.3.repli.5"].values
x2=df_normalized_logr_hap2[df_normalized_logr_hap2["chrom"]=="1"]["bouha.3.repli.5"].values
y1=x1*-1
y2=x2*-1
peaks1, _ = signal.find_peaks(x1,
					height=(None,None),
					threshold=None,
					prominence=0.5,
					distance=2,
					width=1,
					plateau_size=1
					)
valleys1, _ = signal.find_peaks(y1,
					height=(None,None),
					threshold=None,
					prominence=0.5,
					distance=2,
					width=1,
					plateau_size=1
					)
peaks2, _ = signal.find_peaks(x2,
					height=(None,None),
					threshold=None,
					prominence=0.5,
					distance=2,
					width=1,
					plateau_size=1
					)
valleys2, _ = signal.find_peaks(y2,
					height=(None,None),
					threshold=None,
					prominence=0.5,
					distance=2,
					width=1,
					plateau_size=1
					)
plt.plot(x1,c='blue')
plt.plot(x2,c="red")
plt.plot(peaks1, x1[peaks1], "x")
plt.plot(valleys1, -y1[valleys1], "x",c="green")
plt.plot(peaks2, x2[peaks2], "x")
plt.plot(valleys2, -y2[valleys2], "x",c="green")
plt.show()


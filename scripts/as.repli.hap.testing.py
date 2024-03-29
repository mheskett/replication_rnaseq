import csv
import os
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import matplotlib.patheffects as path_effects
import scipy.stats
from scipy.stats import norm
import pickle
import statsmodels.api as sm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Shadow
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
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
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))
df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2 else False,
    axis=1)
df_significant_rows = df[df["significant_deviation"]==True]
df_nonsignificant_rows = df[df["significant_deviation"]==False]
df=df[df["total_reads"]>=10]

###
##############
all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed"]#,
# "/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.500kb.bed"]

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
df_normalized_logr_hap1 = repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap1").add_suffix(".hap1")
df_normalized_logr_hap2 = repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap2").add_suffix(".hap2")
df_normalized_combined = quantile_normalize(pd.concat([df_normalized_logr_hap1,df_normalized_logr_hap2])).reset_index()
df_normalized_combined = df_normalized_combined.set_index(["chrom","start","stop"])
df_normalized_combined["std"] = df_normalized_combined.std(axis=1)
df_normalized_combined["std_gm12878"] = df_normalized_combined.loc[:,[x for x in df_normalized_combined.columns if "gm" in x]].std(axis=1)
df_normalized_combined["std_bouha"] = df_normalized_combined.loc[:,[x for x in df_normalized_combined.columns if "bouha" in x ]].std(axis=1)

zscore = lambda x: (x - x.mean()) / x.std()
df_normalized_combined["std_zscore"] = np.log2(df_normalized_combined["std"]).transform(zscore)
df_normalized_combined["gm_std_zscore"] = np.log2(df_normalized_combined["std_gm12878"]).transform(zscore)
df_normalized_combined["bouha_std_zscore"] = np.log2(df_normalized_combined["std_bouha"]).transform(zscore)

df_normalized_combined = df_normalized_combined.reset_index()
df_normalized_combined["arm"] = df_normalized_combined.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
# sns.kdeplot(df_normalized_combined["std_zscore"],cut=0)
# plt.show()
for j in range(len(chromosomes)):
    f,ax=plt.subplots(figsize=(10,2))
    for i in filenames_repli:
        ax.set_xlim([0,chromosome_length[chromosomes[j]]])
        ax.set_ylim([-4,4])
        ax.plot(df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="p")]["start"],
                    df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="p")][i+".hap1"],
                    c="blue",alpha=1,lw=0.3)
        ax.plot(df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="p")]["start"],
                    df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j])& (df_normalized_combined["arm"]=="p")][i+".hap2"],
                    c="red",alpha=1,lw=0.3)
        ax.plot(df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="q")]["start"],
                    df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="q")][i+".hap1"],
                    c="blue",alpha=1,lw=0.3)
        ax.plot(df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j]) & (df_normalized_combined["arm"]=="q")]["start"],
                    df_normalized_combined[(df_normalized_combined["chrom"]==chromosomes[j])& (df_normalized_combined["arm"]=="q")][i+".hap2"],
                    c="red",alpha=1,lw=0.3)
    for index,row in df_normalized_combined[(df_normalized_combined["gm_std_zscore"]>=2.5) & (df_normalized_combined["chrom"]==chromosomes[j])].iterrows():
        rect=Rectangle((row["start"], -10), width=row["stop"]-row["start"], height=20,
                     facecolor="yellow",fill=True,alpha=0.6)       
        ax.add_patch(rect)
    for index,row in df_normalized_combined[(df_normalized_combined["bouha_std_zscore"]>=2.5) & (df_normalized_combined["chrom"]==chromosomes[j])].iterrows():
        rect=Rectangle((row["start"], -10), width=row["stop"]-row["start"], height=20,
                     facecolor="green",fill=True,alpha=0.6)       
        ax.add_patch(rect)

    plt.savefig("all.samples.each.allele."+str(chromosomes[j])+".png",
            dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()    


all_sig = df_normalized_combined[(df_normalized_combined["gm_std_zscore"]>=2.5) | (df_normalized_combined["bouha_std_zscore"]>=2.5)]
for index,row in all_sig.iterrows():
    f,ax=plt.subplots(figsize=(10,2))
    ax.set_xlim([row["start"]-5000000,row["stop"]+5000000])
    tmp = df_normalized_combined[(df_normalized_combined["start"].between(row["start"]-6000000,row["stop"]+6000000)) & (df_normalized_combined["chrom"]==row["chrom"])]
    ymin = tmp.set_index(["chrom","start","stop"]).min(axis=1).min()
    ymax = tmp.set_index(["chrom","start","stop"]).max(axis=1).max()
    ax.set_ylim([ymin,ymax])
    ax.set_xticks = np.linspace(max(0,row["start"]-5000000),row["stop"]+5000000, 6)

    for i in filenames_repli:
        ax.plot(tmp["start"],
                    tmp[i+".hap1"],
                    c="blue",alpha=1,lw=0.5,linestyle="--" if "bouha" in i else "solid")
        ax.plot(tmp["start"],
                    tmp[i+".hap2"],
                    c="red",alpha=1,lw=0.5,linestyle="--" if "bouha" in i else "solid")
    for index2,row2 in tmp.iterrows():
        if row2["gm_std_zscore"]>=2.5:
            rect=Rectangle((row2["start"], -10), width=row2["stop"]-row2["start"], height=20,
                         facecolor="darkgreen",fill=True,alpha=0.5)       
            ax.add_patch(rect)
        if row2["bouha_std_zscore"]>=2.5:
            rect=Rectangle((row2["start"], -10), width=row2["stop"]-row2["start"], height=20,
                         facecolor="gold",fill=True,alpha=0.5)
            ax.add_patch(rect)
    plt.savefig("all.samples.each.allele."+str(row["chrom"])+"-"+str(row["start"])+".png",
            dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()


#### all samples each allele, but put everything into difference between allele space?

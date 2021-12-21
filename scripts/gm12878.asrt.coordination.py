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
    result["start"] = result["start"].astype(int)
    result["stop"] = result["stop"].astype(int)
    print(result)
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

def remove_blacklisted(df):
    blacklist = pd.read_table("encode.blacklist.final.bed",sep="\t",names=["chrom","start","stop","name","score","strand"])
    blacklist_bed = pybedtools.BedTool.from_dataframe(blacklist)
    df_bed = pybedtools.BedTool.from_dataframe(df)
    result = df_bed.intersect(blacklist_bed,f=0.15,wa=True,v=True).sort(g="human_g1k_v37.fasta.fai")
    result = result.to_dataframe(names=list(df.columns))
    result["chrom"] = result["chrom"].astype(str)
    
    # print(result)
    return result
  
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

### LOAD RTQTLS
rtqtls = pd.read_csv("test123.txt",index_col=False,sep="\t",
    names=["rtQTL","replicatio_variant","chrom","position","-log10_p_value","associated_region_start", 
            "associated_region_end","associate_region_length", "affected_region_start",   
            "affected_region_end", "affected_region_length","distance_to_rep_peak"],
            dtype={"chrom":str,"position":int})
## load vlincs
## loading coding
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))
coding_files=["bouha2.protein.coding.all.counts.bed",
                "bouha3.protein.coding.all.counts.bed",
                "bouha4.protein.coding.all.counts.bed",
                "bouha10.protein.coding.all.counts.bed",
                "bouha13.protein.coding.all.counts.bed",
                "bouha15.protein.coding.all.counts.bed"]

## GM FILES
# coding_files=["gm12878.4.protein.coding.all.counts.bed",
#                 "gm12878.5.protein.coding.all.counts.bed"]
# coding_dfs = []
# for i in range(len(coding_files)):
#     coding_df = pd.read_csv(coding_files[i],sep="\t",
#                             names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
#                             dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
#     coding_df["total_reads"] = coding_df["hap1_counts"] + coding_df["hap2_counts"]
#     coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
#     coding_df["sample"] = coding_files[i][0:9]
#     add_binom_pval(coding_df)
#     coding_dfs += [coding_df]
# df_coding = pd.concat(coding_dfs)
# df_coding = df_coding[df_coding["total_reads"]>=30]
# df_coding = df_coding[df_coding["chrom"]!="X"]
# df_coding["reads_per_kb"] = df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 )
# df_coding  = df_coding[df_coding["reads_per_kb"]>=1]
# df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
#     .reshape(1,-1))*2.5 else False,
#     axis=1)

##### check on FHIT gene for matt
# print("check fhit gene")
# df_coding[(df_coding["chrom"]=="3")& (df_coding["start"]>=56000000) & (df_coding["stop"]<=63000000)].to_csv("fhit.expression.txt",sep="\t")
# peak of 5 reads per kb for coding (informative reads!)
# sns.kdeplot(df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 ),cut=0,clip=(0,100))
# plt.show()

# ## total reads peak is about 100 for coding
# sns.kdeplot(df_coding["total_reads"],cut=0,clip=(0,2000))
# plt.show()
######
#### get all 6 repliseq and lncrna files in this analysis!!!
# vlinc_files=["/Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.bed"]

# # vlinc_files = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.all.bed",
# #                 "/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.all.bed"]
# dfs = []
# for i in range(len(vlinc_files)):
#     df = pd.read_csv(vlinc_files[i],sep="\t",
#                             names= ["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_counts","hap2_counts"],
#                             dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
#     df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
#     df["skew"] = df.apply(helper_func, axis = 1)
#     df["sample"] = os.path.basename(vlinc_files[i])[0:9]
#     add_binom_pval(df)
#     dfs += [df]
# df = pd.concat(dfs)
# unique_genes = list(df["name"].drop_duplicates())
# switchers = [] # list of rows that are switchers
# nonswitchers=[]
# df_significant_rows = df[df["binom_pval"]<=0.001]
# df_nonsignificant_rows = df[df["binom_pval"] >=0.001]
# model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))
# df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2 else False,
#     axis=1)
# df_significant_rows = df[df["significant_deviation"]==True]
# df_nonsignificant_rows = df[df["significant_deviation"]==False]
# df=df[df["total_reads"]>=20]
# # reads per KB
## ok so vlncRNAs have about 0.25 reads median per KB in THIS library
# sns.kdeplot(df["total_reads"] / ((df["stop"] - df["start"]) / 1000 ),cut=0,clip=(0,20))
# plt.show()

# ##median total reads is about 50 with a wide tail.
# sns.kdeplot(df["total_reads"],cut=0,clip=(0,1250))
# plt.show()
###
##############
# all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed"]
# all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.250kb.bed"]

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.250kb.bed"]

filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
repli_li
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
    df_repli["logr_hap1_"+filenames_repli[i]] = np.log2((df_repli["hap1_early"]+1) / (df_repli["hap1_late"]+1))
    df_repli["logr_hap2_"+filenames_repli[i]] = np.log2((df_repli["hap2_early"]+1) / (df_repli["hap2_late"]+1))
    df_repli["low_coverage_"+filenames_repli[i]] = (df_repli["hap1_early"] + df_repli["hap2_early"] + df_repli["hap1_late"] + df_repli["hap2_late"])<=250
    df_repli = df_repli[df_repli["low_coverage_"+filenames_repli[i]]==False]
    repli_li.append(df_repli.loc[:,["logr_hap1_"+filenames_repli[i],"logr_hap2_"+filenames_repli[i]]])

repli_df = pd.concat(repli_li,axis=1)
repli_df_non_normalized = repli_df.reset_index()


# plot before and after normalization
# for i in range(len(all_files_repli)):
#     plt.plot(repli_df_non_normalized[repli_df_non_normalized["chrom"]=="1"]["logr_hap1_"+filenames_repli[i]])
#     plt.plot(repli_df_non_normalized[repli_df_non_normalized["chrom"]=="1"]["logr_hap2_"+filenames_repli[i]])   
# plt.show()

repli_df = quantile_normalize(repli_df.dropna(how="any",axis=0)).reset_index()
repli_df= repli_df[repli_df["chrom"]!="X"]

repli_df = remove_blacklisted(repli_df)
# for i in range(len(all_files_repli)):
#     plt.plot(repli_df[repli_df["chrom"]=="1"]["logr_hap1_"+filenames_repli[i]])
#     plt.plot(repli_df[repli_df["chrom"]=="1"]["logr_hap2_"+filenames_repli[i]])   
# plt.show()

# logr diff raw is positive if hap1 early, negative if hap1 late
for i in range(len(all_files_repli)):
    repli_df["logr_diff_raw_"+filenames_repli[i]] = repli_df["logr_hap1_"+filenames_repli[i]] - repli_df["logr_hap2_"+filenames_repli[i]]# positive if hap1 early, negative if hap2 early
for i in range(len(all_files_repli)):
    repli_df["logr_diff_abs_"+filenames_repli[i]] = abs(repli_df["logr_hap1_"+filenames_repli[i]] - repli_df["logr_hap2_"+filenames_repli[i]])# positive if hap1 early, negative if hap2 early
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
zscore = lambda x: (x - x.mean()) / x.std()
for i in range(len(all_files_repli)):
    repli_df["zscore_logr_diff_abs"+filenames_repli[i]] = repli_df["logr_diff_abs_"+filenames_repli[i]].transform(zscore) # same thing as below

repli_df["std_dev"] = repli_df.filter(like="logr_hap",axis=1).std(axis="columns")

# print(repli_df.filter(like="logr_hap",axis=1).std(axis="columns"))
# sns.kdeplot(repli_df["std_dev"],cut=0)
# plt.show()
# for i in range(len(all_files_repli)):
#     sns.kdeplot((abs(repli_df["logr_diff_raw_"+filenames_repli[i]])),cut=0)
# plt.show()

# for i in range(len(all_files_repli)):
#     sns.kdeplot((abs(repli_df["zscore_logr_diff_abs"+filenames_repli[i]])),cut=0)
# plt.show()

repli_df["any_asrt"] = [True if x>=3 else False for x in repli_df.filter(like="zscore_logr_diff_abs",axis=1).max(axis=1)]

color_dict = {"bouha.4.a":"green",
"bouha.15.":"olivedrab",
"bouha.10.":"red",
"bouha.3.a":"yellow",
"bouha.2.a":"cyan",
"bouha.13.":"plum",
"gm12878.4":"red",
 "gm12878.5":"blue"}

color_dict_coding = {"bouha4.pr":"green",
"bouha15.p":"olivedrab",
"bouha10.p":"red",
"bouha3.pr":"yellow",
"bouha2.pr":"cyan",
"bouha13.p":"plum",
"gm12878.4":"red",
"gm12878.5":"blue"}

color_dict_repli = {"bouha.10.repli.":"red",
"bouha.2.repli.2":"cyan",
"bouha.3.repli.2":"yellow",
"bouha.4.repli.2":"green",
"bouha.13.repli.":"plum",
"bouha.15.repli.":"olivedrab",
"gm12878.4.repli":"red",
"gm12878.5.repli":"blue"}
# #######
# df_coding["color"] = [color_dict_coding[x] for x in df_coding["sample"]]
# df["color"] = [color_dict[x] for x in df["sample"]]
####

## create sum difference hetween all AS1-AS2 pairs, look at distrubtion
## really looking for ASRT in one or more subclones
tmp = repli_df[repli_df["any_asrt"]==True]
tmp = repli_df[repli_df["std_dev"]>=1.25]

a = pybedtools.BedTool.from_dataframe(tmp.drop_duplicates(["chrom","start","stop"]).loc[:,["chrom","start","stop"]])
merged = a.merge().to_dataframe(names=["chrom","start","stop"])
sum_bases = (merged["stop"] - merged["start"]).sum()
print("sumb bases")
print(sum_bases)
# keep.
# for index,row in tmp.drop_duplicates(["chrom","start","stop"]).iterrows():
# # for index,row in thayer_fish_loci.drop_duplicates(["chrom","start","stop"]).iterrows():

#     f, ax = plt.subplots(1,1,figsize=(6,2))
#     plt.rc('xtick', labelsize=3) 
#     plt.rc('ytick', labelsize=8) 
#     start=row["start"]
#     stop=row["stop"]
#     chrom=str(row["chrom"])
#     plt.suptitle(chrom)

#     for i in range(len(filenames_repli)):
#         tmp2 = repli_df[(repli_df["chrom"]==chrom) 
#             &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000)]
#         ax.plot(tmp2["start"],
#                 smooth_vector(tmp2["start"],tmp2["logr_hap1_"+filenames_repli[i]]),
#             c=color_dict_repli[filenames_repli[i]],lw=1)
#         ax.plot(tmp2["start"],
#             smooth_vector(tmp2["start"],tmp2["logr_hap2_"+filenames_repli[i]]),
#             c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1) ## -- line style is haplotype 2
#     ax.set_ylim([-3.5,3.5])
#     ax.set_yticks([-3,-2,-1,0,1,2,3])
#     ax.axhline(y=0,linestyle="--",c="black",lw=0.4)
#     ax.set_xlim([max(0,start-3000000),stop+3000000])
#     ax.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 12))
#     plt.show()


### for each chromosome 
logr_diff_tmp = repli_df.filter(like="zscore_logr_diff_abs",axis=1)

asrt_multi = []
for index,row in logr_diff_tmp.iterrows():
    sorted_vals = row.sort_values(ascending=False)
    if (sorted_vals[0]>=3 and sorted_vals[1]>=3):
        asrt_multi += [True]
    else:
        asrt_multi += [False]
repli_df["asrt_multi"] = asrt_multi

repli_df_asrt_multi_merged = repli_df[repli_df["asrt_multi"]==True]
repli_df_asrt_multi_merged = pybedtools.BedTool.from_dataframe(repli_df_asrt_multi_merged).merge(d=500001).to_dataframe(names=["chrom","start","stop"])
repli_df_asrt_multi_merged["chrom"] = repli_df_asrt_multi_merged["chrom"].astype(str)

for i in range(len(chromosomes)):
    tmp_chrom = repli_df_asrt_multi_merged[(repli_df_asrt_multi_merged["chrom"]==chromosomes[i])]
    subplot_index = 0
    if len(tmp_chrom)==0:
        continue
    plt.rc('xtick', labelsize=3) 
    plt.rc('ytick', labelsize=8) 
    f, (ax,ax_map) = plt.subplots(2,len(tmp_chrom)+1,figsize=(16,2.3),gridspec_kw={'height_ratios': [6, 1]},sharex=False)
    for index,row in tmp_chrom.iterrows():
        plt.rc('xtick', labelsize=3) 
        plt.rc('ytick', labelsize=8) 
        start=row["start"]
        stop=row["stop"]
        chrom=row["chrom"]
        for j in range(len(filenames_repli)):
            ## this is per sample now
            tmp2 = repli_df[(repli_df["chrom"]==chrom) 
                    &(repli_df["start"]>=start-4000000) & (repli_df["stop"]<=stop+4000000)]
            sample_alpha =  1 if len(tmp2[tmp2["zscore_logr_diff_abs"+filenames_repli[j]]>=2.5])>0 else 0.1
            ax[subplot_index].plot(tmp2["start"],
                    smooth_vector(tmp2["start"],tmp2["logr_hap1_"+filenames_repli[j]]),
                     c=color_dict_repli[filenames_repli[j]],lw=1.4,alpha=sample_alpha)
            ax[subplot_index].plot(tmp2["start"],
                    smooth_vector(tmp2["start"],tmp2["logr_hap2_"+filenames_repli[j]]),
                    c=color_dict_repli[filenames_repli[j]],linestyle="--",lw=1.4,alpha=sample_alpha) ## -- line style is haplotype 2            
            ax[subplot_index].set_ylim([-4.5,3.5])
            # ax[subplot_index].set_yticks([-3,-2,-1,0,1,2,3])
            ax[subplot_index].axhline(y=0,linestyle="--",c="black",lw=0.4)
            ax[subplot_index].set_xlim([max(0,start-2000000),stop+2000000])
            ax[subplot_index].set_xticks(np.linspace(max(0,start-2000000),stop+2000000, 6))


        for j in range(len(filenames_repli)):
            xpos = j
            ## widen up the window to consider other samples ASRT, as above
            tmp_exact_pos = repli_df[(repli_df["chrom"]==chrom) 
                    &(repli_df["start"]>=start-500000) & (repli_df["stop"]<=stop+500000)]
            raw_logr_diff_vector = tmp_exact_pos["logr_diff_raw_"+filenames_repli[j]]
            idxmax = abs(raw_logr_diff_vector).idxmax()

            zscore_logr_diff_vector = tmp_exact_pos["zscore_logr_diff_abs"+filenames_repli[j]]
            sample_alpha =  1 if (zscore_logr_diff_vector.max()>=2.5) else 0

            ypos = (0.05,-1) if raw_logr_diff_vector[idxmax]>0 else (-1,0.05)
            rect_hap1=Rectangle((xpos, ypos[0]), width=0.95, height=0.90,
                 facecolor=color_dict_repli[filenames_repli[j]],fill=True,alpha=sample_alpha)
            # rect_hap2=Rectangle((xpos, ypos[1]), width=0.95, height=0.90,
            #      facecolor=color_dict_repli[filenames_repli[j]],fill=True,hatch="////////",alpha=sample_alpha)
            
            ax_map[subplot_index].add_patch(rect_hap1)
            # ax_map[subplot_index].add_patch(rect_hap2)
            ax_map[subplot_index].set_xlim([0,6])
            ax_map[subplot_index].set_ylim([-1,1])

        subplot_index +=1
    plt.suptitle(chromosomes[i])
    plt.savefig("chromosome.coordination.gm12878."+str(chrom)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

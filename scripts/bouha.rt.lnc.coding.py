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
from scipy.stats import truncnorm


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
coding_dfs = []
for i in range(len(coding_files)):
    coding_df = pd.read_csv(coding_files[i],sep="\t",
                            names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    coding_df["total_reads"] = coding_df["hap1_counts"] + coding_df["hap2_counts"]
    coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
    coding_df["sample"] = coding_files[i][0:9]
    add_binom_pval(coding_df)
    coding_dfs += [coding_df]
df_coding = pd.concat(coding_dfs)
df_coding = df_coding[df_coding["total_reads"]>=20]
df_coding = df_coding[df_coding["chrom"]!="X"]
df_coding["reads_per_kb"] = df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 )
df_coding  = df_coding[df_coding["reads_per_kb"]>=1]
df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
    .reshape(1,-1))*2.5 else False,
    axis=1)

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
vlinc_files=["/Users/mike/replication_rnaseq/all.final.data/bouha.2.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.10.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.all.bouha.vlinc.calls.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.all.bouha.vlinc.calls.bed"]

# vlinc_files = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.all.bed",
#                 "/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rep1.rep2.vlincs.all.bed"]
dfs = []
for i in range(len(vlinc_files)):
    df = pd.read_csv(vlinc_files[i],sep="\t",
                            names= ["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
    df["skew"] = df.apply(helper_func, axis = 1)
    df["sample"] = os.path.basename(vlinc_files[i])[0:9]
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
df=df[df["total_reads"]>=20]
# reads per KB
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

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.250kb.bed"]

# all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.250kb.bed"]

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
repli_df["total_reads"] = repli_df["hap1_early"] + repli_df["hap2_early"] + repli_df["hap1_late"] + repli_df["hap2_late"] 
repli_df["low_coverage"] = [True if x <=100 else False for x in repli_df["total_reads"]]
repli_df = repli_df[repli_df["low_coverage"]==False]
# repli_df = repli_df[repli_df["low_coverage"]==False] ## removing low coverage...
## weird peak aroudn 500 reads per window. this will change if RT window size changes.
# sns.kdeplot(repli_df["total_reads"],cut=0,clip=(0,100000))
# plt.show()
## remove X
## plot histogram of coverage for all 6?
# for i in repli_df["sample"].drop_duplicates():
#     sns.kdeplot(repli_df[repli_df["sample"]==i]["total_reads"],cut=0,clip=(0,10000))
# plt.show()

repli_df= repli_df[repli_df["chrom"]!="X"]
repli_df = repli_df.dropna(how="any",axis="index")
#############################
### we want to find any regions that have high variation of individual haplotypes(outlier of std dev per haplotype?), 
### or high difference between all hap1s and all hap2s (difference of mean?)
df_normalized_logr_hap1 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap1")).reset_index()
df_normalized_logr_hap2 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap2")).reset_index()
df_normalized_logr_hap1["std_dev"] = df_normalized_logr_hap1.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap2["std_dev"] = df_normalized_logr_hap2.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap1["arm"] = df_normalized_logr_hap1.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_normalized_logr_hap2["arm"] = df_normalized_logr_hap2.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

####
## calculate the sum of the absolute difference between hap1s and hap2s
zscore = lambda x: (x - x.mean()) / x.std()
abs_diff_of_hap_means = abs(df_normalized_logr_hap1.filter(like="bouha",axis=1).mean(axis="columns") - df_normalized_logr_hap2.filter(like="bouha",axis=1).mean(axis="columns"))
df_normalized_logr_hap1["std_dev_zscore"] = df_normalized_logr_hap1["std_dev"].transform(zscore)
df_normalized_logr_hap2["std_dev_zscore"] = df_normalized_logr_hap2["std_dev"].transform(zscore)
df_normalized_logr_hap1["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as below
df_normalized_logr_hap2["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as above


test = pd.concat([df_normalized_logr_hap1.set_index(["chrom","start","stop"]).add_suffix("_hap1"),
            df_normalized_logr_hap2.set_index(["chrom","start","stop"]).add_suffix("_hap2")],axis=1).filter(like="bouha",axis=1).dropna(how="any",axis="rows")

test["std_dev"] = test.std(axis="columns")
test_locs = tes[test["std_dev"]>=1.06]
# sns.kdeplot(test["std_dev"],cut=0)
# plt.show()
# exit()
# ### just plot all regions with high std dev
# combined = pd.concat([df_normalized_logr_hap1[df_normalized_logr_hap1["std_dev_zscore"]>=2.6],
#                     df_normalized_logr_hap2[df_normalized_logr_hap2["std_dev_zscore"]>=2.6]],axis=0)
#                     # df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3]],axis=0)

# # combined = df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3]
# combined = combined.dropna(how="any",axis="index")
# ## take the combined data frame and cluster each site as a "sample" and each RT value as an observation
# ####


# plt.figure(figsize=(2,2))
# sns.kdeplot(df_normalized_logr_hap1["std_dev"],cut=0,clip=(0,2),lw=2)
# # sns.kdeplot(df_normalized_logr_hap2["std_dev"],cut=0,clip=(0,2),lw=2)
# plt.savefig("bouha.rt.hap.std.dev.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)

print("median of allele std dev across windows: ", np.median(df_normalized_logr_hap1["std_dev"].dropna()))
print("std deviation of allele std dev across windows: ", np.std(df_normalized_logr_hap1["std_dev"]) )

# result=[]
# for index,row in df_normalized_logr_hap1.iterrows():
#     tmp1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==row["chrom"]) & ((df_normalized_logr_hap1["start"]==row["start"]))].dropna(how="any",)
#     tmp2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==row["chrom"]) & ((df_normalized_logr_hap2["start"]==row["start"]))].dropna(how="any")
#     if ~tmp1.empty and ~tmp2.empty:
#         maxi = max(tmp1.loc[:,filenames_repli].max().max(),
#                     tmp2.loc[:,filenames_repli].max().max())
#         mini = min(tmp1.loc[:,filenames_repli].min().min(),
#                     tmp2.loc[:,filenames_repli].min().min())
#         result += [maxi-mini]

# sns.kdeplot(result,cut=0)
# plt.show()
# exit()
color_dict = {"bouha.4.a":"green",
"bouha.15.":"royalblue",
"bouha.10.":"red",
"bouha.3.a":"yellow",
"bouha.2.a":"cyan",
"bouha.13.":"plum",
"gm12878.4":"red",
 "gm12878.5":"blue"}

color_dict_coding = {"bouha4.pr":"green",
"bouha15.p":"royalblue",
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
"bouha.15.repli.":"royalblue",
"gm12878.4.repli":"red",
"gm12878.5.repli":"blue"}
# #######
df_coding["color"] = [color_dict_coding[x] for x in df_coding["sample"]]
df["color"] = [color_dict[x] for x in df["sample"]]
####

## start with ALL repli variation areas, then add lncs and coding geens for all samps
## this will give all sites that have either hap1 or hap2 significant variation
tmp =  pd.concat([df_normalized_logr_hap1[(df_normalized_logr_hap1["std_dev"]>=1)],
                        df_normalized_logr_hap2[df_normalized_logr_hap2["std_dev"]>=1]],axis=0)

tmp = remove_blacklisted(tmp)



tmp_merged_bed = pybedtools.BedTool.from_dataframe(tmp.drop_duplicates(["chrom","start","stop"]).loc[:,["chrom","start","stop"]])
tmp_merged = tmp_merged_bed.merge(d=250001).to_dataframe(names=["chrom","start","stop"])
tmp_merged["chrom"] = tmp_merged["chrom"].astype(str)
print("number of merged windows with >3 std dev ",len(tmp_merged))
exit()
# tmp =  df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3]

tmp=tmp.dropna(how="any",axis="index")
tmp.to_csv("bouha.vert.merged.txt",sep="\t",index=False,header=True)
df = df.dropna(how="any",axis="index")
df_coding = df_coding.dropna(how="any",axis="index")
tmp_lnc = intersect_tables(df,tmp)
tmp_coding = intersect_tables(df_coding,tmp) ## could add more slop here 

thayer_fish_loci = pd.read_csv("thayer.fish.loci.txt",sep="\t",
                        names=["chrom","start","stop"],
                        dtype={"chrom":str,"start":int,"stop":int})
tmp_coding[tmp_coding["significant_deviation"]==True]["name"].drop_duplicates().to_csv("coding.genes.in.vert.regions.bouha.txt",sep="\t",index=False,header=False)



#### chromosome 1, just do the STD DEV of haplotypes for all samps

for j in range(len(chromosomes)):
    f,ax = plt.subplots(figsize=(10,1))
    hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chromosomes[j])]
    hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chromosomes[j])]
    
    ax.axhline(y=1,linestyle="--",lw=0.5,c="black")
    ax.plot(hap1[hap1["arm"]=="p"]["start"],
            smooth_vector(hap1[hap1["arm"]=="p"]["start"],
                        np.maximum(hap1[hap1["arm"]=="p"]["std_dev"],hap2[hap2["arm"]=="p"]["std_dev"])),c="black",lw=2)
    
    ax.plot(hap1[hap1["arm"]=="q"]["start"],
            smooth_vector(hap1[hap1["arm"]=="q"]["start"],
                            np.maximum(hap1[hap1["arm"]=="q"]["std_dev"],hap2[hap2["arm"]=="q"]["std_dev"])),c="black",lw=2)

    for index3,row3 in tmp_merged[tmp_merged["chrom"]==chromosomes[j]].iterrows():
        rect=Rectangle((row3["start"]-250000, 0), width=row3["stop"]-row3["start"]+500000, height=5,
                 facecolor="lightgray",alpha=1,fill=True)
        ax.add_patch(rect)
    ax.set_ylim([0,1.8])
    ax.set_yticks([0,.5,1,1.5,2])
    ax.set_xlim([0,chromosome_length[chromosomes[j]]])
    ax.set_xticks(np.linspace(0,chromosome_length[chromosomes[j]],16))
    plt.savefig("bouha.vert."+str(chromosomes[j])+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()   
################

#### STD DEV of each haplotype for all samples. chromosome view. red and blue for matt
for j in range(len(chromosomes)):
    f,ax = plt.subplots(figsize=(10,1))
    hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chromosomes[j])]
    hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chromosomes[j])]
    
    ax.axhline(y=1,linestyle="--",lw=0.5,c="black")
    ax.plot(hap1[hap1["arm"]=="p"]["start"],
            smooth_vector(hap1[hap1["arm"]=="p"]["start"],hap1[hap1["arm"]=="p"]["std_dev"]),c="red",lw=1)
    ax.plot(hap1[hap1["arm"]=="q"]["start"],
            smooth_vector(hap1[hap1["arm"]=="q"]["start"],hap1[hap1["arm"]=="q"]["std_dev"]),c="red",lw=1)
    ax.plot(hap2[hap2["arm"]=="p"]["start"],
            smooth_vector(hap2[hap2["arm"]=="p"]["start"],hap2[hap2["arm"]=="p"]["std_dev"]),c="blue",lw=1)
    ax.plot(hap2[hap2["arm"]=="q"]["start"],
            smooth_vector(hap2[hap2["arm"]=="q"]["start"],hap2[hap2["arm"]=="q"]["std_dev"]),c="blue",lw=1)
    for index3,row3 in tmp_merged[tmp_merged["chrom"]==chromosomes[j]].iterrows():
        rect=Rectangle((row3["start"]-250000, 0), width=row3["stop"]-row3["start"]+500000, height=5,
                 facecolor="lightgray",alpha=1,fill=True)
        ax.add_patch(rect)
    ax.set_ylim([0,1.8])
    ax.set_yticks([0,.5,1,1.5,2])
    ax.set_xlim([0,chromosome_length[chromosomes[j]]])
    ax.set_xticks(np.linspace(0,chromosome_length[chromosomes[j]],16))
    plt.savefig("bouha.vert.hap"+str(chromosomes[j])+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()   
exit()
# print("starting repli and asar plotting loops")
plt.rc('xtick', labelsize=3) 
plt.rc('ytick', labelsize=8) 

####### keep all this. just commenting out for now.
for index,row in tmp_merged.drop_duplicates(["chrom","start","stop"]).iterrows():
# for index,row in thayer_fish_loci.drop_duplicates(["chrom","start","stop"]).iterrows():
    plt.rc('xtick', labelsize=3) 
    plt.rc('ytick', labelsize=8) 
    f, (ax,ax_peaks) = plt.subplots(2,1,figsize=(2,2.3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1]})

    start=row["start"]
    stop=row["stop"]
    chrom=str(row["chrom"])
    plt.suptitle(chrom)
    ax_lnc = ax.twinx()
    # df is lncs
    for index2, row2 in df[(df["chrom"]==chrom) & (df["start"]>=start-2000000) & (df["stop"]<=stop+2000000) 
                            & (df["significant_deviation"]==True) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.0125), width=row2["stop"]-row2["start"], height=0.025,
                     facecolor=row2["color"], edgecolor="black",fill=True,lw=.4)
        # shadow = Shadow(rect, 10000,-0.0015 )                             
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    for index2, row2 in df[(df["chrom"]==chrom) & (df["start"]>=start-2000000) & (df["stop"]<=stop+2000000) 
                        & (df["significant_deviation"]==False) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.0125), width=row2["stop"]-row2["start"], height=0.025,
                 facecolor=row2["color"], edgecolor="black",fill=True,lw=.4,alpha=0.1)
        # shadow = Shadow(rect, 10000,-0.0015 )                             
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    # df coding is coding genes
    for index5, row5 in df_coding[(df_coding["chrom"]==chrom) & (df_coding["start"]>=start-2000000) & (df_coding["stop"]<=stop+2000000)
                                 & (df_coding["significant_deviation"]==True) ].iterrows():
        rect=Rectangle((row5["start"], row5["skew"]-.0125), width=row5["stop"]-row5["start"], height=0.025,
                     facecolor=row5["color"], edgecolor="black",fill=True,lw=.6,linestyle="dotted")
        # shadow = Shadow(rect, 10000,-0.0015 )                                   
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    for index5, row5 in df_coding[(df_coding["chrom"]==chrom) & (df_coding["start"]>=start-2000000) & (df_coding["stop"]<=stop+2000000)
                                 & (df_coding["significant_deviation"]==False) ].iterrows():
        rect=Rectangle((row5["start"], row5["skew"]-.0125), width=row5["stop"]-row5["start"], height=0.025,
                     facecolor=row5["color"], edgecolor="black",fill=True,lw=.6,alpha=0.1,linestyle="dotted")
        # shadow = Shadow(rect, 10000,-0.0015 )                                   
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    # print("done adding coding")
    # ax_lnc.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax_lnc.set_xlim([max(0,start-3000000),stop+3000000])
    ax_lnc.set_ylim([-0.52,0.52])
    ax_lnc.set_yticks([-0.5,-.25,0,.25,.5])
    # ax_lnc.set_xticks(np.linspace(max(0,start-6000000),stop+6000000, 12))
    ### fix the gray shadowing
    for index3,row3 in df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) & (df_normalized_logr_hap1["std_dev"]>=1)].iterrows():
        rect=Rectangle((row3["start"]-250000, -10), width=row3["stop"]-row3["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)
    for index4,row4 in df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom) & (df_normalized_logr_hap2["std_dev"]>=1)].iterrows():
        rect=Rectangle((row4["start"]-250000, -10), width=row4["stop"]-row4["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)

    # print("done adding asars, shadows, and background highlights")
    hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) 
            &(df_normalized_logr_hap1["start"]>=start-5000000) & (df_normalized_logr_hap1["stop"]<=stop+5000000) ]
    hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom)  
            &(df_normalized_logr_hap2["start"]>=start-5000000) & (df_normalized_logr_hap2["stop"]<=stop+5000000)]
    # print("merely subsetting into hap1 and hap2")
    ax.set_xlim([max(0,start-3000000),stop+3000000])
    ax.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 12))
    # print(hap1)
    # ax.set_ylim([min(hap1.filter(like="bouha",axis=1).values.min(),hap2.filter(like="bouha",axis=1).values.min()),
    #     max(hap1.filter(like="bouha",axis=1).values.max(),hap2.filter(like="bouha",axis=1).values.max())])
    ax.set_ylim([-3.5,3.5])
    ax.set_yticks([-3,-2,-1,0,1,2,3])
    ax.axhline(y=0,linestyle="--",c="black",lw=0.4)
    # print("how could you be getting stuck before this and after the previous....")
    #### normalized repliseq
    ## this needs to be looped and colors fixed?
    for i in range(len(filenames_repli)):
        # print("starting to plot ", filenames_repli[i])
        # print("printing hap1 df", hap1)
        # print("printing hap2 df", hap2)
        ax.plot(hap1["start"],
                smooth_vector(hap1["start"],hap1[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],lw=1)
        ax.plot(hap2["start"],
            smooth_vector(hap2["start"],hap2[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1) ## -- line style is haplotype 2
        # print("plotted sample ",filenames_repli[i])

    ### now plot the variance curve in the subplot
    ax_peaks.plot(hap1["start"],
                np.maximum(hap1["std_dev"],hap2["std_dev"]),
                c="black",
                lw=0.8)
    ax_peaks.set_ylim([0,2])
    ax_peaks.set_yticks([0,1,2])
    ax_peaks.set_xlim([max(0,start-3000000),stop+3000000])
    ax_peaks.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 12))  
    rtqtl_label=False
    rtqtls_to_plot = rtqtls[(rtqtls["chrom"]==chrom) & (rtqtls["position"] >= start-4000000) & (rtqtls["position"]<=stop+4000000)]
    if len(rtqtls_to_plot)>0:
        rtqtl_label = True
        print("rtqtl on chromosome ",chrom,start)
        ax.text(rtqtls_to_plot["affected_region_start"],3.1,"rtQTL region",fontdict = {'family': 'serif','color':  'black','weight': 'normal','size': 6})
        ax.hlines(y=3,xmin=rtqtls_to_plot["affected_region_start"],xmax=rtqtls_to_plot["affected_region_end"],linewidth=1,color="black",zorder=10)
    
    ## plot thayer fish loci
    # ax.text(row["start"],-3.4,"FISH probe",fontdict = {'family': 'serif','color':  'black','weight': 'normal','size': 6})
    # ax.hlines(y=-3,xmin=row["start"],xmax=row["stop"],linewidth=1,color="black",zorder=10)

    if rtqtl_label==True:
        plt.savefig("fig4.epigenetic.coding.lnc.rt.rtQTL."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    else:
        plt.savefig("fig4.epigenetic.coding.lnc.rt."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)        
    plt.close()

####
### making shared figs
####### keep all this. just commenting out for now.
shared_vert = pd.read_table("bouha_gm12878_shared_vert_loci.txt",sep="\t",header=0,index_col=False)
for index,row in shared_vert.iterrows():
# for index,row in thayer_fish_loci.drop_duplicates(["chrom","start","stop"]).iterrows():
    plt.rc('xtick', labelsize=3) 
    plt.rc('ytick', labelsize=8) 
    f, (ax,ax_peaks) = plt.subplots(2,1,figsize=(2,2.3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1]})

    start=row["start"]
    stop=row["stop"]
    chrom=str(row["chrom"])
    plt.suptitle(chrom)
    ax_lnc = ax.twinx()
    # df is lncs
    for index2, row2 in df[(df["chrom"]==chrom) & (df["start"]>=start-2000000) & (df["stop"]<=stop+2000000) 
                            & (df["significant_deviation"]==True) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.0125), width=row2["stop"]-row2["start"], height=0.025,
                     facecolor=row2["color"], edgecolor="black",fill=True,lw=.4)
        # shadow = Shadow(rect, 10000,-0.0015 )                             
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    for index2, row2 in df[(df["chrom"]==chrom) & (df["start"]>=start-2000000) & (df["stop"]<=stop+2000000) 
                        & (df["significant_deviation"]==False) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.0125), width=row2["stop"]-row2["start"], height=0.025,
                 facecolor=row2["color"], edgecolor="black",fill=True,lw=.4,alpha=0.1)
        # shadow = Shadow(rect, 10000,-0.0015 )                             
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    # df coding is coding genes
    for index5, row5 in df_coding[(df_coding["chrom"]==chrom) & (df_coding["start"]>=start-2000000) & (df_coding["stop"]<=stop+2000000)
                                 & (df_coding["significant_deviation"]==True) ].iterrows():
        rect=Rectangle((row5["start"], row5["skew"]-.0125), width=row5["stop"]-row5["start"], height=0.025,
                     facecolor=row5["color"], edgecolor="black",fill=True,lw=.6,linestyle="dotted")
        # shadow = Shadow(rect, 10000,-0.0015 )                                   
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    for index5, row5 in df_coding[(df_coding["chrom"]==chrom) & (df_coding["start"]>=start-2000000) & (df_coding["stop"]<=stop+2000000)
                                 & (df_coding["significant_deviation"]==False) ].iterrows():
        rect=Rectangle((row5["start"], row5["skew"]-.0125), width=row5["stop"]-row5["start"], height=0.025,
                     facecolor=row5["color"], edgecolor="black",fill=True,lw=.6,alpha=0.1,linestyle="dotted")
        # shadow = Shadow(rect, 10000,-0.0015 )                                   
        # ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    # print("done adding coding")
    # ax_lnc.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax_lnc.set_xlim([max(0,start-3000000),stop+3000000])
    ax_lnc.set_ylim([-0.52,0.52])
    ax_lnc.set_yticks([-0.5,-.25,0,.25,.5])
    # ax_lnc.set_xticks(np.linspace(max(0,start-6000000),stop+6000000, 12))
    ### fix the gray shadowing
    for index3,row3 in df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) & (df_normalized_logr_hap1["std_dev"]>=1)].iterrows():
        rect=Rectangle((row3["start"]-500000, -10), width=row3["stop"]-row3["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)
    for index4,row4 in df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom) & (df_normalized_logr_hap2["std_dev"]>=1)].iterrows():
        rect=Rectangle((row4["start"]-500000, -10), width=row4["stop"]-row4["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)

    # print("done adding asars, shadows, and background highlights")
    hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) 
            &(df_normalized_logr_hap1["start"]>=start-5000000) & (df_normalized_logr_hap1["stop"]<=stop+5000000) ]
    hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom)  
            &(df_normalized_logr_hap2["start"]>=start-5000000) & (df_normalized_logr_hap2["stop"]<=stop+5000000)]
    # print("merely subsetting into hap1 and hap2")
    ax.set_xlim([max(0,start-3000000),stop+3000000])
    ax.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 12))
    # print(hap1)
    # ax.set_ylim([min(hap1.filter(like="bouha",axis=1).values.min(),hap2.filter(like="bouha",axis=1).values.min()),
    #     max(hap1.filter(like="bouha",axis=1).values.max(),hap2.filter(like="bouha",axis=1).values.max())])
    ax.set_ylim([-3.5,3.5])
    ax.set_yticks([-3,-2,-1,0,1,2,3])
    ax.axhline(y=0,linestyle="--",c="black",lw=0.4)
    # print("how could you be getting stuck before this and after the previous....")
    #### normalized repliseq
    ## this needs to be looped and colors fixed?
    for i in range(len(filenames_repli)):
        # print("starting to plot ", filenames_repli[i])
        # print("printing hap1 df", hap1)
        # print("printing hap2 df", hap2)
        ax.plot(hap1["start"],
                smooth_vector(hap1["start"],hap1[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],lw=1)
        ax.plot(hap2["start"],
            smooth_vector(hap2["start"],hap2[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1) ## -- line style is haplotype 2
        # print("plotted sample ",filenames_repli[i])

    ### now plot the variance curve in the subplot
    ax_peaks.plot(hap1["start"],
                np.maximum(hap1["std_dev"],hap2["std_dev"]),
                c="black",
                lw=0.8)
    ax_peaks.set_ylim([0,2])
    ax_peaks.set_yticks([0,1,2])
    ax_peaks.set_xlim([max(0,start-3000000),stop+3000000])
    ax_peaks.set_xticks(np.linspace(max(0,start-3000000),stop+3000000, 12))  
    rtqtl_label=False
    rtqtls_to_plot = rtqtls[(rtqtls["chrom"]==chrom) & (rtqtls["position"] >= start-4000000) & (rtqtls["position"]<=stop+4000000)]
    if len(rtqtls_to_plot)>0:
        rtqtl_label = True
        print("rtqtl on chromosome ",chrom,start)
        ax.text(rtqtls_to_plot["affected_region_start"],3.1,"rtQTL region",fontdict = {'family': 'serif','color':  'black','weight': 'normal','size': 6})
        ax.hlines(y=3,xmin=rtqtls_to_plot["affected_region_start"],xmax=rtqtls_to_plot["affected_region_end"],linewidth=1,color="black",zorder=10)
    
    ## plot thayer fish loci
    # ax.text(row["start"],-3.4,"FISH probe",fontdict = {'family': 'serif','color':  'black','weight': 'normal','size': 6})
    # ax.hlines(y=-3,xmin=row["start"],xmax=row["stop"],linewidth=1,color="black",zorder=10)

    if rtqtl_label==True:
        plt.savefig("bouha.shared.vert.epigenetic.coding.lnc.rt.rtQTL."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    else:
        plt.savefig("bouha.shared.vert.epigenetic.coding.lnc.rt."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)        
    plt.close()



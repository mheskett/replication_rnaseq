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


def remove_blacklisted(df):
    blacklist = pd.read_table("encode.blacklist.final.bed",sep="\t",names=["chrom","start","stop","name","score","strand"])
    blacklist_bed = pybedtools.BedTool.from_dataframe(blacklist)
    df_bed = pybedtools.BedTool.from_dataframe(df)
    result = df_bed.intersect(blacklist_bed,f=0.15,wa=True,v=True).sort(g="human_g1k_v37.fasta.fai")
    result = result.to_dataframe(names=list(df.columns))
    result["chrom"] = result["chrom"].astype(str)
    
    # print(result)
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


all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.250kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.250kb.bed"]

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
abs_diff_of_hap_means = abs(df_normalized_logr_hap1.filter(like="bouha|gm12878",axis=1).mean(axis="columns") - df_normalized_logr_hap2.filter(like="bouha|gm12878",axis=1).mean(axis="columns"))
df_normalized_logr_hap1["std_dev_zscore"] = df_normalized_logr_hap1["std_dev"].transform(zscore)
df_normalized_logr_hap2["std_dev_zscore"] = df_normalized_logr_hap2["std_dev"].transform(zscore)
df_normalized_logr_hap1["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as below
df_normalized_logr_hap2["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as above

bouha_vert = pd.read_table("bouha.vert.merged.txt",header=0)
gm12878_vert = pd.read_table("gm12878.vert.merged.txt",header=0)

color_dict_repli = {"bouha.10.repli.":"red",
"bouha.2.repli.2":"cyan",
"bouha.3.repli.2":"yellow",
"bouha.4.repli.2":"green",
"bouha.13.repli.":"plum",
"bouha.15.repli.":"olivedrab",
"gm12878.4.repli":"red",
"gm12878.5.repli":"blue"}

bouha_gm12878_shared_vert = intersect_tables(bouha_vert,gm12878_vert)
bouha_gm12878_shared_vert_loci = pybedtools.BedTool.from_dataframe(bouha_gm12878_shared_vert).merge(d=500001).to_dataframe(names=["chrom","start","stop"])

bouha_gm12878_shared_vert_loci.to_csv("bouha_gm12878_shared_vert_loci.txt",sep="\t",index=False,header=True)



    # for index,row in bouha_gm12878_shared_vert_loci.iterrows():
    #     f, (ax,ax_map) = plt.subplots(3,1,figsize=(2.3,16),gridspec_kw={'height_ratios': [3, 3, 1]},sharex=True)
    #     plt.rc('xtick', labelsize=3) 
    #     plt.rc('ytick', labelsize=8) 
    #     start=row["start"]
    #     stop=row["stop"]
    #     chrom=row["chrom"]
    #     for j in range(len(filenames_repli[0:6])):
    #         ## this is per sample now
    #         tmp2 = repli_df[(repli_df["chrom"]==chrom) 
    #                 &(repli_df["start"]>=start-4000000) & (repli_df["stop"]<=stop+4000000)]
    #         ax[0].plot(tmp2["start"],
    #                 smooth_vector(tmp2["start"],tmp2["logr_hap1_"+filenames_repli[j]]),
    #                  c=color_dict_repli[filenames_repli[j]],lw=1.4,alpha=sample_alpha)
    #         ax[0].plot(tmp2["start"],
    #                 smooth_vector(tmp2["start"],tmp2["logr_hap2_"+filenames_repli[j]]),
    #                 c=color_dict_repli[filenames_repli[j]],linestyle="--",lw=1.4,alpha=sample_alpha) ## -- line style is haplotype 2            
    #         ax[0].set_ylim([-3.5,3.5])
    #         # ax[subplot_index].set_yticks([-3,-2,-1,0,1,2,3])
    #         ax[0].axhline(y=0,linestyle="--",c="black",lw=0.4)
    #         ax[0].set_xlim([max(0,start-2000000),stop+2000000])
    #         ax[0].set_xticks(np.linspace(max(0,start-2000000),stop+2000000, 6))
        
    #     for j in range(len(filenames_repli[6:8])):
    #         ## this is per sample now
    #         tmp2 = repli_df[(repli_df["chrom"]==chrom) 
    #                 &(repli_df["start"]>=start-4000000) & (repli_df["stop"]<=stop+4000000)]
    #         ax[1].plot(tmp2["start"],
    #                 smooth_vector(tmp2["start"],tmp2["logr_hap1_"+filenames_repli[j]]),
    #                  c=color_dict_repli[filenames_repli[j]],lw=1.4,alpha=sample_alpha)
    #         ax[1].plot(tmp2["start"],
    #                 smooth_vector(tmp2["start"],tmp2["logr_hap2_"+filenames_repli[j]]),
    #                 c=color_dict_repli[filenames_repli[j]],linestyle="--",lw=1.4,alpha=sample_alpha) ## -- line style is haplotype 2            
    #         ax[1].set_ylim([-3.5,3.5])
    #         # ax[subplot_index].set_yticks([-3,-2,-1,0,1,2,3])
    #         ax[1].axhline(y=0,linestyle="--",c="black",lw=0.4)
    #         ax[1].set_xlim([max(0,start-2000000),stop+2000000])
    #         ax[1].set_xticks(np.linspace(max(0,start-2000000),stop+2000000, 6))

    #     for j in range(len(filenames_repli[0:6])):
    #         xpos = j
    #         ## widen up the window to consider other samples ASRT, as above
    #         tmp_exact_pos = repli_df[(repli_df["chrom"]==chrom) 
    #                 &(repli_df["start"]>=start-500000) & (repli_df["stop"]<=stop+500000)]
    #         raw_logr_diff_vector = tmp_exact_pos["logr_diff_raw_"+filenames_repli[j]]
    #         idxmax = abs(raw_logr_diff_vector).idxmax()

    #         zscore_logr_diff_vector = tmp_exact_pos["zscore_logr_diff_abs"+filenames_repli[j]]
    #         sample_alpha =  1 if (zscore_logr_diff_vector.max()>=2.5) else 0

    #         ypos = (0.05,-1) if raw_logr_diff_vector[idxmax]>0 else (-1,0.05)
    #         rect_hap1=Rectangle((xpos, ypos[0]), width=0.95, height=0.90,
    #              facecolor=color_dict_repli[filenames_repli[j]],fill=True,alpha=sample_alpha)
    #         # rect_hap2=Rectangle((xpos, ypos[1]), width=0.95, height=0.90,
    #         #      facecolor=color_dict_repli[filenames_repli[j]],fill=True,hatch="////////",alpha=sample_alpha)
            
    #         ax_map[0].add_patch(rect_hap1)
    #         # ax_map[subplot_index].add_patch(rect_hap2)
    #         ax_map[0].set_xlim([0,6])
    #         ax_map[0].set_ylim([-1,1])
    #     for j in range(len(filenames_repli[6:8])):
    #         xpos = j
    #         ## widen up the window to consider other samples ASRT, as above
    #         tmp_exact_pos = repli_df[(repli_df["chrom"]==chrom) 
    #                 &(repli_df["start"]>=start-500000) & (repli_df["stop"]<=stop+500000)]
    #         raw_logr_diff_vector = tmp_exact_pos["logr_diff_raw_"+filenames_repli[j]]
    #         idxmax = abs(raw_logr_diff_vector).idxmax()

    #         zscore_logr_diff_vector = tmp_exact_pos["zscore_logr_diff_abs"+filenames_repli[j]]
    #         sample_alpha =  1 if (zscore_logr_diff_vector.max()>=2.5) else 0

    #         ypos = (0.05,-1) if raw_logr_diff_vector[idxmax]>0 else (-1,0.05)
    #         rect_hap1=Rectangle((xpos, ypos[0]), width=0.95, height=0.90,
    #              facecolor=color_dict_repli[filenames_repli[j]],fill=True,alpha=sample_alpha)
    #         # rect_hap2=Rectangle((xpos, ypos[1]), width=0.95, height=0.90,
    #         #      facecolor=color_dict_repli[filenames_repli[j]],fill=True,hatch="////////",alpha=sample_alpha)
            
    #         ax_map[0].add_patch(rect_hap1)
    #         # ax_map[subplot_index].add_patch(rect_hap2)
    #         ax_map[0].set_xlim([0,6])
    #         ax_map[0].set_ylim([-1,1])
    # plt.suptitle(chromosomes[i])
    # plt.savefig("shared.vert."+str(chrom)+str(start)+".png",
    #     dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    # plt.close()

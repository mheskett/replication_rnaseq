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

#### get all 6 repliseq and lncrna files in this analysis!!!
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
print(df_significant_rows)
###
##############
all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.3.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.13.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.15.repli.500kb.bed"]
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
repli_df = repli_df.dropna(how="any",axis="index")
#############################
### we want to find any regions that have high variation of individual haplotypes(outlier of std dev per haplotype?), 
### or high difference between all hap1s and all hap2s (difference of mean?)
df_normalized_logr_hap1 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap1")).reset_index()
df_normalized_logr_hap2 = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="logr_hap2")).reset_index()

print("df_normalizedlogrihap1",df_normalized_logr_hap1)
## youre accidentally including the chrom coordinates in setd dev calculation here dumbo
df_normalized_logr_hap1["std_dev"] = df_normalized_logr_hap1.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap2["std_dev"] = df_normalized_logr_hap2.filter(like="bouha",axis=1).std(axis="columns")
df_normalized_logr_hap1["arm"] = df_normalized_logr_hap1.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
df_normalized_logr_hap2["arm"] = df_normalized_logr_hap2.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)


plt.subplots(figsize=(2,2))
sns.kdeplot(df_normalized_logr_hap1["std_dev"],cut=0)
plt.savefig("std.dev.hap1.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()

plt.subplots(figsize=(2,2))
sns.kdeplot(df_normalized_logr_hap2["std_dev"],cut=0)
plt.savefig("std.dev.hap2.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()

####
## calculate the sum of the absolute difference between hap1s and hap2s
zscore = lambda x: (x - x.mean()) / x.std()

# sum_difference_gm = abs(df_normalized_logr_hap1["gm12878.4.gm12878.5"]) + abs(df_normalized_logr_hap2["gm12878.4.gm12878.5"])
# sum_difference_eb = abs(df_normalized_logr_hap1["bouha.2.bouha.10"]) + abs(df_normalized_logr_hap2["bouha.2.bouha.10"]) 
# df_normalized_logr_hap1["bouha.epigenetic.difference.raw"] = sum_difference_eb
# df_normalized_logr_hap1["gm.epigenetic.difference.raw"] = sum_difference_gm
abs_diff_of_hap_means = abs(df_normalized_logr_hap1.filter(like="bouha",axis=1).mean(axis="columns") - df_normalized_logr_hap2.filter(like="bouha",axis=1).mean(axis="columns"))
df_normalized_logr_hap1["std_dev_zscore"] = df_normalized_logr_hap1["std_dev"].transform(zscore)
df_normalized_logr_hap2["std_dev_zscore"] = df_normalized_logr_hap2["std_dev"].transform(zscore)
df_normalized_logr_hap1["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as below
df_normalized_logr_hap2["zscore_abs_diff_hap_means"] = abs_diff_of_hap_means.transform(zscore) # same thing as above

# df_normalized_logr_hap1["gm.epigenetic.difference"] = sum_difference_gm.transform(zscore)
###

# plt.subplots(figsize=(2,2))
# sns.kdeplot(df_normalized_logr_hap1["std_dev_zscore"],c="blue",cut=0)
# plt.xlim([-3.5,3.5])
# plt.savefig("std.dev.zscores1.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.close()
# ####
# plt.subplots(figsize=(2,2))
# sns.kdeplot(df_normalized_logr_hap2["std_dev_zscore"],c="orange",cut=0)
# plt.xlim([-3.5,3.5])
# plt.savefig("std.dev.zscores2.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.close()
# ####
# plt.subplots(figsize=(2,2))
# sns.kdeplot(df_normalized_logr_hap1["zscore_abs_diff_hap_means"],c="blue",cut=0)
# plt.xlim([-3.5,3.5])
# plt.savefig("zscore_abs_diff_means.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.close()

color_dict_repli = {"bouha.10.repli.":"red",
"bouha.2.repli.5":"cyan",
"bouha.3.repli.5":"yellow",
"bouha.4.repli.5":"green",
"bouha.13.repli.":"plum",
"bouha.15.repli.":"olivedrab"}
### just plot all regions with high std dev
combined = pd.concat([df_normalized_logr_hap1[df_normalized_logr_hap1["std_dev_zscore"]>=2.6],
                    df_normalized_logr_hap2[df_normalized_logr_hap2["std_dev_zscore"]>=2.6]],axis=0)
                    # df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3]],axis=0)


# combined = df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3]
print("combined",combined)
print(combined[combined.isna()])## dono why so many nans......
print(df_normalized_logr_hap1[df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=3])
combined = combined.dropna(how="any",axis="index")
## take the combined data frame and cluster each site as a "sample" and each RT value as an observation
all_cluster_labels = []
for index,row in combined.sort_values("std_dev_zscore",ascending=False).iterrows():
    X = row.filter(like="bouha").to_numpy().reshape(-1,1)
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
    all_cluster_labels +=[kmeans.labels_]

# sns.clustermap(all_cluster_labels,col_cluster=True,row_cluster=False,figsize=(6,2))
# plt.show()
### also just try clustering samples heirarchically based on features
# df_normalized_logr_hap1 = df_normalized_logr_hap1.dropna(how="any",axis="index")
# sns.clustermap(df_normalized_logr_hap1[df_normalized_logr_hap1["std_dev_zscore"]>=2.6].filter(like="bouha").transpose(),
#     col_cluster=False,
#     row_cluster=True,
#     cmap="RdBu_r",vmax=2,vmin=-2,
#     figsize=(6,2))
# plt.show()
# plt.close()
# plt.subplots(figsize=(12,4))
# df_normalized_logr_hap2 = df_normalized_logr_hap2.dropna(how="any",axis="index")
# sns.clustermap(df_normalized_logr_hap2[df_normalized_logr_hap2["std_dev_zscore"]>=2.6].filter(like="bouha").transpose(),
#     col_cluster=False,
#     row_cluster=True,
#     cmap="RdBu_r",vmax=2,vmin=-2,
#     figsize=(12,4))
# plt.show()
# plt.close()


######
##### keep all this!!!

# for index,row in combined.sort_values("std_dev_zscore",ascending=False).iterrows():
#     start=row["start"]
#     stop=row["stop"]
#     chrom=row["chrom"]
#     # print("start: ",start,"stop",stop,"chrom",chrom)
#     f,ax=plt.subplots(figsize=(10,2))
#     for i in range(len(filenames_repli)):
#         hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) 
#                 &(df_normalized_logr_hap1["start"]>=start-6000000) & (df_normalized_logr_hap1["stop"]<=stop+6000000) ]
#         hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom)  
#             &(df_normalized_logr_hap2["start"]>=start-6000000) & (df_normalized_logr_hap2["stop"]<=stop+6000000)]
#         ax.plot(hap1["start"],
#                 smooth_vector(hap1["start"],hap1[filenames_repli[i]]),
#             c=color_dict_repli[filenames_repli[i]],lw=1,zorder=1)
#         ax.plot(hap2["start"],
#             smooth_vector(hap2["start"],hap2[filenames_repli[i]]),
#             c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1,zorder=2)

#         rect=Rectangle((start-250000, -10), width=stop-start+500000, height=20,
#                  facecolor="lightgray",alpha=1,fill=True) 
#         ax.add_patch(rect)
#         ax.set_xlim([max(0,start-6000000),stop+6000000])
#         ax.set_xticks(np.linspace(max(0,start-6000000),stop+6000000, 12))
#         ax.set_ylim([-3.5,3.5])
#         plt.suptitle("chromosome: "+str(chrom)+". Std. dev. 1: "+str( df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) 
#                                                                     & (df_normalized_logr_hap1["start"]==start)]["std_dev"].values )
#                         +"std. dev 2: "+ str( df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom) 
#                                                                     & (df_normalized_logr_hap2["start"]==start)]["std_dev"].values))
#     plt.savefig("bouha.all.repliseq.outliers."+str(chrom)+"."+str(start)+ ".png",
#             dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
#     plt.close()



####
color_dict = {"bouha.4.":"green",
"bouha.15":"olivedrab",
"bouha.10":"red",
"bouha.3.":"yellow",
"bouha.2.":"cyan",
"bouha.13":"plum"}

color_dict_repli = {"bouha.10.repli.":"red",
"bouha.2.repli.5":"cyan",
"bouha.3.repli.5":"yellow",
"bouha.4.repli.5":"green",
"bouha.13.repli.":"plum",
"bouha.15.repli.":"olivedrab"}
# #######
tmp = intersect_tables(df[df["significant_deviation"]==True],
                df_normalized_logr_hap1[(df_normalized_logr_hap1["std_dev_zscore"]>=2) | (df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=2)])#.sort_values(["bouha.epigenetic.difference1"],ascending=False)

tmp2 = intersect_tables(df[df["significant_deviation"]==True],
                df_normalized_logr_hap2[df_normalized_logr_hap2["std_dev_zscore"]>=2])

tmp = pd.concat([tmp,tmp2])
print("tmp table length ",len(tmp))
# tmp.to_csv("testing.txt",sep="\t")
tmp["color"]=[color_dict[x] for x in tmp["sample"]] # color dict for rna files

print("starting repli and asar plotting loops")
for index,row in tmp.drop_duplicates(["chrom","start","stop"]).iterrows():
    f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
    plt.rc('xtick', labelsize=4) 
    plt.rc('ytick', labelsize=8) 
    start=row["start"]
    stop=row["stop"]
    chrom=str(row["chrom"])
    plt.suptitle(chrom)
    ax_lnc = ax.twinx()
    for index2, row2 in tmp[(tmp["chrom"]==chrom) & (tmp["start"]==start) & (tmp["stop"]==stop) ].iterrows():
        rect=Rectangle((row2["start"], row2["skew"]-.05), width=row2["stop"]-row2["start"], height=0.05,
                     facecolor=row2["color"], edgecolor="black",fill=True,lw=.5)
        shadow = Shadow(rect, 10000,-0.0015 )                                        
        ax_lnc.add_patch(shadow)
        ax_lnc.add_patch(rect)
    print("done adding asars")
    # ax_lnc.axhline(y=0,linestyle="--",lw=0.4,c="black")
    ax_lnc.set_xlim([max(0,start-4000000),stop+4000000])
    ax_lnc.set_ylim([-0.6,0.6])
    ax_lnc.set_yticks([-0.5,-.25,0,.25,.5])
    # ax_lnc.set_xticks(np.linspace(max(0,start-6000000),stop+6000000, 12))
    ### fix the gray shadowing
    for index3,row3 in df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) & (df_normalized_logr_hap1["std_dev_zscore"]>=2)
                                                | (df_normalized_logr_hap1["zscore_abs_diff_hap_means"]>=2)].iterrows():
        rect=Rectangle((row3["start"]-250000, -10), width=row3["stop"]-row3["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)
    for index4,row4 in df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom) & (df_normalized_logr_hap2["std_dev_zscore"]>=2)].iterrows():
        rect=Rectangle((row4["start"]-250000, -10), width=row4["stop"]-row4["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)
    print("done adding asars, shadows, and background highlights")
    hap1 = df_normalized_logr_hap1[(df_normalized_logr_hap1["chrom"]==chrom) 
            &(df_normalized_logr_hap1["start"]>=start-4000000) & (df_normalized_logr_hap1["stop"]<=stop+4000000) ]
    hap2 = df_normalized_logr_hap2[(df_normalized_logr_hap2["chrom"]==chrom)  
            &(df_normalized_logr_hap2["start"]>=start-4000000) & (df_normalized_logr_hap2["stop"]<=stop+4000000)]
    print("merely subsetting into hap1 and hap2")
    ax.set_xlim([max(0,start-4000000),stop+4000000])
    ax.set_xticks(np.linspace(max(0,start-4000000),stop+4000000, 12))
    # print(hap1)
    # ax.set_ylim([min(hap1.filter(like="bouha",axis=1).values.min(),hap2.filter(like="bouha",axis=1).values.min()),
    #     max(hap1.filter(like="bouha",axis=1).values.max(),hap2.filter(like="bouha",axis=1).values.max())])
    ax.set_ylim([-3.5,3.5])
    ax.axhline(y=0,linestyle="--",c="black",lw=0.4)
    print("how could you be getting stuck before this and after the previous....")
    #### normalized repliseq
    ## this needs to be looped and colors fixed?
    for i in range(len(filenames_repli)):
        print("starting to plot ",filenames_repli[i])
        print("printing hap1 df",hap1)
        print("printing hap2 df",hap2)
        ax.plot(hap1["start"],
                smooth_vector(hap1["start"],hap1[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],lw=1,zorder=1)
        ax.plot(hap2["start"],
            smooth_vector(hap2["start"],hap2[filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1,zorder=2)
        print("plotted sample ",filenames_repli[i])
    plt.savefig("fig4.epigenetic.diff3.all.bouha."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

## counter


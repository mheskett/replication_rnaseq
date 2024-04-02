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
from scipy.stats import ttest_ind

### Add two-sided binomial p-value and FDR-BH q-value to a data frame
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

### Calculate the percent allelic expression. Hap1 is positive and hap2 is negative just for convenience.
def helper_func(x):
    if x["total_reads"]==0:
        return 0
    elif x["hap1_counts"] >= x["hap2_counts"]:
        return x["hap1_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["hap2_counts"]  / x["total_reads"] + 0.5
    return

### Calls bedtools intersect and returns pandas data frame with the correct data types for each column
def intersect_tables(df1,df2):
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    b = pybedtools.BedTool.slop(b,b=250000,g="human_g1k_v37.fasta.fai")
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    result["start"] = result["start"].astype(int)
    result["stop"] = result["stop"].astype(int)
    return result

### Adds the arm information for every genomic position. Requires the cytoband file that details the locatino of the arms.
def get_arms(cytoband):
	arm_dict = {}
	for i in range(len(chromosomes)):
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

### Quantile normalization for repli-seq
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

### Get the total number of bases in a dataframe
def sum_bases(df):
    length = df["stop"] - df["start"]

    return length.sum()

### Use Lowess smoothing for repli-seq regions
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

### Remove the ENCODE blacklisted regions. Requires the ENCODE blacklist bed file
def remove_blacklisted(df):
    blacklist = pd.read_table("encode.blacklist.final.bed",sep="\t",names=["chrom","start","stop","name","score","strand"])
    blacklist_bed = pybedtools.BedTool.from_dataframe(blacklist)
    df_bed = pybedtools.BedTool.from_dataframe(df)
    result = df_bed.intersect(blacklist_bed,f=0.15,wa=True,v=True).sort(g="human_g1k_v37.fasta.fai")
    result = result.to_dataframe(names=list(df.columns))
    result["chrom"] = result["chrom"].astype(str)
    return result
  
### List of chromosome names as strings
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
### List of chromosome arm names.
arms = ["p","q"]
#### Cytoband data to skip over centromeres and gaps				
cytoband = pd.read_table("cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)

### Chromosome length dictionary of hg19 chromosomes
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

### rtQTLs from Koren et al 2014, Cell.
rtqtls = pd.read_csv("test123.txt",index_col=False,sep="\t",
    names=["rtQTL","replicatio_variant","chrom","position","-log10_p_value","associated_region_start", 
            "associated_region_end","associate_region_length", "affected_region_start",   
            "affected_region_end", "affected_region_length","distance_to_rep_peak"],
            dtype={"chrom":str,"position":int})

## Loads the regression model to account for heteroskedastic RNA-seq data
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))

### List of allele specific read counts of coding genes for all subclones
coding_files=["bouha2.protein.coding.all.counts.bed",
                "bouha3.protein.coding.all.counts.bed",
                "bouha4.protein.coding.all.counts.bed",
                "bouha10.protein.coding.all.counts.bed",
                "bouha13.protein.coding.all.counts.bed",
                "bouha15.protein.coding.all.counts.bed"]


## Reads in allele-specific coding gene expression files and performs filtering and stats.
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
df_coding = df_coding[df_coding["total_reads"]>=30]
df_coding = df_coding[df_coding["chrom"]!="X"]
df_coding["reads_per_kb"] = df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 )
df_coding  = df_coding[df_coding["reads_per_kb"]>=1]
df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
    .reshape(1,-1))*2.5 else False,
    axis=1)
df_coding["significant_deviation"] = (df_coding["significant_deviation"]==True) & (df_coding["fdr_pval"]<=0.01)

## List of vlinc call files for all subclones
vlinc_files=["bouha.2.all.bouha.vlinc.calls.bed",
"bouha.3.all.bouha.vlinc.calls.bed",
"bouha.4.all.bouha.vlinc.calls.bed",
"bouha.10.all.bouha.vlinc.calls.bed",
"bouha.13.all.bouha.vlinc.calls.bed",
"bouha.15.all.bouha.vlinc.calls.bed"]

## Reads in allele-specific vlinc expression files and performs filtering. 
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
df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
    axis=1)
df["significant_deviation"] = (df["significant_deviation"]==True) & (df["fdr_pval"]<=0.01) & (df["chrom"]!="X")
df=df[df["total_reads"]>=20]

##### Find subset of genes that switch their allelic bias from hap1 to hap2 or hap2 to hap1 between subclones
df_significant_rows = df[(df["significant_deviation"]==True) & (df["fdr_pval"]<=0.01)]
df_nonsignificant_rows = df[(df["significant_deviation"]==False) | (df["fdr_pval"]>=0.01)]
unique_genes = list(df["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
for i in range(len(unique_genes)):
    samples = df_significant_rows[df_significant_rows["name"]==unique_genes[i]]
    # samples = samples[(samples["binom_pval_plus"]<=0.05) | (samples["binom_pval_minus"] <=0.05)]
    if len(samples)<=1:
        continue
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
nonswitchers = pd.concat(nonswitchers)


### switchers algorithm for coding_
unique_genes_coding = list(df_coding["name"].drop_duplicates())
coding_switchers = [] 
coding_nonswitchers=[]
df_coding_significant_rows = df_coding[df_coding["significant_deviation"]==True]
for i in range(len(unique_genes_coding)):
    samples = df_coding_significant_rows[df_coding_significant_rows["name"]==unique_genes_coding[i]]
    if len(samples)<=1:
        continue
    hap1_skew,hap2_skew= False,False
    for index,row in samples.iterrows():
        if (row["skew"]>=0.1):
            hap1_skew = True
        if (row["skew"]<=-0.1):
            hap2_skew = True
    if hap1_skew and hap2_skew:
        coding_switchers += [samples]
    elif hap1_skew ^ hap2_skew:
        coding_nonswitchers += [samples]
coding_switchers = pd.concat(coding_switchers)
coding_nonswitchers = pd.concat(coding_nonswitchers)

## List of allele specific repli-seq bed files
all_files_repli = ["bouha.10.repli.250kb.bed",
"bouha.2.repli.250kb.bed",
"bouha.3.repli.250kb.bed",
"bouha.4.repli.250kb.bed",
"bouha.13.repli.250kb.bed",
"bouha.15.repli.250kb.bed"]
filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []

## Load and filter repli-seq
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
repli_df = quantile_normalize(repli_df.dropna(how="any",axis=0)).reset_index()
repli_df= repli_df[repli_df["chrom"]!="X"]
repli_df = remove_blacklisted(repli_df)

# logr diff raw is positive if hap1 early, negative if hap1 late
# Add standard deviation for repliseq hap1+hap2, hap1 only, and hap2 onl
# Create genome wide distribution of repli-seq standard deviation per window for all haps
for i in range(len(all_files_repli)):
    repli_df["logr_diff_raw_"+filenames_repli[i]] = repli_df["logr_hap1_"+filenames_repli[i]] - repli_df["logr_hap2_"+filenames_repli[i]]# positive if hap1 early, negative if hap2 early
for i in range(len(all_files_repli)):
    repli_df["logr_diff_abs_"+filenames_repli[i]] = abs(repli_df["logr_hap1_"+filenames_repli[i]] - repli_df["logr_hap2_"+filenames_repli[i]])# positive if hap1 early, negative if hap2 early
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
zscore = lambda x: (x - x.mean()) / x.std()
repli_df["std_dev"] = repli_df.filter(like="logr_hap",axis=1).std(axis="columns")
repli_df["std_dev_hap1"] = repli_df.filter(like="logr_hap1_",axis=1).std(axis="columns")
repli_df["std_dev_hap2"] = repli_df.filter(like="logr_hap2_",axis=1).std(axis="columns")
logr_diff_abs_mean = repli_df.filter(like="logr_diff_abs").mean()
logr_diff_abs_std_dev = repli_df.filter(like="logr_diff_abs").std()
## all haps
mean_std_dev = repli_df[repli_df["chrom"]!="X"]["std_dev"].mean()
std_std_dev = repli_df[repli_df["chrom"]!="X"]["std_dev"].std()
threshold = mean_std_dev + 2.5*std_std_dev
## hap1
mean_std_dev1 = repli_df[repli_df["chrom"]!="X"]["std_dev_hap1"].mean()
std_std_dev1 = repli_df[repli_df["chrom"]!="X"]["std_dev_hap1"].std()
threshold1 = mean_std_dev1 + 2.5*std_std_dev1
###hap2
mean_std_dev2 = repli_df[repli_df["chrom"]!="X"]["std_dev_hap2"].mean()
std_std_dev2 = repli_df[repli_df["chrom"]!="X"]["std_dev_hap2"].std()
threshold2 = mean_std_dev2 + 2.5*std_std_dev2

######
tmp = repli_df[repli_df["std_dev"]>=threshold]
tmp_merged_bed = pybedtools.BedTool.from_dataframe(tmp.drop_duplicates(["chrom","start","stop"]).loc[:,["chrom","start","stop"]])
tmp_merged = tmp_merged_bed.merge(d=250001).to_dataframe(names=["chrom","start","stop"])
tmp_merged["chrom"] = tmp_merged["chrom"].astype(str)
tmp=tmp.dropna(how="any",axis="index")
####
switcher_vert_table = intersect_tables(switchers.drop_duplicates(["name"]),tmp_merged)
coding_switcher_vert_table = intersect_tables(coding_switchers.drop_duplicates(["name"]),tmp_merged)
dae_tl_vert = intersect_tables(df[(df["significant_deviation"]==True)& (df["chrom"]!="X")].drop_duplicates(["chrom","start","stop"]),tmp_merged)
dae_coding_vert = intersect_tables(df_coding[(df_coding["significant_deviation"]==True) & (df_coding["chrom"]!="X")].drop_duplicates(["chrom","start","stop"]),tmp_merged)
df = df.dropna(how="any",axis="index")
df_coding = df_coding.dropna(how="any",axis="index")
tmp_lnc = intersect_tables(df,tmp)
tmp_coding = intersect_tables(df_coding,tmp) ## could add more slop here 
tmp_coding[tmp_coding["significant_deviation"]==True]["name"].drop_duplicates().to_csv("coding.genes.in.vert.regions.bouha.txt",sep="\t",index=False,header=False)


#### PLOTTING
#### Various plots that were included in data analysis/interpretation as well as used for the manuscript.
### Create color for plotting
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

#### ASRT view
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(7,1,figsize=(8,1.5),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1,1,1,1]})
    for i in range(len(filenames_repli)):
        hap1 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like="hap1",axis=1).reset_index()
        hap2 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like="hap2",axis=1).reset_index()
        ax[0].axhline(y=0,linestyle="--",lw=0.5,c="black")
        ax[0].plot(hap1[hap1["arm"]=="p"]["start"],
                    smooth_vector(hap1[hap1["arm"]=="p"]["start"],hap1[hap1["arm"]=="p"]["logr_hap1_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],lw=0.5)
        ax[0].plot(hap2[hap2["arm"]=="p"]["start"],
                smooth_vector(hap2[hap2["arm"]=="p"]["start"],hap2[hap2["arm"]=="p"]["logr_hap2_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(hap1[hap1["arm"]=="q"]["start"],
                    smooth_vector(hap1[hap1["arm"]=="q"]["start"],hap1[hap1["arm"]=="q"]["logr_hap1_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],lw=0.5)
        ax[0].plot(hap2[hap2["arm"]=="q"]["start"],
                smooth_vector(hap2[hap2["arm"]=="q"]["start"],hap2[hap2["arm"]=="q"]["logr_hap2_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2
        ax[0].set_ylim([-3.7,3.5])
        ax[0].set_yticks([-3,-2,-1,0,1,2,3])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["logr_diff_abs_"+filenames_repli[i]]>=(logr_diff_abs_mean["logr_diff_abs_"+filenames_repli[i]] + 2.5*logr_diff_abs_std_dev["logr_diff_abs_"+filenames_repli[i]]))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="red" if row3["logr_diff_raw_"+filenames_repli[i]]>=0 else "blue",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[i+1].add_patch(rect)

        #####
        ax[i+1].plot(repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="p")]["start"],
                    repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="p")]["logr_diff_abs_"+filenames_repli[i]],
                    c="black",
                    lw=0.7)
        ax[i+1].plot(repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="q")]["start"],
                    repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="q")]["logr_diff_abs_"+filenames_repli[i]],
                    c="black",
                    lw=0.7)

        ax[i+1].set_ylim([0,2.3])
        ax[i+1].set_yticks([0,1,2])
        ax[i+1].set_xlim([0,chromosome_length[chrom]])
        ax[i+1].set_xticks([]) 


    plt.savefig("bouha.all.clones.asrt."+str(chrom)+".alleles.png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()


## whole view. RT only.
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, (ax,ax_peaks,ax_peaks2,ax_peaks3) = plt.subplots(4,1,figsize=(8,1.5),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1]})
    for i in range(len(filenames_repli)):
        hap1 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like="hap1",axis=1).reset_index()
        hap2 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like="hap2",axis=1).reset_index()
        ax.axhline(y=0,linestyle="--",lw=0.5,c="black")
        ax.plot(hap1[hap1["arm"]=="p"]["start"],
                    smooth_vector(hap1[hap1["arm"]=="p"]["start"],hap1[hap1["arm"]=="p"]["logr_hap1_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],lw=0.5)
        ax.plot(hap2[hap2["arm"]=="p"]["start"],
                smooth_vector(hap2[hap2["arm"]=="p"]["start"],hap2[hap2["arm"]=="p"]["logr_hap2_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax.plot(hap1[hap1["arm"]=="q"]["start"],
                    smooth_vector(hap1[hap1["arm"]=="q"]["start"],hap1[hap1["arm"]=="q"]["logr_hap1_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],lw=0.5)
        ax.plot(hap2[hap2["arm"]=="q"]["start"],
                smooth_vector(hap2[hap2["arm"]=="q"]["start"],hap2[hap2["arm"]=="q"]["logr_hap2_"+filenames_repli[i]]),
                c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev"]>=threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks.add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev"]>=threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks.add_patch(rect)
        ######
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev_hap1"]>=threshold1)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks2.add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev_hap1"]>=threshold1)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks2.add_patch(rect)
        #####
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev_hap2"]>=threshold2)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks3.add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev_hap2"]>=threshold2)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="lightgray",alpha=1,fill=True)
            ax_peaks3.add_patch(rect)
        #####
        ax_peaks.plot(repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="p")]["start"],
                    repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="p")]["std_dev"],
                    c="black",
                    lw=0.7)
        ax_peaks.plot(repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="q")]["start"],
                    repli_df[(repli_df["chrom"]==chrom) & (repli_df["arm"]=="q")]["std_dev"],
                    c="black",
                    lw=0.7)

        ax_peaks2.plot(hap1[hap1["arm"]=="p"]["start"],
                    hap1[hap1["arm"]=="p"]["std_dev_hap1"],
                    c="black",
                    lw=0.7)
        ax_peaks2.plot(hap1[hap1["arm"]=="q"]["start"],
                    hap1[hap1["arm"]=="q"]["std_dev_hap1"],
                    c="black",
                    lw=0.7)

        ax_peaks3.plot(hap2[hap2["arm"]=="p"]["start"],
                    hap2[hap2["arm"]=="p"]["std_dev_hap2"],
                    c="black",
                    linestyle="--",
                    lw=0.7)
        ax_peaks3.plot(hap2[hap2["arm"]=="q"]["start"],
                    hap2[hap2["arm"]=="q"]["std_dev_hap2"],
                    c="black",
                    linestyle="--",
                    lw=0.7)
        ax_peaks.set_ylim([0,2.3])
        ax_peaks.set_yticks([0,1,2])
        ax_peaks.set_xlim([0,chromosome_length[chrom]])
        ax_peaks.set_xticks([]) 

        ax_peaks2.set_ylim([0,2])
        ax_peaks2.set_yticks([0,1,2])
        ax_peaks2.set_xlim([0,chromosome_length[chrom]])
        ax_peaks2.set_xticks([]) 

        ax_peaks3.set_ylim([0,2])
        ax_peaks3.set_yticks([0,1,2])
        ax_peaks3.set_xlim([0,chromosome_length[chrom]])
        ax_peaks3.set_xticks(np.linspace(0,chromosome_length[chrom],16)) 

        ax.set_ylim([-3.7,3.5])
        ax.set_yticks([-3,-2,-1,0,1,2,3])
        ax.set_xlim([0,chromosome_length[chrom]])
        ax.set_xticks([])
    plt.savefig("bouha.all.clones."+str(chrom)+".alleles.png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()


######
for index,row in tmp_merged.drop_duplicates(["chrom","start","stop"]).iterrows():
# for index,row in thayer_fish_loci.drop_duplicates(["chrom","start","stop"]).iterrows():
    plt.rc('xtick', labelsize=3) 
    plt.rc('ytick', labelsize=8) 
    f, (ax,ax_peaks) = plt.subplots(2,1,figsize=(2,2.3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1]})

#### for rtqtls
    # chrom=row["chrom"]
    # start=row["affected_region_start"]
    # stop=row["affected_region_end"]
    ####### 
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
    for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev"]>=threshold)].iterrows():
        rect=Rectangle((row3["start"]-250000, -10), width=row3["stop"]-row3["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)

    # print("done adding asars, shadows, and background highlights")
    hap1 = repli_df[(repli_df["chrom"]==chrom) 
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000) ].set_index(["chrom","start","stop"]).filter(like="hap1",axis=1).reset_index()
    hap2 = repli_df[(repli_df["chrom"]==chrom)  
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000)].set_index(["chrom","start","stop"]).filter(like="hap2",axis=1).reset_index()
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
                smooth_vector(hap1["start"],hap1["logr_hap1_"+filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],lw=1)
        ax.plot(hap2["start"],
            smooth_vector(hap2["start"],hap2["logr_hap2_"+filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1) ## -- line style is haplotype 2
        # print("plotted sample ",filenames_repli[i])

    ### now plot the variance curve in the subplot

    tmp_std = repli_df[(repli_df["chrom"]==chrom) 
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000) ]
    ax_peaks.plot(tmp_std["start"],
                tmp_std["std_dev"],
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
###3
#### FOR SHARED VERT
shared_vert = pd.read_csv("bouha_gm12878_shared_vert_loci.txt",sep="\t",header=0,dtype={"chrom":str,'start':int,"stop":int})
for index,row in shared_vert.iterrows():
####
    plt.rc('xtick', labelsize=3) 
    plt.rc('ytick', labelsize=8) 
    f, (ax,ax_peaks) = plt.subplots(2,1,figsize=(2,2.3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1]})

#### for rtqtls
    # chrom=row["chrom"]
    # start=row["affected_region_start"]
    # stop=row["affected_region_end"]
    ####### 
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
    for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & (repli_df["std_dev"]>=threshold)].iterrows():
        rect=Rectangle((row3["start"]-250000, -10), width=row3["stop"]-row3["start"]+500000, height=20,
                 facecolor="lightgray",alpha=1,fill=True) 
        ax.add_patch(rect)

    # print("done adding asars, shadows, and background highlights")
    hap1 = repli_df[(repli_df["chrom"]==chrom) 
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000) ].set_index(["chrom","start","stop"]).filter(like="hap1",axis=1).reset_index()
    hap2 = repli_df[(repli_df["chrom"]==chrom)  
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000)].set_index(["chrom","start","stop"]).filter(like="hap2",axis=1).reset_index()
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
                smooth_vector(hap1["start"],hap1["logr_hap1_"+filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],lw=1)
        ax.plot(hap2["start"],
            smooth_vector(hap2["start"],hap2["logr_hap2_"+filenames_repli[i]]),
            c=color_dict_repli[filenames_repli[i]],linestyle="--",lw=1) ## -- line style is haplotype 2
        # print("plotted sample ",filenames_repli[i])

    ### now plot the variance curve in the subplot

    tmp_std = repli_df[(repli_df["chrom"]==chrom) 
            &(repli_df["start"]>=start-5000000) & (repli_df["stop"]<=stop+5000000) ]
    ax_peaks.plot(tmp_std["start"],
                tmp_std["std_dev"],
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
        plt.savefig("shared.epigenetic.coding.lnc.rt.rtQTL."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    else:
        plt.savefig("shared.epigenetic.coding.lnc.rt."+str(chrom)+"."+str(start)+ ".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)        
    plt.close()



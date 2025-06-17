import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.stats import norm
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact



def sum_region_length(df):
    diffs = df["stop"] - df["start"]
    return diffs.sum()


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


chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]



arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.hg38.txt",sep="\t",
                            names =["chrom","start","stop","arm","band"])

arm_dict = get_arms(cytoband)

chromosome_length = {"chr1":    248956422,
"chr2":    242193529,
"chr3":    198295559,
"chr4":    190214555,
"chr5":    181538259,
"chr6":    170805979,
"chr7":    159345973,
"chrX":    156040895,
"chr8":    145138636,
"chr9":    138394717,
"chr11":   135086622,
"chr10":   133797422,
"chr12":   133275309,
"chr13":   114364328,
"chr14":   107043718,
"chr15":   101991189,
"chr16":   90338345,
"chr17":   83257441,
"chr18":   80373285,
"chr20":   64444167,
"chr19":   58617616,
"chrY":    57227415,
"chr22":   50818468,
"chr21":   46709983}

ratios = {"chr1":    248956422/248956422,
"chr2":    242193529/248956422,
"chr3":    198295559/248956422,
"chr4":    190214555/248956422,
"chr5":    181538259/248956422,
"chr6":    170805979/248956422,
"chr7":    159345973/248956422,
"chrX":    156040895/248956422,
"chr8":    145138636/248956422,
"chr9":    138394717/248956422,
"chr11":   135086622/248956422,
"chr10":   133797422/248956422,
"chr12":   133275309/248956422,
"chr13":   114364328/248956422,
"chr14":   107043718/248956422,
"chr15":   101991189/248956422,
"chr16":   90338345/248956422,
"chr17":   83257441/248956422,
"chr18":   80373285/248956422,
"chr20":   64444167/248956422,
"chr19":   58617616/248956422,
"chrY":    57227415/248956422,
"chr22":   50818468/248956422,
"chr21":   46709983/248956422}


## TLs
df_eb_tls = pd.read_csv("gm.as.tl.counts.all.txt",sep="\t")

## repli
df_eb_rt = pd.read_csv("../combined.rt.data.set/gm.rt.txt",sep='\t',dtype={"gm_std_dev_both_haps":float})
vert_rt_windows = df_eb_rt[(df_eb_rt["gm_vert"]==True) & (df_eb_rt["chrom"]!="chrX")].loc[:,["chrom","start","stop"]].drop_duplicates()
# vert_rt_windows_bed = pybedtools.BedTool.from_dataframe(vert_rt_windows).sort()#.merge()
all_rt_windows_bed = pybedtools.BedTool.from_dataframe(df_eb_rt[df_eb_rt["chrom"]!="chrX"]).sort().merge()
all_rt_windows_bed.to_dataframe(disable_auto_names=True, header=None).to_csv("all_rt_windows.bed",sep="\t",header=None,index=None)

# make another one without merging
all_rt_windows_bed_unmerged = pybedtools.BedTool.from_dataframe(df_eb_rt[(df_eb_rt["chrom"]!="chrX") ].loc[:,["chrom","start","stop","gm_std_dev_both_haps"]]).sort()


a = pybedtools.BedTool.from_dataframe(vert_rt_windows)
vert = a.merge(d=250001).filter(lambda x: len(x) > 250000).saveas().to_dataframe().reset_index()
vert = vert.drop("index",axis="columns")
vert.columns = ["chrom","start","stop"]
vert = vert.drop_duplicates()
vert["chrom"] = vert["chrom"].astype(str)
vert_rt_windows_bed = pybedtools.BedTool.from_dataframe(vert).sort()

### regular genes
df_eb_genes = pd.read_csv("gm.as.gene.counts.txt",sep="\t")
aei_std_dev_df = df_eb_genes.groupby(["chrom","start","stop","name"])["aei"].std().reset_index()
aei_std_dev_df.columns = ["chrom","start","stop","name","aei_std_dev"]

aei_std_dev_df  = aei_std_dev_df.dropna(subset=["aei_std_dev"])
aei_std_dev_df = aei_std_dev_df[aei_std_dev_df["aei_std_dev"]!=0]
aei_std_dev_df = aei_std_dev_df[aei_std_dev_df["chrom"]!="chrX"]
### test test
# aei_std_dev_df_high = aei_std_dev_df[aei_std_dev_df["aei_std_dev"]>=0.25]

sns.kdeplot(aei_std_dev_df["aei_std_dev"],clip=(0,1),linewidth=4)
plt.show()
plt.close()
# exit()



aei_std_dev_bed = pybedtools.BedTool.from_dataframe(aei_std_dev_df).sort()
gene_aei_intersect_all_rt_windows = aei_std_dev_bed.intersect(all_rt_windows_bed_unmerged,wa=True,wb=True)
# print(gene_aei_intersect_rt_vert)
gene_aei_intersect_all_rt_windows_df = gene_aei_intersect_all_rt_windows.to_dataframe(disable_auto_names=True, header=None)
gene_aei_intersect_all_rt_windows_df.columns = ["chrom","start","stop","name","aei_std_dev",
                                        "chrom_rt","start_rt","stop_rt","gm_std_dev_both_haps"]


### try a fisher exact test
aei_high = gene_aei_intersect_all_rt_windows_df['aei_std_dev'] >= np.percentile(gene_aei_intersect_all_rt_windows_df['aei_std_dev'], 97)  # top 10%
rt_high = gene_aei_intersect_all_rt_windows_df['gm_std_dev_both_haps'] >= np.percentile(gene_aei_intersect_all_rt_windows_df['gm_std_dev_both_haps'], 97)
a = np.sum(aei_high & rt_high)     # both high
print("num aei high rt high",a)
b = np.sum(aei_high & ~rt_high)    # X high, Y low
print("num aei high rt low",b)
c = np.sum(~aei_high & rt_high)    # X low, Y high
print("num aei low rt high",c)
d = np.sum(~aei_high & ~rt_high)   # both low
print("num aei low rt low",d)
## 97th percentile seems highest
oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
print("odds ratio",oddsratio)
print("fisher pval",pval)

gene_aei_intersect_all_rt_windows_df["variable_aei"] = gene_aei_intersect_all_rt_windows_df['aei_std_dev'] >= np.percentile(gene_aei_intersect_all_rt_windows_df['aei_std_dev'], 97)
gene_aei_intersect_all_rt_windows_df["variable_rt"] = gene_aei_intersect_all_rt_windows_df['gm_std_dev_both_haps'] >= np.percentile(gene_aei_intersect_all_rt_windows_df['gm_std_dev_both_haps'], 97)


high_aei_high_rt = gene_aei_intersect_all_rt_windows_df[(gene_aei_intersect_all_rt_windows_df["variable_rt"]==True) & (gene_aei_intersect_all_rt_windows_df["variable_aei"]==True)]

print(high_aei_high_rt.sort_values(["aei_std_dev"],ascending=False).drop_duplicates("name"))
exit()
# print(gene_aei_intersect_all_rt_windows_df) 

# pearson
corr, p_value = pearsonr(gene_aei_intersect_all_rt_windows_df["aei_std_dev"].astype(float), 
    gene_aei_intersect_all_rt_windows_df["gm_std_dev_both_haps"].astype(float))
print("pearson Correlation:", corr)
print("P-value:", p_value)

# spearman
# corr, p_value = spearmanr(gene_aei_intersect_rt_vert_df["aei_std_dev"].astype(float), 
#     gene_aei_intersect_rt_vert_df["gm_std_dev_both_haps"].astype(float))
# print("Spearman Correlation:", corr)
# print("P-value:", p_value)


### try natty log
# corr, p_value = pearsonr(np.log(gene_aei_intersect_rt_vert_df["aei_std_dev"].astype(float)), 
#     np.log(gene_aei_intersect_rt_vert_df["gm_std_dev_both_haps"].astype(float)))
# print("nautral log pearson Correlation:", corr)
# print("P-value:", p_value)

rt_vert_intersect_aei_std_dev = vert_rt_windows_bed.intersect(aei_std_dev_bed,wa=True,wb=True)
rt_vert_intersect_aei_std_dev_df = rt_vert_intersect_aei_std_dev.to_dataframe(disable_auto_names=True, header=None)
rt_vert_intersect_aei_std_dev_df.columns = ["chrom","start","stop","name",
                                        "chrom_rt","start_rt","stop_rt","aei_std_dev"]

non_vert_genes = aei_std_dev_bed.subtract(vert_rt_windows_bed)
non_vert_genes_df = non_vert_genes.to_dataframe(disable_auto_names=True, header=None)
non_vert_genes_df.columns=["chrom","start","stop","name","aei_std_dev"]
# print(non_vert_genes)


fig,ax=plt.subplots(figsize=(6,4))
sns.kdeplot(rt_vert_intersect_aei_std_dev_df["aei_std_dev"],linestyle="--",c="black",linewidth=2,clip=(0,.6))
sns.kdeplot(non_vert_genes_df["aei_std_dev"],c="black",linewidth=2,clip=(0,.6))
plt.show()

# Mann–Whitney U test (two-sided by default)
stat, p = mannwhitneyu(rt_vert_intersect_aei_std_dev_df["aei_std_dev"], 
                        non_vert_genes_df["aei_std_dev"], alternative='two-sided')

print(f"U statistic: {stat}")
print(f"P-value: {p}")




all_genes_bed = pybedtools.BedTool.from_dataframe(df_eb_genes[df_eb_genes["chrom"]!="chrX"]).sort().merge()
sig_genes_bed = pybedtools.BedTool.from_dataframe(df_eb_genes[(df_eb_genes["fdr_reject"]==True) & 
                                    (df_eb_genes["chrom"]!="chrX") &
                                    (df_eb_genes["70_percent_aei"]==True) &
                                    (df_eb_genes["fdr_pval"]<=0.001)].loc[:,["chrom","start","stop"]].drop_duplicates()).sort()#.merge()
all_genes_bed.to_dataframe(disable_auto_names=True, header=None).to_csv("all_genes_gm.bed",sep="\t",header=None,index=None)


# print(df_gm_genes[(df_gm_genes["fdr_reject"]==True) & 
#     (df_gm_genes["chrom"]!="chrX") & 
#     (df_gm_genes["70_percent_aei"]==True)].loc[:,["chrom","start","stop"]].drop_duplicates())
###
real_result = len(vert_rt_windows_bed.intersect(sig_genes_bed,wa=True,wb=True))

vert_rt_windows_bed.intersect(sig_genes_bed,wa=True,wb=True).to_dataframe(disable_auto_names=True, header=None).to_csv("gm.vert.rt.windows.bed",sep="\t",index=None)
sim_results=[]
print(len(sig_genes_bed))
print(len(vert_rt_windows_bed))
print("num real intersections", real_result)

# shuffle vert and genes
# for i in range(5000):
#     shuffle_vert = vert_rt_windows_bed.shuffle(g="../hg38.fa.fai", incl="all_rt_windows.bed",chromFirst=True,noOverlapping=True)
#     shuffle_genes = sig_genes_bed.shuffle(g="../hg38.fa.fai", incl="all_genes.bed",chromFirst=True,noOverlapping=True)
#     sim_results += [len(shuffle_vert.intersect(shuffle_genes,wa=True,wb=True))]

## shuffle genes only
for i in range(10000):
    shuffle_vert = vert_rt_windows_bed#.shuffle(g="../hg38.fa.fai", incl="all_rt_windows.bed",chromFirst=True)
    shuffle_genes = sig_genes_bed.shuffle(g="../hg38.fa.fai", incl="all_genes_gm.bed",chromFirst=True)
    sim_results += [len(shuffle_vert.intersect(shuffle_genes,wa=True,wb=True))]

# ## shuffle vert only
# for i in range(10000):
#     shuffle_vert = vert_rt_windows_bed.shuffle(g="../hg38.fa.fai", incl="all_rt_windows.bed",chromFirst=True)
#     shuffle_genes = sig_genes_bed#.shuffle(g="../hg38.fa.fai", incl="all_genes.bed",chromFirst=True)
#     sim_results += [len(shuffle_vert.intersect(shuffle_genes,wa=True,wb=True))]

# f,ax=plt.subplots(figsize=(2,2),dpi=300)
# sns.kdeplot(sim_results)
# ax.axvline(x=real_result,lw=0.5,linestyle="--",c="black")

# plt.show()
permuted_overlaps = np.array(sim_results)
p_value = np.sum(permuted_overlaps >= real_result) / len(permuted_overlaps)
print("10,000 simulations empirical p_value",p_value)


# Example values
# observed = 25
# null_overlaps = np.array([15, 18, 20, 17, 23, 22, 21, 19, 24, 20, 22, 18, 16, 25, 19])  # your null distribution

# Compute mean and std of null
mu = np.mean(permuted_overlaps)
sigma = np.std(permuted_overlaps, ddof=1)
print("mean and std dev of simulated intersections")
print(mu, sigma)

# Compute z-score
z = (real_result - mu) / sigma

# One-sided p-value (test for enrichment)
p_value = 1 - norm.cdf(z)


print("parametric z score", z)
print("parametric pval",p_value)
exit()


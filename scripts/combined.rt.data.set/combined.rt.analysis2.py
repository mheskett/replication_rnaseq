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
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from scipy.stats import ttest_ind
import glob

def add_binom_pval(df):
    """
    example:
    >>> scipy.stats.binomtest(5,10,p=0.5)
    BinomTestResult(k=5, n=10, alternative='two-sided', statistic=0.5, pvalue=1.0)
    """
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binomtest(row["paternal_counts"],
                                row["paternal_counts"]+row["maternal_counts"],
                            p=0.5,
                            alternative="two-sided").pvalue, # v slow for some reason 
                            axis=1)
    results = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = results[1]
    df["fdr_reject"] = results[0]


def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["paternal_counts"] >= x["maternal_counts"]:
        return x["paternal_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["maternal_counts"]  / x["total_reads"] + 0.5
    return

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
def get_arms(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes)):
        # should be (p end, q end)
        arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
        cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
    return arm_dict

def get_arms_nochr(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes_nochr)):
        # should be (p end, q end)
        arm_dict[chromosomes_nochr[i]] = (cytoband_nochr[(cytoband_nochr["chrom"]==chromosomes_nochr[i]) & (cytoband_nochr["arm"].str.contains("p"))]["stop"].max(),
        cytoband_nochr[(cytoband_nochr["chrom"]==chromosomes_nochr[i]) & (cytoband_nochr["arm"].str.contains("q"))]["stop"].max())
    return arm_dict


def sum_region_length(df):
    diffs = df["stop"] - df["start"]
    return diffs.sum()

### REMOVIONG 15 for BAD DATA!!!!
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

chromosomes_nochr = ["1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X"]

arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])

cytoband_nochr = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)
arm_dict_nochr=get_arms_nochr(cytoband_nochr)

chromosome_length = {"chr1":249250621,
"chr2":243199373,
"chr3":198022430,
"chr4":191154276,
"chr5":180915260,
"chr6":171115067,
"chr7":159138663,
"chr8":146364022,
"chr9":141213431,
"chr10":135534747,
"chr11":135006516,
"chr12":133851895,
"chr13":115169878,
"chr14":107349540,
"chr15":102531392,
"chr16":90354753,
"chr17":81195210,
"chr18":78077248,
"chr19":59128983,
"chr20":63025520,
"chr21":48129895,
"chr22":51304566,
"chrX":155270560}

chromosome_length_nochr = {"1":249250621,
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


chromosome_length_mouse = {"chr1":    195471971,
"chr2":    182113224,
"chrX":    171031299,
"chr3":    160039680,
"chr4":    156508116,
"chr5":    151834684,
"chr6":    149736546,
"chr7":   145441459,
"chr10":   130694993,
"chr8":    129401213,
"chr14":   124902244,
"chr9":   124595110,
"chr11":   122082543,
"chr13":   120421639,
"chr12":   120129022,
"chr15":   104043685,
"chr16":   98207768,
"chr17":   94987271,
"chrY" :   91744698,
"chr18":   90702639,
"chr19":   61431566
}

chromosomes_mouse = {"chr1",
"chr2",
"chrX",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr10",
"chr8",
"chr14",
"chr9",
"chr11",
"chr13",
"chr12",
"chr15",
"chr16",
"chr17",
"chrY",
"chr18",
"chr19"}


##############################
# ACP SAMPLES

##############################

#### ACP6 ACP6 ACP6 ACP6
#### 
#### coding genes coding genes coding genes coding genes coding genes coding genes coding genes coding genes

dfs=[]
acp6_rna = glob.glob("../rna.dec2024/acp6*haplotype.resolved.counts.comprehensive.gene.counts.bed")
for f in acp6_rna:
    samp=os.path.basename(f).split(".")[0].split("Aligned")[0]
    tmp = pd.read_csv(f,sep="\t",names=["chrom","start","stop","gene","score","gene_strand",
                             "tx_start","tx_stop","paternal_counts","maternal_counts","reads_strand"])
    tmp = tmp.set_index(["chrom","start","stop"])
    tmp["total_reads"] = tmp["paternal_counts"] + tmp["maternal_counts"]
    tmp=tmp[tmp["total_reads"]>0]
    tmp["sample"] = samp
    tmp["aei"] =  tmp.apply(helper_func, axis = 1)
    dfs += [tmp]
df_acp6rna = pd.concat(dfs).sort_values(["chrom","start"])
add_binom_pval(df_acp6rna)
df_acp6rna = df_acp6rna.reset_index()
####
####
####
#### Replication timing Replication timing Replication timing Replication timing Replication timing 
acp6_samples = glob.glob("./acp6_113024/acp6*haplotype*counts.windows.bed")
dfs=[]

### MALE SAMPLE NOT FEMALE
# now we just have an early and late hap1 and hap2 within columns
for x in acp6_samples:
    samp = os.path.basename(x).split(".")[0]
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", samp+"_p_counts", samp+"_m_counts"])
    tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 

    tmp=tmp.set_index(["chrom","start","stop"])
    tmp[samp+"_cpm_p_counts"] = (tmp[samp+"_p_counts"] / tmp[samp+"_p_counts"].sum(axis="rows"))*10**6
    tmp[samp+"_cpm_m_counts"] = (tmp[samp+"_m_counts"] / tmp[samp+"_m_counts"].sum(axis="rows"))*10**6
    dfs+=[tmp]
df_acp6 = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
samples=list(set(["_".join(os.path.basename(x).split(".")[0].split("_")[0:2]) for x in acp6_samples]))
print(samples)
for sample in samples:
    df_acp6[sample+"_paternal_logrt"] = np.log2((df_acp6[sample+"_early_cpm_p_counts"]+1) / (df_acp6[sample+"_late_cpm_p_counts"]+1 ))
    df_acp6[sample+"_maternal_logrt"] = np.log2((df_acp6[sample+"_early_cpm_m_counts"]+1) / (df_acp6[sample+"_late_cpm_m_counts"]+1 ))


#### FILTER. FILTER. Drop rows that dont have 100% sample representation. 
#### drop bad regions of the genome that should have already been dropped
print(len(df_acp6))
print("dropped",len(df_acp6.dropna(how="any",axis=0)))
df_acp6=df_acp6.dropna(how="any",axis=0)
######

df_acp6_qn = quantile_normalize(df_acp6.filter(regex="logrt"))
df_acp6_qn["acp6_std_dev_both_haps"] = df_acp6_qn.filter(like="acp6",axis=1).std(axis="columns")
df_acp6_qn = df_acp6_qn.reset_index()
mean_std_dev = df_acp6_qn[df_acp6_qn["chrom"]!="X"]["acp6_std_dev_both_haps"].mean()
std_dev_dev = df_acp6_qn[df_acp6_qn["chrom"]!="X"]["acp6_std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_acp6_qn = df_acp6_qn.set_index(["chrom","start","stop"])
df_acp6_qn["acp6_vert"] =  df_acp6_qn.apply(lambda x:True if x["acp6_std_dev_both_haps"]>=threshold else False,axis=1)

df_acp6_qn = df_acp6_qn.reset_index()
df_acp6_qn["arm"] = df_acp6_qn.apply(lambda x: "q" if (x["stop"] > arm_dict_nochr[x["chrom"]][0]) & (x["stop"] <= arm_dict_nochr[x["chrom"]][1]) else "p", axis=1)


print("acp6 225",df_acp6_qn)
##### make distribution of repliseq std dev. male sample
f,ax=plt.subplots(figsize=(2,2),dpi=300)
sns.kdeplot(df_acp6_qn["acp6_std_dev_both_haps"],clip=(0,20),linewidth=2)
ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 3 * std_dev_dev,lw=0.5,linestyle="--",c="red")


# plt.show()
plt.close()

color_dict_acp6 = {
"acp6_c1":"red", 
"acp6_c2":"cyan",  
"acp6_c5":"yellow",  
"acp6_c6":"green"}

color_dict_acp6_rna = {
"acp6_c1_rna":"red", 
"acp6_c2_rna":"cyan",  
"acp6_c5_rna":"yellow",  
"acp6_c6_rna":"green"}
#####
df_acp6_qn = df_acp6_qn.sort_values(["chrom","start"])

for chrom in chromosomes_nochr:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(2,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_acp6_qn[(df_acp6_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_paternal_logrt").reset_index()
        maternal = df_acp6_qn[(df_acp6_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_maternal_logrt").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")
        # paternal["arm"]=paternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
        # maternal["arm"]=maternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

        for j in ["p","q"]:

            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][samples[i]+"_paternal_logrt"],
                    c=color_dict_acp6[samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][samples[i]+"_maternal_logrt"],
                    c=color_dict_acp6[samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_acp6_qn[(df_acp6_qn["chrom"]==chrom) & 
                            (df_acp6_qn["acp6_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_acp6_qn[(df_acp6_qn["chrom"]==chrom) & (df_acp6_qn["arm"]==k)]["start"],
                        df_acp6_qn[(df_acp6_qn["chrom"]==chrom) & (df_acp6_qn["arm"]==k)]["acp6_std_dev_both_haps"],
                        c="black",
                        lw=0.7)

        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[1].set_xticks([]) 



    plt.savefig("acp6.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()

print(df_acp6_qn)
df_acp6_qn.to_csv("acp6.rt.txt",sep="\t",index=False)
df_acp6_qn.to_csv("acp6.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)


####### zooms zooms zooms zooms zooms zooms 
#######
regions = df_acp6_qn[df_acp6_qn["acp6_vert"]==True]
print("acp6 regions", regions)
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=500001).to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)

print("acp6 region length sum", sum_region_length(regions))

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(2,4),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_paternal_logrt").reset_index()
        maternal = df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_maternal_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-2000000) & (paternal["stop"]<=row["stop"]+2000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-2000000) & (maternal["stop"]<=row["stop"]+2000000)]

        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:
            ## unsmoothed
            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][samples[i]+"_paternal_logrt"],
                    c=color_dict_acp6[samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][samples[i]+"_maternal_logrt"],
                    c=color_dict_acp6[samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

            # #smoothed 
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([row["start"]-2000000,row["stop"]+2000000])
        ax[0].set_xticks([])

        ##### coding genes
        #######
        ## currently pklotting the total and the plus/minus and the antisense. so pick one to plot....
        ax2 = ax[0].twinx()
        rna_tmp = df_acp6rna[(df_acp6rna["chrom"]==row["chrom"]) & 
                            (df_acp6rna["start"]>=row["start"]-2000000) & 
                            (df_acp6rna["stop"]<=row["stop"]+2000000) & 
                            (df_acp6rna["total_reads"]>=10) & 
                            (df_acp6rna["reads_strand"].isin(["plus","minus"]))]

        for index2, row2 in rna_tmp.iterrows():
            rect=Rectangle((row2["start"], row2["aei"]-.0125), 
                            width=row2["stop"]-row2["start"], height=0.025,
                            facecolor=color_dict_acp6_rna[row2["sample"]], edgecolor="black",
                            fill=True,lw=.2)

            ax2.add_patch(rect)
        ax2.set_ylim([-0.52,0.52])
        ax2.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1, 0, .1, .2, .3, .4, 0.5])
        ######
        ######
        #### highlighting allele variant regions
        for index3,row3 in df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"]) & 
                            (df_acp6_qn["acp6_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"]) & (df_acp6_qn["arm"]==k)]["start"],
                        df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"]) & (df_acp6_qn["arm"]==k)]["acp6_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-2000000,row["stop"]+2000000])
        ax[1].set_xticks(np.linspace(row["start"]-2000000,row["stop"]+2000000,5)) 
    plt.savefig("acp6.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()


######
######
######
######
### do the STD DEV work before combining

##### ACP 7 ACP7 ACP7 ACP7
######
##### ACP7 RNA ACP7 RNA ACP7 RNA ACP7 RNA ACP7 RNA ACP7 RNA 
dfs=[]
acp7_rna = glob.glob("../rna.dec2024/acp7*haplotype.resolved.counts.comprehensive.gene.counts.bed")
for f in acp7_rna:
    samp=os.path.basename(f).split(".")[0].split("Aligned")[0]
    tmp = pd.read_csv(f,sep="\t",names=["chrom","start","stop","gene","score","gene_strand",
                             "tx_start","tx_stop","paternal_counts","maternal_counts","reads_strand"])
    tmp = tmp.set_index(["chrom","start","stop"])
    tmp["total_reads"] = tmp["paternal_counts"] + tmp["maternal_counts"]
    tmp=tmp[tmp["total_reads"]>0]
    tmp["sample"] = samp
    tmp["aei"] =  tmp.apply(helper_func, axis = 1)
    dfs += [tmp]
df_acp7rna = pd.concat(dfs).sort_values(["chrom","start"])
add_binom_pval(df_acp7rna)
df_acp7rna = df_acp7rna.reset_index()

#####
####
####
acp7_samples = glob.glob("./acp7_120324/acp7*haplotype*counts.windows.bed")

dfs={}
dfs=[]
# now we just have an early and late hap1 and hap2 within columns
for x in acp7_samples:
    samp = os.path.basename(x).split(".")[0]
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", samp+"_p_counts", samp+"_m_counts"])
    # tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 

    tmp=tmp.set_index(["chrom","start","stop"])
    tmp[samp+"_cpm_p_counts"] = (tmp[samp+"_p_counts"] / tmp[samp+"_p_counts"].sum(axis="rows"))*10**6
    tmp[samp+"_cpm_m_counts"] = (tmp[samp+"_m_counts"] / tmp[samp+"_m_counts"].sum(axis="rows"))*10**6
    dfs+=[tmp]
df_acp7 = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
samples=list(set(["_".join(os.path.basename(x).split(".")[0].split("_")[0:2]) for x in acp7_samples]))
print(samples)
for sample in samples:
    df_acp7[sample+"_paternal_logrt"] = np.log2((df_acp7[sample+"_early_cpm_p_counts"]+1) / (df_acp7[sample+"_late_cpm_p_counts"]+1 ))
    df_acp7[sample+"_maternal_logrt"] = np.log2((df_acp7[sample+"_early_cpm_m_counts"]+1) / (df_acp7[sample+"_late_cpm_m_counts"]+1 ))



df_acp7=df_acp7.dropna(how="any",axis=0)

df_acp7_qn = quantile_normalize(df_acp7.filter(regex="logrt"))
df_acp7_qn["acp7_std_dev_both_haps"] = df_acp7_qn.filter(like="acp7",axis=1).std(axis="columns")
df_acp7_qn = df_acp7_qn.reset_index()
mean_std_dev = df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"].mean()
std_dev_dev = df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_acp7_qn = df_acp7_qn.set_index(["chrom","start","stop"])
df_acp7_qn["acp7_vert"] =  df_acp7_qn.apply(lambda x:True if x["acp7_std_dev_both_haps"]>=threshold else False,axis=1)
df_acp7_qn = df_acp7_qn.reset_index()
df_acp7_qn["arm"] = df_acp7_qn.apply(lambda x: "q" if (x["stop"] > arm_dict_nochr[x["chrom"]][0]) & (x["stop"] <= arm_dict_nochr[x["chrom"]][1]) else "p", axis=1)


##### make distribution of repliseq std dev
f,ax=plt.subplots(figsize=(2,1),dpi=200)
sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=2)
sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]=="X"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=2,c="red")
# ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
# ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# plt.show()
plt.close()
#### normal transform the data, then compare to a perfect normal
f,ax=plt.subplots(figsize=(2,2),dpi=200)
sns.kdeplot((scipy.stats.zscore(df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"])))
sns.kdeplot(np.random.normal(0, 1, 100000))

# plt.show()
print("skew",scipy.stats.skewtest(a=df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"]))
print("kurtosis",scipy.stats.kurtosistest(a=df_acp7_qn[df_acp7_qn["chrom"]!="X"]["acp7_std_dev_both_haps"]))

### skew test

## kurtosis test



color_dict_acp7 = {
"acp7_c2":"red", 
"acp7_c4":"green",  
"acp7_c5":"blue"}

color_dict_acp7_rna = {
"acp7_c1_rna":"yellow", 
"acp7_c2_rna":"red", 
"acp7_c4_rna":"green",  
"acp7_c5_rna":"blue"}
#####
df_acp7_qn = df_acp7_qn.sort_values(["chrom","start"])

for chrom in chromosomes_nochr:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(2,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_acp7_qn[(df_acp7_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_paternal_logrt").reset_index()
        maternal = df_acp7_qn[(df_acp7_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_maternal_logrt").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")
        # paternal["arm"]=paternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
        # maternal["arm"]=maternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

        for j in ["p","q"]:

            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][samples[i]+"_paternal_logrt"],
                    c=color_dict_acp7[samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][samples[i]+"_maternal_logrt"],
                    c=color_dict_acp7[samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_acp7_qn[(df_acp7_qn["chrom"]==chrom) & 
                            (df_acp7_qn["acp7_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_acp7_qn[(df_acp7_qn["chrom"]==chrom) & (df_acp7_qn["arm"]==k)]["start"],
                        df_acp7_qn[(df_acp7_qn["chrom"]==chrom) & (df_acp7_qn["arm"]==k)]["acp7_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[1].set_xticks([]) 
    plt.savefig("acp7.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()

df_acp7_qn.to_csv("acp7.rt.txt",sep="\t",index=False)
df_acp7_qn.to_csv("acp7.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)
### zooms zooms zooms zooms zooms zooms 
####
####
####
####
print("acp7",df_acp7_qn[(df_acp7_qn["acp7_vert"]==True) & (df_acp7_qn["chrom"]!="X")])
regions = df_acp7_qn[(df_acp7_qn["acp7_vert"]==True) & (df_acp7_qn["chrom"]!="X")]
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=500001).to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)
print(regions)

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(2,4),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_paternal_logrt").reset_index()
        maternal = df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_maternal_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-3000000) & (paternal["stop"]<=row["stop"]+3000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-3000000) & (maternal["stop"]<=row["stop"]+3000000)]

        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:
            ## unsmoothed
            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][samples[i]+"_paternal_logrt"],
                    c=color_dict_acp7[samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][samples[i]+"_maternal_logrt"],
                    c=color_dict_acp7[samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2


            # #smoothed 
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ##### coding genes
        #######
        ## currently pklotting the total and the plus/minus and the antisense. so pick one to plot....
        ax2 = ax[0].twinx()
        rna_tmp = df_acp7rna[(df_acp7rna["chrom"]==row["chrom"]) & 
                            (df_acp7rna["start"]>=row["start"]-2000000) & 
                            (df_acp7rna["stop"]<=row["stop"]+2000000) & 
                            (df_acp7rna["total_reads"]>=10) & 
                            (df_acp7rna["reads_strand"].isin(["plus","minus"]))]

        for index2, row2 in rna_tmp.iterrows():
            rect=Rectangle((row2["start"], row2["aei"]-.0125), 
                            width=row2["stop"]-row2["start"], height=0.025,
                            facecolor=color_dict_acp7_rna[row2["sample"]], edgecolor="black",
                            fill=True,lw=.2)

            ax2.add_patch(rect)
        ax2.set_ylim([-0.52,0.52])
        ax2.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1, 0, .1, .2, .3, .4, 0.5])


        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"]) & 
                            (df_acp7_qn["acp7_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"]) & (df_acp7_qn["arm"]==k)]["start"],
                        df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"]) & (df_acp7_qn["arm"]==k)]["acp7_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[1].set_xticks(np.linspace(row["start"]-3000000,row["stop"]+3000000,5)) 
    plt.savefig("acp7.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()

##############################
###############################
###############################
## Mouse data

###############################
bad_samples = [
"GSM3756321_pre-b_clone3.1.CAST.bedGraph", 
"GSM3756323_pre-b_clone3.3.CAST.bedGraph",      
"GSM3756322_pre-b_clone3.2.CAST.bedGraph",
"GSM3756321_pre-b_clone3.1.C57BL.bedGraph",      
"GSM3756323_pre-b_clone3.3.C57BL.bedGraph",      
"GSM3756322_pre-b_clone3.2.C57BL.bedGraph"]
 
# this has mulitple clones AND replicates 123

all_files_repli = glob.glob("GSM*")
all_files_repli = [x for x in all_files_repli if x not in bad_samples]
all_files_repli = [x for x in all_files_repli if "pure" not in x]
cast = [x for x in all_files_repli if "CAST" in x]
c57 = [x for x in all_files_repli if "C57" in x]
# pure = [x for x in all_files_repli if "pure" in x]

repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",header=0,
                        names= ["chrom","start","stop",all_files_repli[i]],
                        dtype = {"chrom":str,"start":int,"stop":int, all_files_repli[i]:float} )
    # df_repli=df_repli[df_repli["chrom"]!="chrX"]
    # df_repli=df_repli[df_repli["chrom"]!="chrY"]
    repli_li += [df_repli.set_index(["chrom","start","stop"])]

df_mouse = pd.concat(repli_li,axis=1)#.reset_index() 

df_mouse = df_mouse.reset_index()
df_mouse["mouse_std_dev_both_haps"] = df_mouse.filter(like="GSM",axis=1).std(axis="columns")
mean_std_dev = df_mouse[df_mouse["chrom"]!="chrX"]["mouse_std_dev_both_haps"].mean()
std_dev_dev = df_mouse[df_mouse["chrom"]!="chrX"]["mouse_std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_mouse = df_mouse.sort_values(["chrom","start","stop"])
df_mouse = df_mouse.set_index(["chrom","start","stop"])
df_mouse["mouse_vert"] =  df_mouse.apply(lambda x:True if x["mouse_std_dev_both_haps"]>=threshold else False,axis=1)

##### make distribution of repliseq std dev
df_mouse = df_mouse.reset_index()
f,ax=plt.subplots(figsize=(2,2),dpi=200)
sns.kdeplot(df_mouse[df_mouse["chrom"]!="chrX"]["mouse_std_dev_both_haps"],clip=(0,20),linewidth=2)
sns.kdeplot(df_mouse[df_mouse["chrom"]=="chrX"]["mouse_std_dev_both_haps"],clip=(0,20),linewidth=2,c="red")
# plt.show()

print(df_mouse)
print(all_files_repli)
samples = [".".join(x.split(".")[0:2]) for x in all_files_repli]
print(samples)
# color_dict_mouse={
#  'GSM3756320_pre-b_clone2.3', 
#  'GSM3756326_pre-b_clone8.3', 
#  'GSM3756326_pre-b_clone8.3', 
#  'GSM3756324_pre-b_clone8.1', 
#  'GSM3756316_pre-b_clone_e5.2', 
#  'GSM3756319_pre-b_clone2.2', 
#  'GSM3756317_pre-b_clone_e5.3', 
#  'GSM3756316_pre-b_clone_e5.2', 
#  'GSM3756315_pre-b_clone_e5.1',
# 'GSM3756318_pre-b_clone2.1', 
# 'GSM3756319_pre-b_clone2.2',
# 'GSM3756315_pre-b_clone_e5.1',
#     'GSM3756318_pre-b_clone2.1', 
#     'GSM3756325_pre-b_clone8.2', 
#     'GSM3756317_pre-b_clone_e5.3',
#      'GSM3756320_pre-b_clone2.3', 
#      'GSM3756325_pre-b_clone8.2', 
#      'GSM3756324_pre-b_clone8.1'  

# }
#### mouse plots

for chrom in chromosomes_mouse:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(2,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        cast = df_mouse[(df_mouse["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=samples[i]+".CAST.bedGraph").reset_index()
        c57 = df_mouse[(df_mouse["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=samples[i]+".C57BL.bedGraph").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        ax[0].plot(cast["start"],
                cast[samples[i]+".CAST.bedGraph"],
                c="red",lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(c57["start"],
                c57[samples[i]+".C57BL.bedGraph"],
                c="blue",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-2.2,2.2])
        ax[0].set_yticks([-2,-1,0,1,2])
        ax[0].set_xlim([0,chromosome_length_mouse[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_mouse[(df_mouse["chrom"]==chrom) & 
                            (df_mouse["mouse_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)


        ## std dev 
        ax[1].plot(df_mouse[(df_mouse["chrom"]==chrom) ]["start"],
                    df_mouse[(df_mouse["chrom"]==chrom)]["mouse_std_dev_both_haps"],
                    c="black",
                    lw=0.7)

        ax[1].set_ylim([0,1])
        ax[1].set_yticks([0,1])
        ax[1].set_xlim([0,chromosome_length_mouse[chrom]])
        ax[1].set_xticks([]) 
    plt.savefig("mouse.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()


### zooms zooms zooms zooms zooms zooms 
####
####
####
####
regions = df_mouse[(df_mouse["mouse_vert"]==True) & (df_mouse["chrom"]!="chrX")]
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=500001).to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)
print(regions)

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(2,4),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        cast = df_mouse[(df_mouse["chrom"]==row["chrom"])].set_index(["chrom","start","stop"]).filter(like=samples[i]+".CAST.bedGraph").reset_index()
        c57 = df_mouse[(df_mouse["chrom"]==row["chrom"])].set_index(["chrom","start","stop"]).filter(like=samples[i]+".C57BL.bedGraph").reset_index()

        cast = cast[(cast["chrom"]==row["chrom"]) & (cast["start"]>=row["start"]-3000000) & (cast["stop"]<=row["stop"]+3000000)]
        c57 = c57[(c57["chrom"]==row["chrom"]) & (c57["start"]>=row["start"]-3000000) & (c57["stop"]<=row["stop"]+3000000)]

        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        ## unsmoothed
        ax[0].plot(cast["start"],
                cast[samples[i]+".CAST.bedGraph"],
                c="red",lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(c57["start"],
                c57[samples[i]+".C57BL.bedGraph"],
                c="blue",lw=0.5) ## -- line style is haplotype 2


            # #smoothed 
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ##### coding genes
        #######
        ## currently pklotting the total and the plus/minus and the antisense. so pick one to plot....
        ax[0].set_ylim([-2.2,2.2])
        ax[0].set_yticks([2,-1,0,1,2])
        ax[0].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_mouse[(df_mouse["chrom"]==row["chrom"]) & 
                            (df_mouse["mouse_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)



        ax[1].plot(df_mouse[(df_mouse["chrom"]==row["chrom"])]["start"],
                    df_mouse[(df_mouse["chrom"]==row["chrom"])]["mouse_std_dev_both_haps"],
                    c="black",
                    lw=0.7)
        ax[1].set_ylim([0,1])
        ax[1].set_yticks([0,1])
        ax[1].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[1].set_xticks(np.linspace(row["start"]-3000000,row["stop"]+3000000,5)) 
    plt.savefig("mouse.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()



df_mouse.to_csv("mouse.rt.txt",sep="\t",index=False)
df_mouse.to_csv("mouse.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)





## with the output of the mouse file youll have to use the UCSC tool to 
## get hg19 coordinates
####################################
####################################
### eb samples eb samples eb samples eb samples eb samples eb samples eb samples eb samples
#########################
#######################
#######################


eb_samples = glob.glob("eb*rt*windows.bed")
dfs=[]
for x in eb_samples:
    x_samp="_".join(x.split("_")[0:4])
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", x_samp+"_hap1_counts", x_samp+"_hap2_counts"])
    # tmp=tmp[tmp["chrom"]!="X"]###removing chrom x 
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 

    tmp=tmp.set_index(["chrom","start","stop"])
    tmp[x_samp+"_cpm_hap1_counts"] = (tmp[x_samp+"_hap1_counts"] ) / tmp[x_samp+"_hap1_counts"].sum(axis="rows")*10**6
    tmp[x_samp+"_cpm_hap2_counts"] = (tmp[x_samp+"_hap2_counts"] ) / tmp[x_samp+"_hap2_counts"].sum(axis="rows")*10**6
    dfs+=[tmp]

df_eb = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
eb_samples = ["_".join(x.split("_")[0:3]) for x in eb_samples]
for sample in eb_samples:
    df_eb[sample+"_hap1_logrt"] = np.log2((df_eb[sample+"_early_cpm_hap1_counts"]+1) / (df_eb[sample+"_late_cpm_hap1_counts"]+1 ))
    df_eb[sample+"_hap2_logrt"] = np.log2((df_eb[sample+"_early_cpm_hap2_counts"]+1) / (df_eb[sample+"_late_cpm_hap2_counts"]+1 ))

print(df_eb)

df_eb=df_eb.dropna(how="any",axis=0)

df_eb_qn = quantile_normalize(df_eb.filter(regex="logrt"))
# df_eb_qn = df_eb_qn.reset_index()
### do the STD DEV work before combining
print("asdf",df_eb_qn)

df_eb_qn["eb_std_dev_both_haps"] = df_eb_qn.filter(like="eb",axis=1).std(axis="columns")
df_eb_qn = df_eb_qn.reset_index()
print("ffff",df_eb_qn)

mean_std_dev = df_eb_qn[df_eb_qn["chrom"]!="X"]["eb_std_dev_both_haps"].mean()
std_dev_dev = df_eb_qn[df_eb_qn["chrom"]!="X"]["eb_std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_eb_qn = df_eb_qn.set_index(["chrom","start","stop"])
print("dddd",df_eb_qn)

df_eb_qn["eb_vert"] =  df_eb_qn.apply(lambda x:True if x["eb_std_dev_both_haps"]>=threshold else False,axis=1)
df_eb_qn = df_eb_qn.sort_values(["chrom","start"])

df_eb_qn = df_eb_qn.reset_index()
df_eb_qn["arm"] = df_eb_qn.apply(lambda x: "q" if (x["stop"] > arm_dict_nochr[x["chrom"]][0]) & (x["stop"] <= arm_dict_nochr[x["chrom"]][1]) else "p", axis=1)


f,ax=plt.subplots(figsize=(2,1),dpi=200)
sns.kdeplot(df_eb_qn[df_eb_qn["chrom"]!="X"]["eb_std_dev_both_haps"],clip=(0,20),linewidth=2)
sns.kdeplot(df_eb_qn[df_eb_qn["chrom"]=="X"]["eb_std_dev_both_haps"],clip=(0,20),linewidth=2,c="red")
# ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
# ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# plt.show()

color_dict_eb = {'eb3_2_clone15':"red", 
'eb3_2_clone4':"cyan", 
'eb3_2_clone13':"yellow", 
'eb3_2_clone10':"green", 
'eb3_2_clone2':"blue", 
'eb3_2_clone3':"purple",}
print("qqq",df_eb_qn)
df_eb_qn.to_csv("eb32.rt.vert.txt",sep="\t",index=False)
df_eb_qn.to_csv("eb32.rt.vert.bed",sep="\t",header=None,index=False)
os.system("sort -k1,1 -k2,2n eb32.rt.vert.bed > eb32.rt.vert.sorted.bed")
os.system("bedtools map -a eb32.rt.vert.sorted.bed -b ucsc.known.gene.hg19.txn.start.stop.bed.cds.only.first.isoform.nochr.sorted.bed -o distinct -c 4 > eb32.rt.vert.intersect.coding.bed ")

###
for chrom in chromosomes_nochr:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(eb_samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_eb_qn[(df_eb_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=eb_samples[i]+"_hap1_logrt").reset_index()
        maternal = df_eb_qn[(df_eb_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=eb_samples[i]+"_hap2_logrt").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:

            #unsmoothed
            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"],
                    c=color_dict_eb[eb_samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"],
                    c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

            # ## smoothed testing
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2
        ax[0].set_ylim([-3.2,3.2])
        ax[0].set_yticks([-3,-2,-1,0,1,2,3])
        ax[0].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_eb_qn[(df_eb_qn["chrom"]==chrom) & 
                            (df_eb_qn["eb_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_eb_qn[(df_eb_qn["chrom"]==chrom) & (df_eb_qn["arm"]==k)]["start"],
                        df_eb_qn[(df_eb_qn["chrom"]==chrom) & (df_eb_qn["arm"]==k)]["eb_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[1].set_xticks(np.linspace(0,chromosome_length_nochr[chrom],25)) 
    plt.savefig("eb.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()
#### Now do zooms
#####
##### zooms zooms zooms zoosm zooms zoosms zooms zooms zooms 
#####
####
regions = df_eb_qn[df_eb_qn["eb_vert"]==True]
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=500001).to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(2,4),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(eb_samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_eb_qn[(df_eb_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=eb_samples[i]+"_hap1_logrt").reset_index()
        maternal = df_eb_qn[(df_eb_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=eb_samples[i]+"_hap2_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-3000000) & (paternal["stop"]<=row["stop"]+3000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-3000000) & (maternal["stop"]<=row["stop"]+3000000)]

        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:
            ## unsmoothed
            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"],
                    c=color_dict_eb[eb_samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"],
                    c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2


            # #smoothed 
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][eb_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][eb_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_eb[eb_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2


        ax[0].set_ylim([-3.2,3.2])
        ax[0].set_yticks([-3,-2,-1,0,1,2,3])
        ax[0].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_eb_qn[(df_eb_qn["chrom"]==row["chrom"]) & 
                            (df_eb_qn["eb_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_eb_qn[(df_eb_qn["chrom"]==row["chrom"]) & (df_eb_qn["arm"]==k)]["start"],
                        df_eb_qn[(df_eb_qn["chrom"]==row["chrom"]) & (df_eb_qn["arm"]==k)]["eb_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[1].set_xticks(np.linspace(row["start"]-3000000,row["stop"]+3000000,5)) 
    plt.savefig("eb.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()

### intersect VERTs with coding genes


#####
#####
#####
#####
################################
#################################
## GM samples
###################################
##################################
#################################
################################

gm_samples = glob.glob("gm*rt*windows.bed")
dfs=[]
for x in gm_samples:
    x_samp="_".join(x.split("_")[0:3])
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", x_samp+"_hap1_counts", x_samp+"_hap2_counts"])
    # tmp=tmp[tmp["chrom"]!="X"]###removing chrom x 
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 

    tmp=tmp.set_index(["chrom","start","stop"])
    tmp[x_samp+"_cpm_hap1_counts"] = (tmp[x_samp+"_hap1_counts"] ) / tmp[x_samp+"_hap1_counts"].sum(axis="rows")*10**6
    tmp[x_samp+"_cpm_hap2_counts"] = (tmp[x_samp+"_hap2_counts"] ) / tmp[x_samp+"_hap2_counts"].sum(axis="rows")*10**6
    dfs+=[tmp]

df_gm = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
gm_samples = ["_".join(x.split("_")[0:2]) for x in gm_samples]
for sample in gm_samples:
    df_gm[sample+"_hap1_logrt"] = np.log2((df_gm[sample+"_early_cpm_hap1_counts"]+1) / (df_gm[sample+"_late_cpm_hap1_counts"]+1 ))
    df_gm[sample+"_hap2_logrt"] = np.log2((df_gm[sample+"_early_cpm_hap2_counts"]+1) / (df_gm[sample+"_late_cpm_hap2_counts"]+1 ))


df_gm = df_gm.dropna(how="any",axis=0)
df_gm_qn = quantile_normalize(df_gm.filter(regex="logrt"))
df_gm_qn = df_gm_qn.reset_index()
df_gm_qn["gm_std_dev_both_haps"] = df_gm_qn.filter(like="gm",axis=1).std(axis="columns")
mean_std_dev = df_gm_qn[df_gm_qn["chrom"]!="X"]["gm_std_dev_both_haps"].mean()
std_dev_dev = df_gm_qn[df_gm_qn["chrom"]!="X"]["gm_std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_gm_qn = df_gm_qn.sort_values(["chrom","start"]) # this needed for the bedtools map
df_gm_qn = df_gm_qn.set_index(["chrom","start","stop"])
df_gm_qn["gm_vert"] =  df_gm_qn.apply(lambda x:True if x["gm_std_dev_both_haps"]>=threshold else False,axis=1)
df_gm_qn = df_gm_qn.reset_index()
df_gm_qn["arm"] = df_gm_qn.apply(lambda x: "q" if (x["stop"] > arm_dict_nochr[x["chrom"]][0]) & (x["stop"] <= arm_dict_nochr[x["chrom"]][1]) else "p", axis=1)


df_gm_qn.to_csv("gm.rt.vert.txt",sep="\t",index=False)
df_gm_qn.to_csv("gm.rt.vert.bed",sep="\t",header=None,index=False)
os.system("sort -k1,1 -k2,2n gm.rt.vert.bed > gm.rt.vert.sorted.bed")
os.system("bedtools map -a gm.rt.vert.sorted.bed -b ucsc.known.gene.hg19.txn.start.stop.bed.cds.only.first.isoform.nochr.sorted.bed -o distinct -c 4 > gm.rt.vert.intersect.coding.bed ")



color_dict_gm = {'gm12878_clone5':"red", 'gm12878_clone4' :"blue"}

###
for chrom in chromosomes_nochr:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(gm_samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_gm_qn[(df_gm_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=gm_samples[i]+"_hap1_logrt").reset_index()
        maternal = df_gm_qn[(df_gm_qn["chrom"]==chrom)].set_index(["chrom","start","stop","arm"]).filter(like=gm_samples[i]+"_hap2_logrt").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:

            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][gm_samples[i]+"_hap1_logrt"],
                    c=color_dict_gm[gm_samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][gm_samples[i]+"_hap2_logrt"],
                    c=color_dict_gm[gm_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([0,chromosome_length_nochr[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_gm_qn[(df_gm_qn["chrom"]==chrom) & 
                            (df_gm_qn["gm_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_gm_qn[(df_gm_qn["chrom"]==chrom) & (df_gm_qn["arm"]==k)]["start"],
                        df_gm_qn[(df_gm_qn["chrom"]==chrom) & (df_gm_qn["arm"]==k)]["gm_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([0, chromosome_length_nochr[chrom]])
        ax[1].set_xticks(np.linspace(0, chromosome_length_nochr[chrom], 25))
    plt.savefig("gm.rt."+chrom+".png",
        dpi=400,transparent=False)
    plt.close()

#### Now do zooms
#####
##### zooms zooms zooms zoosm zooms zoosms zooms zooms zooms 
#####
####
regions = df_gm_qn[df_gm_qn["gm_vert"]==True]
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=500001).to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(2,1,figsize=(2,4),sharex=False,
                         gridspec_kw={'height_ratios': [6,1]})

    for i in range(len(gm_samples)):
        # print("ordering in plots top to bottom:",i,clones[i])
        paternal = df_gm_qn[(df_gm_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=gm_samples[i]+"_hap1_logrt").reset_index()
        maternal = df_gm_qn[(df_gm_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=gm_samples[i]+"_hap2_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-3000000) & (paternal["stop"]<=row["stop"]+3000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-3000000) & (maternal["stop"]<=row["stop"]+3000000)]

        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")

        for j in ["p","q"]:
            ## unsmoothed
            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][gm_samples[i]+"_hap1_logrt"],
                    c=color_dict_gm[gm_samples[i]],lw=0.5) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][gm_samples[i]+"_hap2_logrt"],
                    c=color_dict_gm[gm_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2


            # #smoothed 
            # ax[0].plot(paternal[paternal["arm"]==j]["start"],
            #         smooth_vector(paternal[paternal["arm"]==j]["start"],
            #             paternal[paternal["arm"]==j][gm_samples[i]+"_hap1_logrt"]),
            #         c=color_dict_gm[gm_samples[i]],lw=0.3) ## -- line style is haplotype 2

            # ax[0].plot(maternal[maternal["arm"]==j]["start"],
            #         smooth_vector(maternal[maternal["arm"]==j]["start"],
            #             maternal[maternal["arm"]==j][gm_samples[i]+"_hap2_logrt"]),
            #         c=color_dict_gm[gm_samples[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in df_gm_qn[(df_gm_qn["chrom"]==row["chrom"]) & 
                            (df_gm_qn["gm_std_dev_both_haps"]>= threshold)].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_gm_qn[(df_gm_qn["chrom"]==row["chrom"]) & (df_gm_qn["arm"]==k)]["start"],
                        df_gm_qn[(df_gm_qn["chrom"]==row["chrom"]) & (df_gm_qn["arm"]==k)]["gm_std_dev_both_haps"],
                        c="black",
                        lw=0.7)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-3000000,row["stop"]+3000000])
        ax[1].set_xticks(np.linspace(row["start"]-3000000,row["stop"]+3000000,5)) 
    plt.savefig("gm.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()


### do the STD DEV work before combining
exit()
###############
print(df_eb_qn,df_gm_qn,df_acp6_qn,df_acp7_qn)

df_combined = pd.concat([df_eb_qn.set_index(["chrom","start","stop","arm"]),
                        df_gm_qn.set_index(["chrom","start","stop","arm"]),
                        df_acp6_qn.set_index(["chrom","start","stop","arm"]),
                        df_acp7_qn.set_index(["chrom","start","stop","arm"])],axis=1).sort_values(["chrom","start"])
df_combined.to_csv("rt.combined.eb.gm.acp6.acp7.txt",sep="\t",index=False)
df_combined.to_csv("rt.combined.eb.gm.acp6.acp7.bed",sep="\t",header=None,na_rep="NaN",index=False)

print(df_combined)


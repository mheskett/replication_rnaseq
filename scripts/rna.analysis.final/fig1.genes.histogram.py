import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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


acp6=pd.read_csv("acp6.with.x.as.gene.counts.txt",sep="\t")#.set_index(["chrom","start","stop"])
acp7=pd.read_csv("acp7.with.x.as.gene.counts.txt",sep="\t")#.set_index(["chrom","start","stop"])
gm=pd.read_csv("gm.with.x.as.gene.counts.txt",sep="\t")#.set_index(["chrom","start","stop"])
eb=pd.read_csv("eb.with.x.as.gene.counts.txt",sep="\t")#.set_index(["chrom","start","stop"])

## eb is in hg19 and has no chr prefix
autosomes = pd.concat([abs(acp6[acp6["chrom"]!="chrX"]["aei"]),
          abs(acp7[acp7["chrom"]!="chrX"]["aei"]),
          abs(gm[gm["chrom"]!="chrX"]["aei"]),
          abs(eb[eb["chrom"]!="X"]["aei"])])


chrx=  pd.concat([abs(acp7[acp7["chrom"]=="chrX"]["aei"]),
       abs(gm[gm["chrom"]=="chrX"]["aei"]), abs(eb[eb["chrom"]=="X"]["aei"])])


# print(abs(acp7[acp7["chrom"]=="chrX"]["aei"]))
# print(abs(gm[gm["chrom"]=="chrX"]["aei"]))
# print(abs(eb[eb["chrom"]=="X"]["aei"]))

print(chrx)
# print(autosomes)
# print(gm[gm["chrom"]=="chrX"])

#### try histogram of counts instead of kde
f,ax=plt.subplots(figsize=(3,1),dpi=300)
# plt.hist(test,bins=75,color="black")#,log=True)
sns.kdeplot(autosomes,color="black",clip=(0,0.5))
sns.kdeplot(chrx,color="red",clip=(0,0.5))
# ax.axvline(x=mean_std_dev + 2.25 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# ax.set_ylim([0,0.75])
# ax.set_yticks([0,0.75])
# ax.set_xlim([0,1.5])
ax.set_xlim([0,0.5])
ax.set_xticks([0,.1,.2,.3,.4,.5])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig("combined.genes.aei.png",bbox_inches='tight')
plt.close()

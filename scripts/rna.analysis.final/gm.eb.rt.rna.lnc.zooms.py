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








## RT

### RNA coding genes
df_eb_genes = pd.read_csv("eb.as.gene.counts.hg38.lifted.txt",sep="\t")

## TLs
df_eb_tls = pd.read_csv("eb.as.tl.counts.all.hg38.lifted.txt",sep="\t")

## repli
df_eb_rt = pd.read_csv("../combined.rt.data.set/eb.rt.hg38.lifted.txt",sep='\t')
# regions = df_eb_rt[(df_eb_rt["eb_vert"]==True) & (df_eb_rt["chrom"]!="chrX")]
# a = pybedtools.BedTool.from_dataframe(regions)
# regions = a.merge(d=250001).filter(lambda x: len(x) > 250000).saveas().to_dataframe().reset_index()
# regions = regions.drop("index",axis="columns")
# regions.columns = ["chrom","start","stop"]
# regions["chrom"] = regions["chrom"].astype(str)
## RT

### RNA coding genes
df_gm_genes = pd.read_csv("gm.as.gene.counts.txt",sep="\t")

## TLs
df_gm_tls = pd.read_csv("gm.as.tl.counts.all.txt",sep="\t")

## repli
df_gm_rt = pd.read_csv("../combined.rt.data.set/gm.rt.txt",sep='\t')
regions = df_gm_rt[(df_gm_rt["gm_vert"]==True) & (df_gm_rt["chrom"]!="chrX")]
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=250001).filter(lambda x: len(x) > 250000).saveas().to_dataframe().reset_index()
regions = regions.drop("index",axis="columns")
regions.columns = ["chrom","start","stop"]
regions["chrom"] = regions["chrom"].astype(str)

color_dict_gm = {
"gm12878_clone4":"red", 
"gm12878_clone5":"blue"}

color_dict_acp6_rna = {
"gm12878_clone4_rnaAligned":"red", 
"gm12878_clone5_rnaAligned":"blue"}

color_dict_acp6_tls = {
"gm12878_clone4":"red", 
"gm12878_clone5":"blue"}


rt_samples_gm = ["gm12878_clone4","gm12878_clone5"]
# rna_samples_gm = ["gm12878_clone4_rnaAligned","gm12878_clone5_rnaAligned"]


rt_samples_eb= ["eb3_2_clone2",
"eb3_2_clone3",
"eb3_2_clone4",
"eb3_2_clone10",
"eb3_2_clone13",
"eb3_2_clone15"]

for index,row in regions.iterrows():
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=5) 
    f, ax = plt.subplots(1,1,figsize=(2,4),sharex=False)
    rect=Rectangle((row["start"], -5), width=row["stop"]-row["start"], height=10,
             facecolor="gray",alpha=0.5,fill=True) 
    ax.add_patch(rect)
    for i in range(len(rt_samples_gm)):

        paternal = df_gm_rt[(df_gm_rt["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=rt_samples_gm[i]+"_paternal_logrt").reset_index()
        maternal = df_gm_rt[(df_gm_rt["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=rt_samples_gm[i]+"_maternal_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-2000000) & (paternal["stop"]<=row["stop"]+2000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-2000000) & (maternal["stop"]<=row["stop"]+2000000)]

        ax.axhline(y=0,lw=0.5, c="black",linestyle="--")

        for j in ["p","q"]:
            ## unsmoothed
            ax.plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][rt_samples_gm[i]+"_paternal_logrt"],
                    c="magenta",lw=0.75) ## -- line style is haplotype 2

            ax.plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][rt_samples_gm[i]+"_maternal_logrt"],
                    c="magenta",linestyle="--",lw=0.75) ## -- line style is haplotype 2
        ax.set_ylim([-4,4])
        ax.set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax.set_xlim([row["start"]-1000000,row["stop"]+1000000]) # needs to be matched with the below xlim
        ax.set_xticks(np.linspace(row["start"]-1000000,row["stop"]+1000000,5)) 
    for i in range(len(rt_samples_eb)):

        paternal = df_eb_rt[(df_eb_rt["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=rt_samples_eb[i]+"_paternal_logrt").reset_index()
        maternal = df_eb_rt[(df_eb_rt["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=rt_samples_eb[i]+"_maternal_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-2000000) & (paternal["stop"]<=row["stop"]+2000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-2000000) & (maternal["stop"]<=row["stop"]+2000000)]

        ax.axhline(y=0,lw=0.5, c="black",linestyle="--")

        ## unsmoothed
        ax.plot(paternal["start"],
                paternal[rt_samples_eb[i]+"_paternal_logrt"],
                c="red",lw=0.75) ## -- line style is haplotype 2

        ax.plot(maternal["start"],
                maternal[rt_samples_eb[i]+"_maternal_logrt"],
                c="red",linestyle="--",lw=0.75) ## -- line style is haplotype 2        



### coding genes

        # ax2 = ax.twinx()
        # rna_tmp = df_acp6_genes[(df_acp6_genes["chrom"]==row["chrom"]) & 
        #                     (df_acp6_genes["start"]>=row["start"]-300000) & 
        #                     (df_acp6_genes["stop"]<=row["stop"]+300000) & 
        #                     (df_acp6_genes["total_reads"]>=10) & 
        #                     (df_acp6_genes["strand_reads"].isin(["plus","minus"]))]

        # for index2, row2 in rna_tmp.iterrows():
        #     rect=Rectangle((row2["start"], row2["aei"]-.0125), 
#                             width=row2["stop"]-row2["start"], height=0.025,
#                             facecolor=color_dict_acp6_rna[row2["sample"]], edgecolor="black",
#                             linestyle='dotted',
#                             fill=True,lw=.5)

#             ax2.add_patch(rect)
# ### TLS
#         tl_tmp = df_acp6_tls[(df_acp6_tls["chrom"]==row["chrom"]) & 
#                             (df_acp6_tls["start"]>=row["start"]-300000) & 
#                             (df_acp6_tls["stop"]<=row["stop"]+300000) & 
#                             (df_acp6_tls["total_reads"]>=10) &
#                             (df_acp6_tls["strand_reads"].isin(["plus","minus"]))]

#         for index2, row2 in tl_tmp.iterrows():
#             rect=Rectangle((row2["start"], row2["aei"]-.0125), 
#                             width=row2["stop"]-row2["start"], height=0.025,
#                             facecolor=color_dict_acp6_tls[row2["sample"]], edgecolor="black",
#                             fill=True,lw=.5)

#             ax2.add_patch(rect)
#         ax2.set_ylim([-0.52,0.52])
#         ax2.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1, 0, .1, .2, .3, .4, 0.5])


        # ax.set_ylim([0,1.5])
        # ax.set_yticks([0,1.5])
        # ax.set_xlim([row["start"]-1200000,row["stop"]+1200000]) # needs to be matched with above xlim
        # ax.set_xticks(np.linspace(row["start"]-1200000,row["stop"]+1200000,5)) 
    # plt.savefig("acp6.rt.genes.tls."+str(row["chrom"])+"-"+str(row["start"])+".png",
    #     dpi=400,transparent=False,bbox_inches='tight')
    plt.suptitle("gm.eb."+str(row["chrom"])+":"+str(row["start"]))

    plt.show()
    plt.close()



# print(df_acp6_genes, df_acp6_tls, df_acp6_rt)








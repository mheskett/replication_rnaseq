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
import glob

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

### REMOVIONG 15 for BAD DATA!!!!
chromosomes = ["chr1","chr2","chr3","chr4","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
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
arm_dict = get_arms(cytoband)
all_files_repli = ["acp7_c2_m_el.log2rt.as.only.bed",
                   "acp7_c2_p_el.log2rt.as.only.bed",
                   "acp7_c4_m_el.log2rt.as.only.bed", 
                   "acp7_c4_p_el.log2rt.as.only.bed", 
                   "acp7_c5_m_el.log2rt.as.only.bed", 
                   "acp7_c5_p_el.log2rt.as.only.bed"]

filenames_repli=[os.path.basename(x)[0:12] for x in all_files_repli]
clones = np.unique([x[0:7] for x in filenames_repli])
repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",
                        names= ["chrom","start","stop","log2rt"],
                        dtype = {"chrom":str,"start":int,"stop":int,"log2rt":float})
    df_repli = df_repli.set_index(["chrom","start","stop"])
    df_repli = df_repli.reset_index()
    df_repli["sample"] = filenames_repli[i]
    repli_li.append(df_repli)

repli_df = pd.concat(repli_li)
repli_df["haplotype"] = ["paternal" if "_p_" in x else "maternal" for x in repli_df["sample"]]
zscore = lambda x: (x - x.mean()) / x.std()
###3

#### BAD DATA REMOVE THIS LATER
repli_df=repli_df[repli_df["chrom"]!="chr5"]
#####
####
###
# print(repli_df)
# print(repli_df[repli_df["haplotype"]=="paternal"].pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any").reset_index())
# print(repli_df[repli_df["haplotype"]=="maternal"].pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(axis="index",how="any").reset_index())
# print(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").reset_index().rename_axis(None, axis=1) )
### .dropna(axis="index",how="any") explore without this, but do it later maybe


# cant normalize with NaNs
repli_df = quantile_normalize(repli_df.pivot(index=["chrom","start","stop"],columns="sample",values="log2rt").dropna(how="any",axis="index")).reset_index().rename_axis(None, axis=1)
###repli_df = repli_df.dropna(how="any",axis=0)
repli_df["std_dev_both_haps"] = repli_df.filter(like="acp",axis=1).std(axis="columns")
repli_df = repli_df[repli_df["chrom"]!="chrX"]
print(repli_df)
mean_std_dev = repli_df["std_dev_both_haps"].mean()
std_dev_dev = repli_df["std_dev_both_haps"].std()
threshold = mean_std_dev + 2.5*std_dev_dev
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)


#### now do the ASRT analysis
samples = []
sample_pair_dict = {}

for i in clones:
    sample_pair_dict[i+"_p_el"]=i+"_m_el"
for i in clones:
    repli_df[i+"_logr_diff_raw"] = repli_df[i+"_p_el"] - repli_df[i+"_m_el"] 
    repli_df[i+"_logr_diff_abs"] = abs(repli_df[i+"_p_el"] - repli_df[i+"_m_el"])

#samples.sort()
#####
####
####
color_dict_repli = {  
"acp7_c2":"cyan",
"acp7_c5":"red",  
"acp7_c4":"green"}
#####
####
####

## get the STD of the averaged out replicates
print(repli_df.filter(like="_p_el"))
repli_df["paternal_std_dev"] = repli_df.filter(like="_p_el",axis=1).std(axis=1)
repli_df["maternal_std_dev"] = repli_df.filter(like="_m_el",axis=1).std(axis=1)
# repli_df["averaged_std_dev"] = repli_df.loc[:,["clone8_cast","clone8_c57","clonee5_cast","clonee5_c57","clone2_cast","clone2_c57"]].std(axis=1)
print(repli_df)
####
####
####
### big all clones RT plot on top, then small subpolots top to bottom: all clones std dev, hap1 std dev, hap2 std dev
### gotta do P Q separation for this
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(4,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1]})
    for i in range(len(clones)):
        print("ordering in plots top to bottom:",i,clones[i])
        paternal = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=clones[i]+"_p_el").reset_index()
        maternal = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=clones[i]+"_m_el").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")
        paternal["arm"]=paternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
        maternal["arm"]=maternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

        for j in ["p","q"]:

            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][clones[i]+"_p_el"],
                    c=color_dict_repli[clones[i]],lw=0.3) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][clones[i]+"_m_el"],
                    c=color_dict_repli[clones[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["std_dev_both_haps"]>=(repli_df["std_dev_both_haps"].mean() + 2*repli_df["std_dev_both_haps"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["paternal_std_dev"]>=(repli_df["paternal_std_dev"].mean() + 2*repli_df["paternal_std_dev"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[2].add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["maternal_std_dev"]>=(repli_df["maternal_std_dev"].mean() + 2*repli_df["maternal_std_dev"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[3].add_patch(rect)

        #####
        for k in ["p","q"]:
            ax[1].plot(repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["start"],
                        repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["std_dev_both_haps"],
                        c="black",
                        lw=0.7)
            ax[2].plot(repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["start"],
                        repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["paternal_std_dev"],
                        c="black",
                        lw=0.7)
            ax[3].plot(repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["start"],
                        repli_df[(repli_df["chrom"]==chrom)&(repli_df["arm"]==k)]["maternal_std_dev"],
                        c="black",
                        lw=0.7)


        ax[1].set_ylim([0,2.5])
        ax[1].set_yticks([0,2])
        ax[1].set_xlim([0,chromosome_length[chrom]])
        ax[1].set_xticks([]) 

        ax[2].set_ylim([0,2.5])
        ax[2].set_yticks([0,2])
        ax[2].set_xlim([0,chromosome_length[chrom]])
        ax[2].set_xticks([]) 

        ax[3].set_ylim([0,2.5])
        ax[3].set_yticks([0,2])
        ax[3].set_xlim([0,chromosome_length[chrom]])
        ax[3].set_xticks(np.linspace(0,chromosome_length[chrom],16))
    plt.savefig("acp7.all.clones.vert.alleles."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()
repli_df.to_csv("acp7.all.clones.txt",sep="\t")


#### this is the ASRT plot i believe
### using average of 3 replicates for this one

###
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(4,1,figsize=(30,3),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1, 1, 1]})
    for i in range(len(clones)):
        print("ordering in plots top to bottom:",i,clones[i])
        paternal = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=clones[i]+"_p_el").reset_index()
        maternal = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=clones[i]+"_m_el").reset_index()
        ax[0].axhline(y=0, linestyle="--" ,lw=0.5, c="black")
        paternal["arm"]=paternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
        maternal["arm"]=maternal.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

        for j in ["p","q"]:

            ax[0].plot(paternal[paternal["arm"]==j]["start"],
                    paternal[paternal["arm"]==j][clones[i]+"_p_el"],
                    c=color_dict_repli[clones[i]],lw=0.3) ## -- line style is haplotype 2

            ax[0].plot(maternal[maternal["arm"]==j]["start"],
                    maternal[maternal["arm"]==j][clones[i]+"_m_el"],
                    c=color_dict_repli[clones[i]],linestyle="--",lw=0.3) ## -- line style is haplotype 2

        ax[0].set_ylim([-4,4])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df[clones[i]+"_logr_diff_abs"]>=(repli_df[clones[i]+"_logr_diff_abs"].mean() + 2*repli_df[clones[i]+"_logr_diff_abs"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="red" if row3[clones[i]+"_logr_diff_raw"]>=0 else "blue",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[i+1].add_patch(rect)

        #####
        ax[i+1].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)][clones[i]+"_logr_diff_abs"],
                    c="black",
                    lw=0.7)

        ax[i+1].set_ylim([0,2.5])
        ax[i+1].set_yticks([0,2])
        ax[i+1].set_xlim([0,chromosome_length[chrom]])
        ax[i+1].set_xticks([]) 


    plt.savefig("acp7.all.clones.avg.asrt.alleles."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()

exit()


#######
#######
## i think noot needed for ACPs here
#######
# 9 sample pairs asrt
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(10,1,figsize=(8,1.5),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1,1,1,1,1,1,1]})
    for i in range(len(samples)):
        hap_cast = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=samples[i]+".CAST.bedGraph").reset_index()
        hap_c57 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=samples[i]+".C57BL.bedGraph").reset_index()
        ax[0].axhline(y=0,linestyle="--",lw=0.5,c="black")
        # print(hap_cast)
        # print(hap_c57)
        ax[0].plot(hap_cast["start"],
                hap_cast[samples[i]+".CAST.bedGraph"],
                c=color_dict_repli[samples[i]],lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(hap_c57["start"],
                hap_c57[samples[i]+".C57BL.bedGraph"],
                c=color_dict_repli[samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-2.5,2.5])
        ax[0].set_yticks([-2,-1,0,1,2])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df[samples[i]+"_logr_diff_abs"]>=(repli_df[samples[i]+"_logr_diff_abs"].mean() + 2.5*repli_df[samples[i]+"_logr_diff_abs"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="red" if row3[samples[i]+"_logr_diff_raw"]>=0 else "blue",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[i+1].add_patch(rect)

        #####
        ax[i+1].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)][samples[i]+"_logr_diff_abs"],
                    c="black",
                    lw=0.7)

        ax[i+1].set_ylim([0,2.5])
        ax[i+1].set_yticks([0,2])
        ax[i+1].set_xlim([0,chromosome_length[chrom]])
        ax[i+1].set_xticks([]) 


    plt.savefig("mouse.f1.israeli.all.clones.asrt.alleles."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()
exit()



#### all samples by chrom, colored by CAST/C57
for j in range(len(chromosomes)):
    f,ax = plt.subplots(figsize=(12,2))

    plt.suptitle(chromosomes[j])
    for i in range(len(cast)):
         ax.plot(repli_df[(repli_df["chrom"]==chromosomes[j])]["start"],
                    repli_df[(repli_df["chrom"]==chromosomes[j])][cast[i]],lw=0.5,c="red")
    for i in range(len(c57)):
        ax.plot(repli_df[(repli_df["chrom"]==chromosomes[j])]["start"],
                    repli_df[(repli_df["chrom"]==chromosomes[j])][c57[i]],lw=0.5,c="blue")
    # for i in range(len(pure)):
    #     ax.plot(repli_df[(repli_df["chrom"]==chromosomes[j])]["start"],
    #                 repli_df[(repli_df["chrom"]==chromosomes[j])][pure[i]],lw=0.5,c="green")
    for index3,row3 in repli_df[(repli_df["chrom"]==chromosomes[j]) & (repli_df["std_dev"]>=threshold)].iterrows():
        rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                 facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
        ax.add_patch(rect)
        # ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_ylim([-2.5,2.5])
    ax.set_yticks([-2,-1,0,1,2])
    ax.set_xlim([0,chromosome_length[chromosomes[j]]])
    ax.set_xticks(np.linspace(0,chromosome_length[chromosomes[j]],16))        # ax.set_xticks([])

    plt.savefig("f1.mouse.israel."+chromosomes[j]+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)



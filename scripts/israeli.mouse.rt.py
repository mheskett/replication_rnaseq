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


chromosome_length ={
"chr1":    197195432,
"chr2":    181748087,
"chr3":    159599783,
"chr4":    155630120,
"chr5":    152537259,
"chr6":    149517037,
"chr7":    152524553,
"chr8":    131738871,
"chr9":    124076172,
"chr10":   129993255,
"chr11":   121843856,
"chr12":   121257530,
"chr13":   120284312,
"chr14":   125194864,
"chr15":   103494974,
"chr16":   98319150,
"chr17":   95272651,
"chr18":   90772031,
"chr19":   61342430,
"chrX":   166650296,
"chrY":    15902555}
chromosomes = ["chr1",
"chr10",
"chr11",
"chr12",
"chr13",
"chr14",
"chr15",
"chr16",
"chr17",
"chr18",
"chr19",
"chr2",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr8",
"chr9",
"chrX"]
# these samples contain NANs in many places
bad_samples = [
"GSM3756321_pre-b_clone3.1.CAST.bedGraph", 
"GSM3756323_pre-b_clone3.3.CAST.bedGraph",      
"GSM3756322_pre-b_clone3.2.CAST.bedGraph",
"GSM3756321_pre-b_clone3.1.C57BL.bedGraph",      
"GSM3756323_pre-b_clone3.3.C57BL.bedGraph",      
"GSM3756322_pre-b_clone3.2.C57BL.bedGraph"]

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
    repli_li += [df_repli.set_index(["chrom","start","stop"])]

repli_df = pd.concat(repli_li,axis=1).reset_index()

# for chrom in chromosomes:
#     print("chrom", chrom)
#     print("chrom", chrom)
#     print("chrom", chrom)
#     print(repli_df[repli_df["chrom"]==chrom].isnull().sum(axis = 0).sort_values(ascending=False))

# print(repli_df.groupby("chrom").count())
repli_df = repli_df.dropna(how="any",axis=0)
# print(repli_df.groupby("chrom").count())
repli_df["std_dev"] = repli_df.filter(like="GSM",axis=1).std(axis="columns")
repli_df_auto = repli_df[repli_df["chrom"]!="chrX"]

mean_std_dev = repli_df_auto["std_dev"].mean()
print(mean_std_dev)
std_dev_dev = repli_df_auto["std_dev"].std()
print(std_dev_dev)
threshold = mean_std_dev + 2.5*std_dev_dev
print(threshold)


#### now do the ASRT analysis
samples = []
sample_pair_dict = {}
repli_df["clone8_cast"] = repli_df.filter(like="GSM3756326_pre-b_clone8",axis=1).filter(like="CAST").mean(axis="columns")
repli_df["clone8_c57"] = repli_df.filter(like="GSM3756326_pre-b_clone8",axis=1).filter(like="C57BL").mean(axis="columns")
repli_df["clonee5_cast"] = repli_df.filter(like="GSM3756315_pre-b_clone_e5",axis=1).filter(like="CAST").mean(axis="columns")
repli_df["clonee5_c57"] = repli_df.filter(like="GSM3756315_pre-b_clone_e5",axis=1).filter(like="C57BL").mean(axis="columns")
repli_df["clone2_cast"] = repli_df.filter(like="GSM3756318_pre-b_clone2.1",axis=1).filter(like="CAST").mean(axis="columns")
repli_df["clone2_c57"] = repli_df.filter(like="GSM3756318_pre-b_clone2.1",axis=1).filter(like="C57BL").mean(axis="columns")
print(repli_df)
for i in cast:
    samples += ['.'.join(i.split(".")[0:2])]
    pair = '.'.join(i.split(".")[0:2])+".C57BL.bedGraph"
    sample_pair_dict[i]=pair
for i in samples:
    repli_df[i+"_logr_diff_raw"] = repli_df[i+".CAST.bedGraph"] - repli_df[i+".C57BL.bedGraph"] 
    repli_df[i+"_logr_diff_abs"] = abs(repli_df[i+".CAST.bedGraph"] - repli_df[i+".C57BL.bedGraph"])
samples.sort()
#####
####
####
color_dict_repli = {
'GSM3756326_pre-b_clone8.3':"blue",
'GSM3756325_pre-b_clone8.2':"blue",
'GSM3756324_pre-b_clone8.1':"blue",
'GSM3756315_pre-b_clone_e5.1':"green",
'GSM3756316_pre-b_clone_e5.2':"green",
'GSM3756317_pre-b_clone_e5.3':"green",
'GSM3756318_pre-b_clone2.1':"red",
'GSM3756319_pre-b_clone2.2':"red",
'GSM3756320_pre-b_clone2.3':"red"
}

average_color_dict_repli = {
"clone8":"blue","clonee5":"green","clone2":"red"
}
#####
####
####
#### lastly do VERT with averaged out samples and hap1hap2 vert, hap1 vert, hap2 vert
average_samples = ["clone8","clonee5","clone2"]
for i in average_samples:
    repli_df[i+"_logr_diff_raw"] = repli_df[i+"_cast"] - repli_df[i+"_c57"] 
    repli_df[i+"_logr_diff_abs"] = abs(repli_df[i+"_cast"] - repli_df[i+"_c57"])

## get the STD of the averaged out replicates
print(repli_df.filter(like="_cast"))
repli_df["cast_std_dev"] = repli_df.filter(like="_cast",axis=1).std(axis=1)
repli_df["c57_std_dev"] = repli_df.filter(like="_c57",axis=1).std(axis=1)
print(repli_df.loc[:,["clone8_cast","clone8_c57","clonee5_cast","clonee5_c57","clone2_cast","clone2_c57"]])
repli_df["averaged_std_dev"] = repli_df.loc[:,["clone8_cast","clone8_c57","clonee5_cast","clonee5_c57","clone2_cast","clone2_c57"]].std(axis=1)
print(repli_df)
###
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(4,1,figsize=(8,1.5),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1]})
    for i in range(len(average_samples)):
        hap_cast = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=average_samples[i]+"_cast").reset_index()
        hap_c57 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=average_samples[i]+"_c57").reset_index()
        ax[0].axhline(y=0,linestyle="--",lw=0.5,c="black")
        # print(hap_cast)
        # print(hap_c57)
        ax[0].plot(hap_cast["start"],
                hap_cast[average_samples[i]+"_cast"],
                c=average_color_dict_repli[average_samples[i]],lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(hap_c57["start"],
                hap_c57[average_samples[i]+"_c57"],
                c=average_color_dict_repli[average_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-2.2,2.2])
        ax[0].set_yticks([-2,-1,0,1,2])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["averaged_std_dev"]>=(repli_df["averaged_std_dev"].mean() + 2.5*repli_df["averaged_std_dev"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[1].add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["cast_std_dev"]>=(repli_df["cast_std_dev"].mean() + 2.5*repli_df["cast_std_dev"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[2].add_patch(rect)
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df["c57_std_dev"]>=(repli_df["c57_std_dev"].mean() + 2.5*repli_df["c57_std_dev"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[3].add_patch(rect)

        #####
        ax[1].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)]["averaged_std_dev"],
                    c="black",
                    lw=0.7)
        ax[2].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)]["cast_std_dev"],
                    c="black",
                    lw=0.7)
        ax[3].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)]["c57_std_dev"],
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
    plt.savefig("mouse.f1.israeli.all.clones.avg.vert.alleles."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

repli_df.to_csv("mouse.israeli.processed.data.clones.txt",sep="\t")



### using average of 3 replicates for this one
average_samples = ["clone8","clonee5","clone2"]
for i in average_samples:
    repli_df[i+"_logr_diff_raw"] = repli_df[i+"_cast"] - repli_df[i+"_c57"] 
    repli_df[i+"_logr_diff_abs"] = abs(repli_df[i+"_cast"] - repli_df[i+"_c57"])
###
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(4,1,figsize=(8,1.5),sharex=False,
                         gridspec_kw={'height_ratios': [6, 1,1,1]})
    for i in range(len(average_samples)):
        hap_cast = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=average_samples[i]+"_cast").reset_index()
        hap_c57 = repli_df[(repli_df["chrom"]==chrom)].set_index(["chrom","start","stop"]).filter(like=average_samples[i]+"_c57").reset_index()
        ax[0].axhline(y=0,linestyle="--",lw=0.5,c="black")
        # print(hap_cast)
        # print(hap_c57)
        ax[0].plot(hap_cast["start"],
                hap_cast[average_samples[i]+"_cast"],
                c=average_color_dict_repli[average_samples[i]],lw=0.5) ## -- line style is haplotype 2

        ax[0].plot(hap_c57["start"],
                hap_c57[average_samples[i]+"_c57"],
                c=average_color_dict_repli[average_samples[i]],linestyle="--",lw=0.5) ## -- line style is haplotype 2

        ax[0].set_ylim([-2.5,2.5])
        ax[0].set_yticks([-2,-1,0,1,2])
        ax[0].set_xlim([0,chromosome_length[chrom]])
        ax[0].set_xticks([])

        #### highlighting allele variant regions
        for index3,row3 in repli_df[(repli_df["chrom"]==chrom) & 
                            (repli_df[average_samples[i]+"_logr_diff_abs"]>=(repli_df[average_samples[i]+"_logr_diff_abs"].mean() + 2.5*repli_df[average_samples[i]+"_logr_diff_abs"].std()))].iterrows():
            rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
                     facecolor="red" if row3[average_samples[i]+"_logr_diff_raw"]>=0 else "blue",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
            ax[i+1].add_patch(rect)

        #####
        ax[i+1].plot(repli_df[(repli_df["chrom"]==chrom)]["start"],
                    repli_df[(repli_df["chrom"]==chrom)][average_samples[i]+"_logr_diff_abs"],
                    c="black",
                    lw=0.7)

        ax[i+1].set_ylim([0,2.5])
        ax[i+1].set_yticks([0,2])
        ax[i+1].set_xlim([0,chromosome_length[chrom]])
        ax[i+1].set_xticks([]) 


    plt.savefig("mouse.f1.israeli.all.clones.avg.asrt.alleles."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()


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



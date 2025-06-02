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



chromosomes_hg38 = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

chromosome_length_mm39 = {
"chr1":    195154279,
"chr2":    181755017,
"chrX":    169476592,
"chr3":    159745316,
"chr4":    156860686,
"chr5":    151758149,
"chr6":    149588044,
"chr7":    144995196,
"chr10":   130530862,
"chr8":    130127694,
"chr14":   125139656,
"chr9":    124359700,
"chr11":   121973369,
"chr13":   120883175,
"chr12":   120092757,
"chr15":   104073951,
"chr16":   98008968,
"chr17":   95294699,
"chrY":    91455967,
"chr18":   90720763,
"chr19":   61420004,}

chromosome_length_hg38 = {"chr1":    248956422,
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

## mouse here
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
# bad_samples = [
# "GSM3756321_pre-b_clone3.1.CAST.bedGraph", 
# "GSM3756323_pre-b_clone3.3.CAST.bedGraph",      
# "GSM3756322_pre-b_clone3.2.CAST.bedGraph",
# "GSM3756321_pre-b_clone3.1.C57BL.bedGraph",      
# "GSM3756323_pre-b_clone3.3.C57BL.bedGraph",      
# "GSM3756322_pre-b_clone3.2.C57BL.bedGraph"]

all_files_repli = glob.glob("GSM*")
# all_files_repli = [x for x in all_files_repli if x not in bad_samples]
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

print(repli_df)

repli_df = repli_df.set_index(["chrom","start","stop"])
samples = []
sample_pair_dict = {}
repli_df["clone8_cast"] = repli_df.filter(like="pre-b_clone8",axis=1).filter(like="CAST").mean(axis="columns",skipna=True)
repli_df["clone8_c57"] = repli_df.filter(like="pre-b_clone8",axis=1).filter(like="C57BL").mean(axis="columns",skipna=True)
repli_df["clonee5_cast"] = repli_df.filter(like="pre-b_clone_e5",axis=1).filter(like="CAST").mean(axis="columns",skipna=True)
repli_df["clonee5_c57"] = repli_df.filter(like="pre-b_clone_e5",axis=1).filter(like="C57BL").mean(axis="columns",skipna=True)
repli_df["clone2_cast"] = repli_df.filter(like="pre-b_clone2",axis=1).filter(like="CAST").mean(axis="columns",skipna=True)
repli_df["clone2_c57"] = repli_df.filter(like="pre-b_clone2",axis=1).filter(like="C57BL").mean(axis="columns",skipna=True)
repli_df["clone3_c57"] = repli_df.filter(like="pre-b_clone3",axis=1).filter(like="C57BL").mean(axis="columns",skipna=True)
repli_df["clone3_cast"] = repli_df.filter(like="pre-b_clone3",axis=1).filter(like="CAST").mean(axis="columns",skipna=True)


repli_df= repli_df.loc[:,["clone8_cast","clone8_c57","clonee5_cast","clonee5_c57","clone2_cast","clone2_c57","clone3_cast","clone3_c57"]].sort_values(["chrom","start"])



repli_df["mouse_std_dev"] = repli_df.std(axis=1)
repli_df["mouse_std_dev_ln"] = np.log(repli_df["mouse_std_dev"])
repli_df= repli_df.reset_index()
mean_std_dev = repli_df[repli_df["chrom"]!="chrX"]["mouse_std_dev_ln"].mean()
std_dev_dev = repli_df[repli_df["chrom"]!="chrX"]["mouse_std_dev_ln"].std()
threshold = mean_std_dev + 2*std_dev_dev

###
f,ax=plt.subplots(figsize=(2,2),dpi=300)
sns.kdeplot(repli_df["mouse_std_dev"],clip=(0,20),linewidth=2)
ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="black")
## plt.show()

repli_df["mouse_vert"] =  repli_df.apply(lambda x:True if x["mouse_std_dev_ln"]>=threshold else False,axis=1)

### mm9 to hg19 then lift hg19 to hg38
### make a liftover python function that inputs a df and outputs a df?
repli_df.to_csv("mouse.rt.mm9.bed",sep="\t",header=None,na_rep="NaN",index=False)
repli_df.to_csv("mouse.rt.mm9.txt",sep="\t",na_rep="NaN",index=False)

## liftover
os.system("/Users/michaelheskett/replication_rnaseq/scripts/liftover_files/liftOver mouse.rt.mm9.bed /Users/michaelheskett/replication_rnaseq/scripts/liftover_files/mm9ToHg19.over.chain mouse.rt.hg19.lifted.bed mouse.mm9.hg19.unmapped.bed -bedPlus=3")
os.system("/Users/michaelheskett/replication_rnaseq/scripts/liftover_files/liftOver mouse.rt.mm9.bed /Users/michaelheskett/replication_rnaseq/scripts/liftover_files/mm9ToMm39.over.chain mouse.rt.mm39.lifted.bed mouse.mm9.mm39.unmapped.bed -bedPlus=3")
os.system("/Users/michaelheskett/replication_rnaseq/scripts/liftover_files/liftOver mouse.rt.mm39.lifted.bed /Users/michaelheskett/replication_rnaseq/scripts/liftover_files/mm39ToHg38.over.chain mouse.rt.mm39.hg38.lifted.bed mouse.mm39.hg38.unmapped.bed -bedPlus=3")

os.system("grep -v \"#\" mouse.mm39.hg38.unmapped.bed > mouse.mm39.hg38.unmapped.fixed.bed")
repli_df_hg38_lifted = pd.read_csv("mouse.mm39.hg38.unmapped.fixed.bed",sep="\t",names=repli_df.columns)
repli_df_hg38_lifted[repli_df_hg38_lifted["mouse_vert"]==True].loc[:,["chrom","start","stop"]].to_csv("mouse.rt.vert.sorted.4col.bed",index=False,header=False,sep="\t")
os.system("bedtools merge -i mouse.rt.vert.sorted.4col.bed > mouse.rt.vert.sorted.4col.merged.bed")

repli_df_mm39_lifted = pd.read_csv("mouse.rt.mm39.lifted.bed",sep="\t",names=repli_df.columns)
repli_df_mm39_lifted.to_csv("mouse.rt.mm39.lifted.txt",sep="\t",index=None)
### get mouse genes 
os.system("sort -k1,1 -k2,2n mouse.rt.mm9.bed | grep True | awk '$1!=\"chrX\"{print $0}' > mouse.rt.vert.sorted.bed")
os.system("bedtools map -a mouse.rt.vert.sorted.bed -b mm9.refseq.cds.only.first.isoform.bed -o distinct -c 4 > mouse.rt.vert.intersect.coding.bed ")
os.system("awk '{print $15}'  mouse.rt.vert.intersect.coding.bed | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v Rik | grep -v ^Gm | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > mouse.vert.genes.txt")

## get mouse genes mm39
os.system("sort -k1,1 -k2,2n mouse.rt.mm39.lifted.bed | grep True | awk '$1!=\"chrX\"{print $0}' > mouse.rt.vert.sorted.mm39.bed")
os.system("bedtools map -a mouse.rt.vert.sorted.mm39.bed -b ../ucsc.refseq.mm39.txn.whole.gene.sorted.bed -o distinct -c 4 > mouse.rt.vert.intersect.coding.mm39.bed ")
os.system("awk '{print $15}'  mouse.rt.vert.intersect.coding.mm39.bed | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v Rik | grep -v ^Gm | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > mouse.vert.genes.mm39.txt")

### plot the actual repli-seq tracks
plot_samples=repli_df_mm39_lifted.filter(like="clone",axis=1).columns
for chrom in chromosomes:

    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp=repli_df_mm39_lifted[repli_df_mm39_lifted["chrom"]==chrom]
    # tmp["color"]=["red" if "cast" in x else "blue" for x in tmp["mouse_std_dev_ln"]]
    for sample in plot_samples:
        print(sample)

        ax.plot(tmp["start"],tmp[sample],c="red" if "cast" in sample else "blue",lw=1)

    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length_mm39[chrom],25)) 

    # ax.set_ylim([0,1.2])
    plt.savefig("mouse.rt.mm38."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()


exit()
print(repli_df_mm39_lifted)
## plot in native mouse chromosomes
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp=repli_df_mm39_lifted[repli_df_mm39_lifted["chrom"]==chrom]
    tmp["color"]=["darkorange" if x>=threshold else "lightgray" for x in tmp["mouse_std_dev_ln"]]
    
    ax.scatter(tmp["start"],
        tmp["mouse_std_dev"],c=tmp["color"],
        s=15,edgecolor="black",lw=0.2)

    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length_mm39[chrom],25)) 

    ax.set_ylim([0,1.2])
    plt.savefig("mouse.std.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()


## plot in lifted over human chromosomes
###
###
###
for chrom in chromosomes_hg38:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp=repli_df_hg38_lifted[repli_df_hg38_lifted["chrom"]==chrom]
    tmp["color"]=["darkorange" if x>=threshold else "lightgray" for x in tmp["mouse_std_dev_ln"]]

    ax.scatter(tmp["start"],
        tmp["mouse_std_dev"],c=tmp["color"],
        s=15,edgecolor="black",lw=0.2)

    ax.set_xlim([0,chromosome_length_hg38[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length_hg38[chrom],25)) 

    ax.set_ylim([0,2.5])
    plt.savefig("mouse.to.hg38.std.rt."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()


exit()






















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



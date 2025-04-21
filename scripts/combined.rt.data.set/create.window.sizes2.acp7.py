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

# df_haps = pd.read_csv("acp6_joint.fb30.all.chrom.no.missing.wh.phased.final.hets.bed",sep="\t",
# 	names=["chrom","start","stop","paternal","maternal"])
# # call bedtools makewindows and do map count
# # os.system("bedtools makewindows -g hg38.fa.fai -w 250000 -s 250000 > hg38.windows.w250.s250.bed")

# # os.system("bedtools intersect -a hg38.windows.w250.s250.bed \
# # 	-b acp6_joint.fb30.all.chrom.no.missing.wh.phased.final.hets.bed -c > hg38.windows.intersect.acp6.snps.bed")

# df2 = pd.read_csv("hg38.windows.intersect.acp6.snps.bed",sep="\t",
# 	names=["chrom","start","stop","count"],dtype={"count":int})
# df2=df2.sort_values(["chrom","start"])
# print(df2)
# # sns.kdeplot(df2["count"],clip=(0,1000))
# # plt.show()

# dfs=[]
# acp6_samples = glob.glob("acp6*hg38*counts.cpms.rmv.blck.bed")
# for x in acp6_samples:
#     samp = os.path.basename(x).split(".")[0]
#     # print(samp)
#     tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", "paternal","maternal",samp+"_p_counts", samp+"_m_counts"])
#     tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
#     tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 
#     tmp=tmp.drop(["paternal","maternal"],axis=1)
#     tmp = tmp[(tmp[samp+"_p_counts"]!=".") & (tmp[samp+"_m_counts"]!=".")]
#     tmp[samp+"_p_counts"] = tmp[samp+"_p_counts"].astype(float)
#     tmp[samp+"_m_counts"] = tmp[samp+"_m_counts"].astype(float)
#     # tmp = tmp[tmp[samp+"_cpm_p_counts"]>=1]
#     # tmp = tmp[tmp[samp+"_cpm_m_counts"]>=1]
#     # tmp=tmp.set_index(["chrom","start","stop"])
#     dfs+=[tmp]
# df = dfs[0]
# samples = [os.path.basename(x).split(".")[0] for x in acp6_samples]

# ## merge the DFs and set index
# for d in dfs[1:]:
#     df = df.merge(d, how="outer",left_on=["chrom","start","stop"],right_on=["chrom","start","stop"])
# df=df.set_index(["chrom","start","stop"])

# ## replace NaNs with 0s
# df = df.replace(np.nan, 0)

# ## remove rows with all 0's
# df = df[df.sum(axis=1)>0]
# df=df.reset_index()
# df.to_csv("acp6_hg38.all.clones.repli.bed",sep="\t",header=None,index=None)

## loop through files and do intersections
# bedtools sort -g hg38.fa.fai -i $1 | 
#               bedtools map -a hg38.windows.w250.s250.bed\
#               -b stdin -o sum,sum -c 6,7  > $filename.allele.counts.windows.s125.bed

#####
#### Replication timing Replication timing Replication timing Replication timing Replication timing 
acp7_samples = glob.glob("./acp7*hg38*counts.rmv.blck.cpms.windows.w250.s250.bed")

dfs=[]
for x in acp7_samples:
    samp = os.path.basename(x).split(".")[0]
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", samp+"_cpm_p_counts", samp+"_cpm_m_counts"])
    tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 
    # tmp = tmp[(tmp[samp+"_cpm_p_counts"]!=".") & (tmp[samp+"_cpm_m_counts"]!=".")]
    tmp = tmp.replace(".", 0)
    tmp[samp+"_cpm_p_counts"] = tmp[samp+"_cpm_p_counts"].astype(float)
    tmp[samp+"_cpm_m_counts"] = tmp[samp+"_cpm_m_counts"].astype(float)
    # tmp = tmp[tmp[samp+"_cpm_p_counts"]>=1]
    # tmp = tmp[tmp[samp+"_cpm_m_counts"]>=1]
    tmp=tmp.set_index(["chrom","start","stop"])
    dfs+=[tmp]
df_acp7 = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
samples=list(set(["_".join(os.path.basename(x).split(".")[0].split("_")[0:2]) for x in acp7_samples]))

#sort for plotting
df_acp7=df_acp7.sort_values(["chrom","start"])

## replace NaNs with 0s
df_acp7 = df_acp7.replace(".", 0)

## remove rows with all 0's
df_acp7 = df_acp7[df_acp7.sum(axis=1)>0]

# plot counts per window
print(df_acp7)
sns.kdeplot(df_acp7.sum(axis=1),clip=(0,5000))
# plt.show()

## according to the plot, remove all windows with less than ~150 LSM counts
df_acp7 = df_acp7[df_acp7.sum(axis=1)>150]

print(samples)
## now do E/L within the windows
for sample in samples:
    df_acp7[sample+"_paternal_logrt"] = np.log2((df_acp7[sample+"_early_cpm_p_counts"]+1) / (df_acp7[sample+"_late_cpm_p_counts"]+1 ))
    df_acp7[sample+"_maternal_logrt"] = np.log2((df_acp7[sample+"_early_cpm_m_counts"]+1) / (df_acp7[sample+"_late_cpm_m_counts"]+1 ))

df_acp7_repli = df_acp7.filter(like="logrt")
print(df_acp7_repli)

### process and plot
df_acp7_qn = quantile_normalize(df_acp7.filter(regex="logrt"))
df_acp7_qn["acp7_std_dev_both_haps"] = df_acp7_qn.filter(like="acp7",axis=1).std(axis="columns")
df_acp7_qn["acp7_std_dev_both_haps_ln"] = np.log(df_acp7_qn["acp7_std_dev_both_haps"])
df_acp7_qn = df_acp7_qn.reset_index()
mean_std_dev = df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps_ln"].mean()
std_dev_dev = df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps_ln"].std()
threshold = mean_std_dev + 2.25 *std_dev_dev
df_acp7_qn = df_acp7_qn.set_index(["chrom","start","stop"])
df_acp7_qn["acp7_vert"] =  df_acp7_qn.apply(lambda x:True if x["acp7_std_dev_both_haps_ln"]>=threshold else False,axis=1)

df_acp7_qn = df_acp7_qn.reset_index()
df_acp7_qn["arm"] = df_acp7_qn.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

print(threshold)
print("acp7 225",df_acp7_qn)
##### make distribution of repliseq std dev. female sample
f,ax=plt.subplots(figsize=(2,2),dpi=300)
sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=1.5,c="black")
sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]=="chrX"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=1.5,c="black",linestyle="--")
# ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
# ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")
plt.savefig("acp7.vert.dist.png")

# plt.show()
plt.close()
#### more informative histograms for paper
####
####


#### try histogram of counts instead of kde
f,ax=plt.subplots(figsize=(3,1),dpi=300)
plt.hist(df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps_ln"],bins=50,color="black")#,log=True)
ax.axvline(x=mean_std_dev + 2.25 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# ax.set_ylim([0,0.75])
# ax.set_yticks([0,0.75])
# ax.set_xlim([0,1.5])
ax.set_xlim([-3.5,1.5])
ax.set_xticks([-3,-2,-1,0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig("acp7.vert.dist.autosomes.png",bbox_inches='tight')
plt.close()

f,ax=plt.subplots(figsize=(3,1),dpi=300)
plt.hist(df_acp7_qn[df_acp7_qn["chrom"]=="chrX"]["acp7_std_dev_both_haps_ln"],bins=50,color="red")#,log=True)
ax.axvline(x=mean_std_dev + 2.25 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# ax.set_ylim([0,0.75])
# ax.set_yticks([0,0.75])
ax.set_xlim([-3.5,1.5])
ax.set_xticks([-3,-2,-1,0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig("acp7.vert.dist.x.png",bbox_inches='tight')
plt.close()

# ### test auto and X on same
# f,ax=plt.subplots(figsize=(3,1),dpi=300)
# plt.hist(df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps_ln"],bins=50,color="black",log=True)
# plt.hist(df_acp7_qn[df_acp7_qn["chrom"]=="chrX"]["acp7_std_dev_both_haps_ln"],bins=50,color="red",log=True)
# # ax.set_ylim([0,0.75])
# # ax.set_yticks([0,0.75])
# # ax.set_xlim([0,1.5])
# ax.spines[['right', 'top']].set_visible(False)
# plt.savefig("acp7.vert.dist.x.auto.png",bbox_inches='tight')
# plt.close()

# #### autosomes only
# f,ax=plt.subplots(figsize=(3,1),dpi=300)
# sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=1.5,c="black")
# ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# ax.set_ylim([0,0.75])
# ax.set_yticks([0,0.75])
# ax.set_xlim([0,4])
# plt.savefig("acp7.vert.dist.autosomes.png",bbox_inches='tight')
# plt.close()

# ## top half broken axis
# f,ax=plt.subplots(figsize=(3,.5),dpi=300)
# sns.kdeplot(df_acp7_qn[df_acp7_qn["chrom"]!="chrX"]["acp7_std_dev_both_haps"],clip=(0,20),linewidth=1.5,c="black")
# ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")
# ax.set_ylim([2.5,3.0])
# ax.set_yticks([2.5,3])
# ax.set_xlim([0,4])
# ax.set_xticks([])
# plt.savefig("acp7.vert.dist.autosomes.top.png",bbox_inches='tight')
# plt.close()



color_dict_acp7 = {
"acp7_c2":"red", 
"acp7_c4":"green",  
"acp7_c5":"blue"}

df_acp7_qn = df_acp7_qn.sort_values(["chrom","start"])

#acp7 should be green
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp=df_acp7_qn[df_acp7_qn["chrom"]==chrom]
    tmp["color"]=["lime" if x>=threshold else "lightgray" for x in tmp["acp7_std_dev_both_haps_ln"]]
    for j in ["p","q"]:
         ax.scatter(tmp[tmp["arm"]==j]["start"],
			tmp[tmp["arm"]==j]["acp7_std_dev_both_haps"],c=tmp[tmp["arm"]==j]["color"],
			s=15,edgecolor="black",lw=0.2)

    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length[chrom],25)) 

    ax.set_ylim([0,2.1])
    ax.set_yticks([0,.5,1,1.5,2])
    plt.savefig("acp7.std.rt."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

df_acp7_qn.to_csv("acp7.rt.txt",sep="\t",index=False)
df_acp7_qn.to_csv("acp7.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)
## decide whether to keep coding and noncoding or not...locsxxx, SNERPAs, lncrnas, mirnas, etc
os.system("sort -k1,1 -k2,2n acp7.rt.bed | grep True | awk '$1!=\"chrX\"{print $0}' > acp7.rt.vert.sorted.bed")
os.system("bedtools map -a acp7.rt.vert.sorted.bed -b ucsc.refseq.hg38.txn.whole.gene.sorted.bed -o distinct -c 4 > acp7.rt.vert.intersect.coding.bed")
os.system("awk '{print $14}'  acp7.rt.vert.intersect.coding.bed | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > acp7.vert.genes.txt")


df_acp7_qn[df_acp7_qn["acp7_vert"]==True].loc[:,["chrom","start","stop","acp7_vert"]].to_csv("acp7.rt.vert.sorted.4col.bed",sep="\t",header=False,index=None)
os.system("bedtools merge -i acp7.rt.vert.sorted.4col.bed > acp7.rt.vert.sorted.4col.merged.bed")





# zooms zooms zooms
regions = df_acp7_qn[df_acp7_qn["acp7_vert"]==True]
print("acp6 regions", regions)
a = pybedtools.BedTool.from_dataframe(regions)
regions = a.merge(d=250001).filter(lambda x: len(x) > 250000).saveas().to_dataframe().reset_index()
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
        paternal = df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_paternal_logrt").reset_index()
        maternal = df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"])].set_index(["chrom","start","stop","arm"]).filter(like=samples[i]+"_maternal_logrt").reset_index()

        paternal = paternal[(paternal["chrom"]==row["chrom"]) & (paternal["start"]>=row["start"]-2000000) & (paternal["stop"]<=row["stop"]+2000000)]
        maternal = maternal[(maternal["chrom"]==row["chrom"]) & (maternal["start"]>=row["start"]-2000000) & (maternal["stop"]<=row["stop"]+2000000)]

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

        ax[0].set_ylim([-4.2,4.2])
        ax[0].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
        ax[0].set_xlim([row["start"]-1500000,row["stop"]+1500000]) # needs to be matched with the below xlim
        ax[0].set_xticks([])

        ##### coding genes
        #######
        ## currently pklotting the total and the plus/minus and the antisense. so pick one to plot....
        # ax2 = ax[0].twinx()
        # rna_tmp = df_acp6rna[(df_acp6rna["chrom"]==row["chrom"]) & 
        #                     (df_acp6rna["start"]>=row["start"]-300000) & 
        #                     (df_acp6rna["stop"]<=row["stop"]+300000) & 
        #                     (df_acp6rna["total_reads"]>=10) & 
        #                     (df_acp6rna["reads_strand"].isin(["plus","minus"]))]

        # for index2, row2 in rna_tmp.iterrows():
        #     rect=Rectangle((row2["start"], row2["aei"]-.0125), 
        #                     width=row2["stop"]-row2["start"], height=0.025,
        #                     facecolor=color_dict_acp6_rna[row2["sample"]], edgecolor="black",
        #                     fill=True,lw=.2)

            # ax2.add_patch(rect)
        # ax2.set_ylim([-0.52,0.52])
        # ax2.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1, 0, .1, .2, .3, .4, 0.5])
        ######
        ######
        #### highlighting allele variant regions
        #### no. this should be the regions. 

        rect=Rectangle((row["start"], -5), width=row["stop"]-row["start"], height=10,
                 facecolor="gray",alpha=1,fill=True) 
        ax[1].add_patch(rect)

        for k in ["p","q"]:

            ax[1].plot(df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"]) & (df_acp7_qn["arm"]==k)]["start"],
                        df_acp7_qn[(df_acp7_qn["chrom"]==row["chrom"]) & (df_acp7_qn["arm"]==k)]["acp7_std_dev_both_haps"],
                        c="black",
                        lw=0.5)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-1200000,row["stop"]+1200000]) # needs to be matched with above xlim
        ax[1].set_xticks(np.linspace(row["start"]-1200000,row["stop"]+1200000,5)) 
    plt.savefig("acp7.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()
exit()


# fig,ax=plt.subplots(1,1)
# ax.scatter(df2["count"],df_acp6.sum(axis=1),s=10)
# ax.set_xlim([0,500])
# ax.set_ylim([0,3500])
# plt.show()

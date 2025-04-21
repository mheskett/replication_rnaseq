import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.stats.multitest as mt

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
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]


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



# ## 
#####
#### Replication timing Replication timing Replication timing Replication timing Replication timing 
acp6_samples = glob.glob("./acp6*hg38*counts.rmv.blck.cpms.windows.w250.s250.bed")

dfs=[]
for x in acp6_samples:
    samp = os.path.basename(x).split(".")[0]
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", samp+"_cpm_p_counts", samp+"_cpm_m_counts"])
    # tmp=tmp[tmp["chrom"]!="chrX"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="chrY"]###removing chrom x 
    # tmp = tmp[(tmp[samp+"_cpm_p_counts"]!=".") & (tmp[samp+"_cpm_m_counts"]!=".")]
    tmp = tmp.replace(".", 0)
    tmp[samp+"_cpm_p_counts"] = tmp[samp+"_cpm_p_counts"].astype(float)
    tmp[samp+"_cpm_m_counts"] = tmp[samp+"_cpm_m_counts"].astype(float)
    # tmp = tmp[tmp[samp+"_cpm_p_counts"]>=1]
    # tmp = tmp[tmp[samp+"_cpm_m_counts"]>=1]
    tmp=tmp.set_index(["chrom","start","stop"])
    dfs+=[tmp]
df_acp6 = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
samples=list(set(["_".join(os.path.basename(x).split(".")[0].split("_")[0:2]) for x in acp6_samples]))

#sort for plotting
df_acp6=df_acp6.sort_values(["chrom","start"])

## replace NaNs with 0s
df_acp6 = df_acp6.replace(".", 0)

## remove rows with all 0's
df_acp6 = df_acp6[df_acp6.sum(axis=1)>0]

# plot counts per window
print(df_acp6)
# sns.kdeplot(df_acp6.sum(axis=1),clip=(0,5000))
# plt.show()

## according to the plot, remove all windows with less than ~150 LSM counts
df_acp6 = df_acp6[df_acp6.sum(axis=1)>150]

print(samples)
## now do E/L within the windows
for sample in samples:
    df_acp6[sample+"_paternal_logrt"] = np.log2((df_acp6[sample+"_early_cpm_p_counts"]+1) / (df_acp6[sample+"_late_cpm_p_counts"]+1 ))
    df_acp6[sample+"_maternal_logrt"] = np.log2((df_acp6[sample+"_early_cpm_m_counts"]+1) / (df_acp6[sample+"_late_cpm_m_counts"]+1 ))

df_acp6_repli = df_acp6.filter(like="logrt")
print(df_acp6_repli)

### process and plot
df_acp6_qn = quantile_normalize(df_acp6.filter(regex="logrt"))
df_acp6_qn["acp6_std_dev_both_haps"] = df_acp6_qn.filter(like="acp6",axis=1).std(axis="columns")
df_acp6_qn["acp6_std_dev_both_haps_ln"] = np.log(df_acp6_qn["acp6_std_dev_both_haps"])

df_acp6_qn["acp6_std_dev_paternal"] = df_acp6_qn.filter(like="_paternal_logrt",axis=1).std(axis="columns")
df_acp6_qn["acp6_std_dev_maternal"] = df_acp6_qn.filter(like="_maternal_logrt",axis=1).std(axis="columns")




### reset index
df_acp6_qn = df_acp6_qn.reset_index()
mean_std_dev = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_both_haps_ln"].mean()
std_dev_dev = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_both_haps_ln"].std()

# mean_std_dev_paternal = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_paternal"].mean()
# std_dev_dev_paternal = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_paternal"].std()
# mean_std_dev_maternal = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_maternal"].mean()
# std_dev_dev_maternal = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]["acp6_std_dev_maternal"].std()
threshold = mean_std_dev + 2.25 *std_dev_dev
# threshold_paternal = mean_std_dev_paternal + 2.25 * std_dev_dev_paternal
# threshold_maternal = mean_std_dev_maternal + 2.25 * std_dev_dev_maternal

#### set index
df_acp6_qn = df_acp6_qn.set_index(["chrom","start","stop"])
df_acp6_qn["acp6_vert"] =  df_acp6_qn.apply(lambda x:True if x["acp6_std_dev_both_haps_ln"]>=threshold else False,axis=1)
# df_acp6_qn["acp6_paternal_vert"] = df_acp6_qn.apply(lambda x:True if x["acp6_std_dev_paternal"]>=threshold_paternal else False,axis=1)
# df_acp6_qn["acp6_maternal_vert"] = df_acp6_qn.apply(lambda x:True if x["acp6_std_dev_maternal"]>=threshold_maternal else False,axis=1)

### reset index
df_acp6_qn = df_acp6_qn.reset_index()
df_acp6_qn["arm"] = df_acp6_qn.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

#### asynchronous per window PER sample
for sample in samples:
	df_acp6_qn[sample+"_rt_diff"] = df_acp6_qn[sample+"_paternal_logrt"] - df_acp6_qn[sample+"_maternal_logrt"]
	df_acp6_qn[sample+"_abs_rt_diff"] = abs(df_acp6_qn[sample+"_paternal_logrt"] - df_acp6_qn[sample+"_maternal_logrt"])

### calculate threshold and significant
for sample in samples:
	autosomal_tmp = df_acp6_qn[df_acp6_qn["chrom"]!="chrX"]
	mean_abs_rt_diff = autosomal_tmp[sample+"_abs_rt_diff"].mean()
	std_dev_abs_rt_diff = autosomal_tmp[sample+"_abs_rt_diff"].std()
	threshold_abs_rt_diff = mean_abs_rt_diff + 2 * std_dev_abs_rt_diff
	df_acp6_qn[sample+"_asynchronous_rt"] =  df_acp6_qn.apply(lambda row:True if row[sample+"_abs_rt_diff"]>=threshold_abs_rt_diff else False,axis=1)


# f,ax=plt.subplots(figsize=(2,2),dpi=300)
# sns.kdeplot(abs(autosomal_tmp[sample+"_paternal_logrt"] - autosomal_tmp[sample+"_maternal_logrt"]))
# ax.axvline(x=mean_abs_rt_diff + 2.5 * std_dev_abs_rt_diff,lw=0.5,linestyle="--",c="red")
# plt.show()
# plt.close()
print(df_acp6_qn)

##### make distribution of repliseq std dev. male sample
f,ax=plt.subplots(figsize=(1,1),dpi=300)
sns.kdeplot(df_acp6_qn["acp6_std_dev_both_haps"],clip=(0,20),linewidth=2)
ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
# ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="red")
# ax.axvline(x=mean_std_dev + 3 * std_dev_dev,lw=0.5,linestyle="--",c="red")
plt.savefig("acp6.vert.dist.png")
# plt.show()
plt.close()


### switchers algorithm
# switchers = []
# for index,row in df_acp6_qn.iterrows():
# 	paternal_early, paternal_late = False, False
# 	for sample in samples:
# 		if (row[sample+"_paternal_logrt"] > row[sample+"_maternal_logrt"]) and (row[sample+"_asynchronous_rt"]==True):
# 			paternal_early = True
# 		if (row[sample+"_paternal_logrt"] < row[sample+"_maternal_logrt"]) and (row[sample+"_asynchronous_rt"]==True):
# 			paternal_late = True

# 	if paternal_early and paternal_late:
# 			switchers += [True]
# 	else:
# 		switchers += [False]

# df_acp6_qn["switcher"]=switchers
# print(df_acp6_qn[(df_acp6_qn["switcher"]==True) & (df_acp6_qn["acp6_vert"]==True)])

#######

color_dict_acp6 = {
"acp6_c1":"red", 
"acp6_c2":"green",  
"acp6_c5":"blue",  
"acp6_c6":"yellow"}

color_dict_acp6_rna = {
"acp6_c1_rna":"red", 
"acp6_c2_rna":"green",  
"acp6_c5_rna":"blue",  
"acp6_c6_rna":"yellow"}
#####
df_acp6_qn = df_acp6_qn.sort_values(["chrom","start"])


df_imprinted = pd.read_csv("hg38.imprinted.genes.txt",sep="\t",
							names=["chrom","start","stop","name","score","strand"])
# acp6 should be blue
for chrom in chromosomes:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)#,height_ratios=[1,1,1])
    tmp=df_acp6_qn[df_acp6_qn["chrom"]==chrom]
    tmp["color"]=["blue" if x>=threshold else "lightgray" for x in tmp["acp6_std_dev_both_haps_ln"]]
    # tmp["color_paternal"]=["blue" if x>=threshold_paternal else "lightgray" for x in tmp["acp6_std_dev_paternal"]]
    # tmp["color_maternal"]=["blue" if x>=threshold_maternal else "lightgray" for x in tmp["acp6_std_dev_maternal"]]


    #### highlighting imprinted genes
    # for index,row in df_imprinted[(df_imprinted["chrom"]==chrom)].iterrows():
    #     rect=Rectangle((row["start"]-5000, 0), width=row["stop"]-row["start"]+10000, height=10,
    #              facecolor="gray",alpha=0.5,fill=True) ## red if hap1 early, blue if hap2 early
    #     ax[0].add_patch(rect)
    # for index,row in df_imprinted[(df_imprinted["chrom"]==chrom)].iterrows():
    #     rect=Rectangle((row["start"]-5000, 0), width=row["stop"]-row["start"]+10000, height=10,
    #          facecolor="gray",alpha=0.5,fill=True) ## red if hap1 early, blue if hap2 early
    #     ax[1].add_patch(rect)
    # for index,row in df_imprinted[(df_imprinted["chrom"]==chrom)].iterrows():
    #     rect=Rectangle((row["start"]-5000, 0), width=row["stop"]-row["start"]+10000, height=10,
    #              facecolor="gray",alpha=0.5,fill=True) ## red if hap1 early, blue if hap2 early
    #     ax[2].add_patch(rect)

    ax.scatter(tmp["start"],
		tmp["acp6_std_dev_both_haps"],c=tmp["color"],
		s=15,edgecolor="black",lw=0.2)
    # ax[1].scatter(tmp["start"],
	# 	tmp["acp6_std_dev_paternal"],c=tmp["color_paternal"],
	# 	s=15,edgecolor="black",lw=0.2)
    # ax[2].scatter(tmp["start"],
	# 	tmp["acp6_std_dev_maternal"],c=tmp["color_maternal"],
	# 	s=15,edgecolor="black",lw=0.2)    


        # ax[1].add_patch(rect)
        # ax[2].add_patch(rect)

    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length[chrom],25)) 
    ax.set_ylim([0,2.1])
    ax.set_yticks([0,.5,1,1.5,2])

    # ax[0].set_xlim([0,chromosome_length[chrom]])
    # ax[0].set_xticks(np.linspace(0,chromosome_length[chrom],25)) 
    # ax[0].set_ylim([0,2.5])
    # ax[1].set_xlim([0,chromosome_length[chrom]])
    # ax[1].set_xticks(np.linspace(0,chromosome_length[chrom],25)) 
    # ax[1].set_ylim([0,2.5])
    # ax[2].set_xlim([0,chromosome_length[chrom]])
    # ax[2].set_xticks(np.linspace(0,chromosome_length[chrom],25)) 
    # ax[2].set_ylim([0,2.5])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig("acp6.std.rt."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

df_acp6_qn.to_csv("acp6.rt.txt",sep="\t",index=False)
df_acp6_qn.to_csv("acp6.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)
## decide whether to keep coding and noncoding or not...locsxxx, SNERPAs, lncrnas, mirnas, etc
os.system("sort -k1,1 -k2,2n acp6.rt.bed | grep True | awk '$1!=\"chrX\"{print $0}' > acp6.rt.vert.sorted.bed")
os.system("bedtools map -a acp6.rt.vert.sorted.bed -b ucsc.refseq.hg38.txn.whole.gene.sorted.bed -o distinct -c 4 > acp6.rt.vert.intersect.coding.bed ")
os.system("awk '{print $16}'  acp6.rt.vert.intersect.coding.bed | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > acp6.vert.genes.txt")
## output chrom, start, stop, acp6_vert
df_acp6_qn[df_acp6_qn["acp6_vert"]==True].loc[:,["chrom","start","stop","acp6_vert"]].to_csv("acp6.rt.vert.sorted.4col.bed",sep="\t",header=False,index=None)
os.system("bedtools merge -i acp6.rt.vert.sorted.4col.bed > acp6.rt.vert.sorted.4col.merged.bed")



##### plot switchers
# tmp = df_acp6_qn[(df_acp6_qn["switcher"]==True) & (df_acp6_qn["acp6_vert"]==True)]
# tmp["unique_pos"] = [row["chrom"]+":"+str(row["start"])+"-"+str(row["stop"]) for index,row in tmp.iterrows()]
# f,ax=plt.subplots(1,1,figsize=(12,2))
# for sample in samples:
# 	ax.scatter(tmp["unique_pos"],tmp[sample+"_paternal_logrt"],
# 		c=color_dict_acp6[sample],s=15,edgecolor="black",lw=0.1,zorder=3)
# 	ax.scatter(tmp["unique_pos"],tmp[sample+"_maternal_logrt"],
# 		c=color_dict_acp6[sample],s=15,edgecolor="black",lw=0.1,zorder=3)

# for index,row in tmp.iterrows():
# 	ax.axvline(x=row["unique_pos"],linestyle="--",lw=0.4,c="black")
# ax.axhline(y=0,linestyle="--")

# plt.show()



# zooms zooms zooms
regions = df_acp6_qn[df_acp6_qn["acp6_vert"]==True]
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

            ax[1].plot(df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"]) & (df_acp6_qn["arm"]==k)]["start"],
                        df_acp6_qn[(df_acp6_qn["chrom"]==row["chrom"]) & (df_acp6_qn["arm"]==k)]["acp6_std_dev_both_haps"],
                        c="black",
                        lw=0.5)
        ax[1].set_ylim([0,1.5])
        ax[1].set_yticks([0,1.5])
        ax[1].set_xlim([row["start"]-1200000,row["stop"]+1200000]) # needs to be matched with above xlim
        ax[1].set_xticks(np.linspace(row["start"]-1200000,row["stop"]+1200000,5)) 
    plt.savefig("acp6.rt."+str(row["chrom"])+"-"+str(row["start"])+".png",
        dpi=400,transparent=False)
    plt.close()

exit()

# fig,ax=plt.subplots(1,1)
# ax.scatter(df2["count"],df_acp6.sum(axis=1),s=10)
# ax.set_xlim([0,500])
# ax.set_ylim([0,3500])
# plt.show()

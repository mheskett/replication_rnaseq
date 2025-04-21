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
import statsmodels.stats.multitest as mt
import glob

def get_switchers(df):
    # unique_genes = np.unique(df["name"])
    unique_genes = np.unique(df["name_strand"])
    switchers = [] 
    df_significant_rows = df[(df["fdr_pval"]<=0.05) & (abs(df["aei"])>=0.2)]
    ### switchers algorithm
    for i in range(len(unique_genes)):
        samples = df_significant_rows[df_significant_rows["name_strand"]==unique_genes[i]]
        if len(samples)<=1:
            continue

        hap1_skew, hap2_skew = False, False
        for index, row in samples.iterrows():
            if (row["aei"] >= 0.2):
                hap1_skew = True
            if (row["aei"] <= -0.2):
                hap2_skew = True

        if hap1_skew and hap2_skew:
            switchers += [samples]

    switchers = pd.concat(switchers)

    return switchers

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
                                alpha=0.05,
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

#####
#RNA
gm_samples = glob.glob("./gm*_c*_rna*rmv.blck.comprehensive.gene.counts.bed")

dfs=[]
dfs_with_x=[]
for x in gm_samples:
    samp = os.path.basename(x).split(".")[0]
    df = pd.read_csv(x,sep="\t",
                            names= ["chrom","start","stop","name","score","strand","paternal_counts","maternal_counts","strand_reads"],
                            dtype = {"chrom":str,"strand_reads":str,"start":int,"stop":int,"paternal_counts":int,"maternal_counts":int})
    df["total_reads"] = df["paternal_counts"] + df["maternal_counts"]
    df["aei"] = df.apply(helper_func, axis = 1)
    df["sample"] = samp
    # get rid of total
    df = df[df["strand_reads"]!= "total"]
    df = df[df["strand_reads"]!= "antisense"]
    df["informative_reads_per_kb"] = df["total_reads"] / ((df["stop"] - df["start"])  / 1000)
    df = df[df["total_reads"] >= 8]
    df_with_x = df.copy()
    df = df[df["chrom"]!="chrX"]
    # filter out chrX to lower P-values on autosomes
    add_binom_pval(df)
    add_binom_pval(df_with_x)

    dfs_with_x +=[df_with_x]
    dfs += [df]
df_gm = pd.concat(dfs)
df_gm_with_x = pd.concat(dfs_with_x)
df_gm=df_gm.sort_values(["chrom","start","sample"])
df_gm_with_x=df_gm_with_x.sort_values(["chrom","start","sample"])

df_gm["name_strand"] = df_gm["name"]+"_"+df_gm["strand_reads"]

df_gm["70_percent_aei"] = abs(df_gm["aei"]) >= 0.20
df_gm["80_percent_aei"] = abs(df_gm["aei"]) >= 0.30
df_gm["90_percent_aei"] = abs(df_gm["aei"]) >= 0.40
df_gm["95_percent_aei"] = abs(df_gm["aei"]) >= 0.45

switchers=get_switchers(df_gm)
# print(switchers["name"].unique())

pd.Series(switchers["name"].unique()).to_csv("gm.rna.switchers.txt",sep="\t",index=False,header=False)

print(switchers)
df_gm.to_csv("gm.as.gene.counts.txt",sep="\t",index=None)
df_gm_with_x.to_csv("gm.with.x.as.gene.counts.txt",sep="\t",index=None)

##
color_dict_gm = {'gm12878_clone5_rnaAligned':"red", 'gm12878_clone4_rnaAligned' :"blue"}


#### make scatter plot dots including faded dots
tmp = df_gm[df_gm["name_strand"].isin(switchers["name_strand"])]
tmp["color"]= [color_dict_gm[x] for x in tmp["sample"]]
tmp["alpha"] = tmp.apply(lambda row: 1 if (row["fdr_reject"]==True) and (abs(row["aei"])>=0.2) else 0.1, axis=1)
tmp["unique_pos"] = [row["chrom"]+":"+row["name_strand"] for index,row in tmp.iterrows()]
f,ax=plt.subplots(1,1,figsize=(3,1))
ax.scatter(tmp["unique_pos"],tmp["aei"],c=tmp["color"],s=15,edgecolor="black",lw=0.1,zorder=3,alpha=tmp["alpha"])
for index,row in tmp.drop_duplicates(["unique_pos"]).iterrows():
    ax.axvline(x=row["unique_pos"],linestyle="--",lw=0.4,c="black")
plt.xticks(rotation = 280,fontsize=3)
ax.margins(x=.015,y=0)
ax.set_ylim([-0.53,.53])
ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
ax.set_yticks([-0.5,-.25,0,.25,.5])
# f.subplots_adjust(bottom=2)
# f.subplots_adjust(left=0.09, bottom=0.2, right=0.1, top=0.0)
# f.subplots_adjust(right=0.7) 
plt.savefig("gm.rna.switchers.png",
        dpi=300,transparent=True,bbox_inches="tight",pad_inches=0)
plt.close()


exit()









#sort for plotting
# df_acp6=df_acp6.sort_values(["chrom","start"])

# plot counts per gene
print(df_acp6)
plt.hist(np.log2(df_acp6["total_reads"]+1),bins=100)
# plt.show()
plt.close()
df_acp6['skew_std_dev'] = df_acp6.groupby('name')['skew'].transform('std')
# exit()
skew_std_dev_df = df_acp6.loc[:,["chrom","start","stop","name","skew_std_dev"]].drop_duplicates()
mean_std_dev = skew_std_dev_df["skew_std_dev"].mean()
std_dev_dev = skew_std_dev_df["skew_std_dev"].std()
threshold = mean_std_dev + 2.5 *std_dev_dev
df_acp6["acp6_vee"] =  df_acp6.apply(lambda x:True if x["skew_std_dev"]>=threshold else False,axis=1)
skew_std_dev_df["acp6_vee"] = skew_std_dev_df.apply(lambda x:True if x["skew_std_dev"]>=threshold else False,axis=1)
print(skew_std_dev_df)

skew_std_dev_df.sort_values("skew_std_dev",ascending=False).to_csv("acp6.rna.vee.txt",sep="\t")
f,ax=plt.subplots(figsize=(3,1),dpi=300)
plt.hist(skew_std_dev_df["skew_std_dev"],bins=100,color="black")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="black")

# plt.show()
plt.close()


skew_std_dev_df = skew_std_dev_df.sort_values(["chrom","start"])

#acp6 should be green
for chrom in chromosomes:
    tmp = skew_std_dev_df[skew_std_dev_df["chrom"]==chrom]
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp["color"]=["blue" if x>=threshold else "lightgray" for x in tmp["skew_std_dev"]]
    
    ax.scatter(tmp["start"],
		tmp["skew_std_dev"],c=tmp["color"],
		s=15,edgecolor="black",lw=0.2)

    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length[chrom],25)) 

    ax.set_ylim([0,0.6])
    plt.savefig("acp6.rna.skew.std."+chrom+".png",
        dpi=400,transparent=False, bbox_inches='tight', pad_inches = 0)
    plt.close()


exit()





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



# fig,ax=plt.subplots(1,1)
# ax.scatter(df2["count"],df_acp6.sum(axis=1),s=10)
# ax.set_xlim([0,500])
# ax.set_ylim([0,3500])
# plt.show()

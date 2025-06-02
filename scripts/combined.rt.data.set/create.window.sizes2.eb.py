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

# def get_arms_nochr(cytoband):
#     ## given a data frame with genome elements, add the arm information to a new column
#     arm_dict = {}
#     for i in range(len(chromosomes)):
#         # should be (p end, q end)
#         arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
#         cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
#     return arm_dict

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



chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X"]


arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband= pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
arm_dict=get_arms(cytoband)
print(arm_dict)

chromosomes_hg38 = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

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

#hg19
chromosome_length = {"1":249250621,
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

ratios={"1":249250621/249250621,
"2":243199373/249250621,
"3":198022430/249250621,
"4":191154276/249250621,
"5":180915260/249250621,
"6":171115067/249250621,
"7":159138663/249250621,
"8":146364022/249250621,
"9":141213431/249250621,
"10":135534747/249250621,
"11":135006516/249250621,
"12":133851895/249250621,
"13":115169878/249250621,
"14":107349540/249250621,
"15":102531392/249250621,
"16":90354753/249250621,
"17":81195210/249250621,
"18":78077248/249250621,
"19":59128983/249250621,
"20":63025520/249250621,
"21":48129895/249250621,
"22":51304566/249250621,
"X":155270560/249250621}

#####
#### Replication timing Replication timing Replication timing Replication timing Replication timing 
eb_samples = glob.glob("eb*sorted.markdup.allele.counts.haplotype.resolved.counts.rmv.blck.cpms.windows.w250.s250.bed")

print(eb_samples)
dfs=[]
for x in eb_samples:
    samp = os.path.basename(x).split(".")[0]
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", samp+"_cpm_p_counts", samp+"_cpm_m_counts"])
    # tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 
    # tmp = tmp[(tmp[samp+"_cpm_p_counts"]!=".") & (tmp[samp+"_cpm_m_counts"]!=".")]
    tmp = tmp.replace(".", 0)
    tmp[samp+"_cpm_p_counts"] = tmp[samp+"_cpm_p_counts"].astype(float)
    tmp[samp+"_cpm_m_counts"] = tmp[samp+"_cpm_m_counts"].astype(float)
    # tmp = tmp[tmp[samp+"_cpm_p_counts"]>=1]
    # tmp = tmp[tmp[samp+"_cpm_m_counts"]>=1]
    tmp=tmp.set_index(["chrom","start","stop"])
    dfs+=[tmp]
df_eb = pd.concat(dfs,axis=1).sort_values(["chrom","start"])
samples=list(set(["_".join(os.path.basename(x).split(".")[0].split("_")[0:3]) for x in eb_samples]))

#sort for plotting
df_eb=df_eb.sort_values(["chrom","start"])

## replace NaNs with 0s
df_eb = df_eb.replace(".", 0)

## remove rows with all 0's
df_eb = df_eb[df_eb.sum(axis=1)>0]

# plot counts per window
print(df_eb)
# sns.kdeplot(df_eb.sum(axis=1),clip=(0,5000))
# plt.show()

## according to the plot, remove all windows with less than ~150 LSM counts
df_eb = df_eb[df_eb.sum(axis=1)>125]

print(samples)
## now do E/L within the windows
for sample in samples:
    df_eb[sample+"_paternal_logrt"] = np.log2((df_eb[sample+"_early_rt_cpm_p_counts"]+1) / (df_eb[sample+"_late_rt_cpm_p_counts"]+1 ))
    df_eb[sample+"_maternal_logrt"] = np.log2((df_eb[sample+"_early_rt_cpm_m_counts"]+1) / (df_eb[sample+"_late_rt_cpm_m_counts"]+1 ))

df_eb_repli = df_eb.filter(like="logrt")
print(df_eb_repli)

### process and plot
df_eb_qn = quantile_normalize(df_eb.filter(regex="logrt"))
df_eb_qn["eb_std_dev_both_haps"] = df_eb_qn.filter(like="eb",axis=1).std(axis="columns")
df_eb_qn["eb_std_dev_both_haps_ln"] = np.log(df_eb_qn["eb_std_dev_both_haps"])
df_eb_qn = df_eb_qn.reset_index()
mean_std_dev = df_eb_qn[df_eb_qn["chrom"]!="X"]["eb_std_dev_both_haps_ln"].mean()
std_dev_dev = df_eb_qn[df_eb_qn["chrom"]!="X"]["eb_std_dev_both_haps_ln"].std()
threshold = mean_std_dev + 2.25 *std_dev_dev
df_eb_qn = df_eb_qn.set_index(["chrom","start","stop"])
df_eb_qn["eb_vert"] =  df_eb_qn.apply(lambda x:True if x["eb_std_dev_both_haps_ln"]>=threshold else False,axis=1)

df_eb_qn = df_eb_qn.reset_index()
df_eb_qn["arm"] = df_eb_qn.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

print(threshold)
print("eb 225",df_eb_qn)
##### make distribution of repliseq std dev. male sample
f,ax=plt.subplots(figsize=(2,2),dpi=300)
sns.kdeplot(df_eb_qn["eb_std_dev_both_haps_ln"],clip=(0,20),linewidth=2)
ax.axvline(x=mean_std_dev,lw=0.5,linestyle="--",c="black")
# ax.axvline(x=mean_std_dev + 2 * std_dev_dev,lw=0.5,linestyle="--",c="red")
ax.axvline(x=mean_std_dev + 2.5 * std_dev_dev,lw=0.5,linestyle="--",c="red")
# ax.axvline(x=mean_std_dev + 3 * std_dev_dev,lw=0.5,linestyle="--",c="red")
plt.savefig("eb.vert.dist.png")
# plt.show()
plt.close()

color_dict_eb = {'eb3_2_clone15':"red", 
'eb3_2_clone4':"cyan", 
'eb3_2_clone13':"yellow", 
'eb3_2_clone10':"green", 
'eb3_2_clone2':"blue", 
'eb3_2_clone3':"purple"}
#####
df_eb_qn = df_eb_qn.sort_values(["chrom","start"])

### hg19 to hg38 liftover for plotting
### make a liftover python function that inputs a df and outputs a df?
df_eb_qn.to_csv("eb.rt.hg19.bed",sep="\t",header=None,na_rep="NaN",index=False)
## liftover chains require chr prefix, but the official hg19 has no chr prefix
os.system("awk 'OFS=\"\t\"{print \"chr\"$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' eb.rt.hg19.bed > eb.rt.hg19.chr.bed")
os.system("/Users/michaelheskett/replication_rnaseq/scripts/liftover_files/liftOver eb.rt.hg19.chr.bed /Users/michaelheskett/replication_rnaseq/scripts/liftover_files/hg19ToHg38.over.chain eb.rt.hg38.lifted.bed eb.rt.hg38.unmapped.bed -bedPlus=3")
# now remove the chrom from lifted
os.system("sed 's/chr//g' eb.rt.hg38.lifted.bed > eb.rt.hg38.lifted.nochr.bed")
### maybe make sure this works?!?!?!
df_eb_qn_hg38=pd.read_csv("eb.rt.hg38.lifted.bed",sep="\t",names=df_eb_qn.columns)
print(df_eb_qn_hg38)
df_eb_qn_hg38.to_csv("eb.rt.hg38.lifted.txt",sep="\t",na_rep="NaN",index=False)

## EB should be red dots
## should be using hg38 coords
for chrom in chromosomes_hg38:
    plt.rc('xtick', labelsize=5) 
    plt.rc('ytick', labelsize=3) 
    f, ax = plt.subplots(1,1,figsize=(15,2),dpi=300)
    tmp=df_eb_qn_hg38[df_eb_qn_hg38["chrom"]==chrom]
    tmp["color"]=["red" if x>=threshold else "lightgray" for x in tmp["eb_std_dev_both_haps_ln"]]
    
    ax.scatter(tmp["start"],tmp["eb_std_dev_both_haps"],c=tmp["color"],
		s=15,edgecolor="black",lw=0.2)

    ax.set_xlim([0,chromosome_length_hg38[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length_hg38[chrom],25)) 

    ax.set_ylim([0,2.1])
    ax.set_yticks([0,.5,1,1.5,2])
    plt.savefig("eb.std.rt."+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
    plt.close()

df_eb_qn.to_csv("eb.rt.txt",sep="\t",index=False)
df_eb_qn.to_csv("eb.rt.bed",sep="\t",header=None,na_rep="NaN",index=False)

## just translate the whole thing to hg38
## decide whether to keep coding and noncoding or not...locsxxx, SNERPAs, lncrnas, mirnas, etc
## get rid of X genes for female sample.
os.system("sort -k1,1 -k2,2n eb.rt.bed | grep True | awk '$1!=\"X\"{print $0}' > eb.rt.vert.sorted.bed")
os.system("bedtools map -a eb.rt.vert.sorted.bed -b ucsc.known.gene.hg19.txn.start.stop.bed.cds.only.first.isoform.nochr.sorted.bed -o distinct -c 4 > eb.rt.vert.intersect.coding.bed ")
os.system("awk '{print $20}'  eb.rt.vert.intersect.coding.bed | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > eb.vert.genes.txt")


# do with hg38
os.system("sort -k1,1 -k2,2n eb.rt.hg38.lifted.bed | grep True | awk '$1!=\"chrX\"{print $0}' > eb.rt.vert.hg38.lifted.sorted.bed")
os.system("bedtools map -a eb.rt.vert.hg38.lifted.sorted.bed -b ucsc.refseq.hg38.txn.whole.gene.sorted.bed -o distinct -c 4 > eb.hg38.rt.vert.intersect.coding.bed ")
os.system("awk '{print $19}'  eb.hg38.rt.vert.intersect.coding.bed  | awk '{$1=$1} 1' FS=, OFS='\\n'| sort | uniq | grep -v ^LOC | grep -v LINC | grep -v MIR | grep -v SNORD > eb.hg38.vert.genes.txt")

df_eb_qn_hg38[df_eb_qn_hg38["eb_vert"]==True].loc[:,["chrom","start","stop","eb_vert"]].to_csv("eb.rt.vert.hg38.lifted.sorted.4col.bed",sep="\t",header=False,index=None)
os.system("bedtools merge -d 250001 -i eb.rt.vert.hg38.lifted.sorted.4col.bed | awk '$3-$2>250000{print $0}' | awk 'OFS=\"\\t\"{print $1,$2,$3}' > eb.rt.vert.hg38.lifted.sorted.4col.merged.bed")

## w



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

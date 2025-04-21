import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
## goes through each sample and generates window sizes based on snp density in the phasing

# step 1. plot reads per SNP and remove any with super low coverage across samples
# step 2. plot SNPs per 1kb, snps per 10kb, snps per 100kb within the haplotype
# step 3. library size normalization
# step 4. algorithm to create windows with sufficient reads and snps

# haplotye = pd.read_csv("acp6_joint.fb30.all.chrom.no.missing.wh.phased.final.hets.bed",
# 						sep="\t",names=["chrom","start","stop","hap1","hap2"])

cytoband = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
dfs=[]
acp6_samples = glob.glob("acp6*hg38*counts.bed")
for x in acp6_samples:
    samp = os.path.basename(x).split(".")[0]
    # print(samp)
    tmp=pd.read_csv(x,sep="\t",names=["chrom" ,"start", "stop", "paternal","maternal",samp+"_p_counts", samp+"_m_counts"])
    tmp=tmp[tmp["chrom"]!="X"]###removing chrom x
    tmp=tmp[tmp["chrom"]!="Y"]###removing chrom x 
    tmp=tmp.drop(["paternal","maternal"],axis=1)
    tmp = tmp[(tmp[samp+"_p_counts"]!=".") & (tmp[samp+"_m_counts"]!=".")]
    tmp[samp+"_p_counts"] = tmp[samp+"_p_counts"].astype(float)
    tmp[samp+"_m_counts"] = tmp[samp+"_m_counts"].astype(float)
    # tmp = tmp[tmp[samp+"_cpm_p_counts"]>=1]
    # tmp = tmp[tmp[samp+"_cpm_m_counts"]>=1]
    # tmp=tmp.set_index(["chrom","start","stop"])
    dfs+=[tmp]
df = dfs[0]
samples = [os.path.basename(x).split(".")[0] for x in acp6_samples]

## merge the DFs and set index
for d in dfs[1:]:
    df = df.merge(d, how="outer",left_on=["chrom","start","stop"],right_on=["chrom","start","stop"])
df=df.set_index(["chrom","start","stop"])

## replace NaNs with 0s
df = df.replace(np.nan, 0)

## remove rows with all 0's
df = df[df.sum(axis=1)>0]

## do LSM
for sample in samples:
	total_reads = (df[sample+"_p_counts"]+df[sample+"_m_counts"]).sum()
	df[sample+"_p_counts_lsm"] = (df[sample+"_p_counts"]/total_reads)*10**6
	df[sample+"_m_counts_lsm"] = (df[sample+"_m_counts"]/total_reads)*10**6

## make the lsm dataframe separated from the original df
df_lsm = df.filter(regex="lsm",axis=1)

# create column that has read sums at each position
# df_lsm["position_sum"] = df_lsm.sum(axis=1)

# call bedtools makewindows and do map count
os.system("bedtools makewindows -g hg38.fa.fai -w 200000 -s 200000 > hg38.windows.w200.s200.bed")

print(df_lsm)

# sns.kdeplot(df_lsm.sum(axis=1),clip=(0,60))
# plt.show()
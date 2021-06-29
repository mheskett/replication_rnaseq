import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
import seaborn as sns
import scipy.stats
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
import glob

l1_df = pd.read_csv("/Users/mike/replication_rnaseq/scripts/ucsc.genes.cds.l1.fraction.final.bed",
	sep="\t",header=None,index_col=None,names=["chrom","start","stop","name","score","strand", "l1_fraction"],
	dtype={"chrom":str,"start":int,"stop":int,"score":float,"strand":str,"l1_fraction":float})

l1_df = l1_df.drop_duplicates(subset=["chrom","start","stop"])
l1_auto = l1_df[l1_df["chrom"]!="X"]
l1_x = l1_df[l1_df["chrom"]=="X"]
# sns.kdeplot(data=l1_auto["l1_fraction"],cut=0,lw=4)
# sns.kdeplot(data=l1_x["l1_fraction"],cut=0,lw=4)

# plt.show()




#["/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed",
			#"/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.10Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"]


# # all_files = list(glob.glob("/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression" + "/*vlinc.discovery.all.bed"))
# all_files = ["/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.2.all.vlincs.and.repliseq.bed",
# 			"/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.10.all.vlincs.and.repliseq.bed",]


all_files = ["/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/gm12878.4.all.vlincs.and.repliseq.bed",
"/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/gm12878.5.all.vlincs.and.repliseq.bed"]
# # all_files.sort()
filenames=[os.path.basename(x)[0:15] for x in all_files]
li=[]
number=[]
skew=[]
f,x = plt.subplots(figsize=(2,2))
for i in range(len(all_files)):
	df = pd.read_csv(all_files[i],sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew","hap1_early","hap2_early","hap1_late","hap2_late"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"reject":str,"hap1_counts":int,"hap2_counts":int,
						"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})
	df["sample"] = filenames[i]
	tmp = df.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
	tmp = tmp.astype(int)
	df.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
	###########
	df.loc[:,"logr_hap1"] = df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
	df.loc[:,"logr_hap2"] = df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
	df.loc[:,"logr_diff"] = abs(df["logr_hap1"] - df["logr_hap2"])
	li.append(df)
	number += [[filenames[i],len(df)]]
	skew += [df["skew"]]
	# df["skew"].plot.kde(ind=np.linspace(-.5,.5,1000),label=filenames[i])
	sns.kdeplot(data=df[df["chrom"]!="X"]["skew"],cut=0,label=filenames[i],lw=4,c="mediumblue" if (filenames[i]=="bouha.2.all.vli" or filenames[i]=="gm12878.4.all.v") else "orange")
plt.xlim([-0.5,0.5])
plt.xticks([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,.2,.3,.4,0.5])
# plt.legend(bbox_to_anchor=(1, 0.5))
plt.savefig("skew_histo.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
##################
# f,x = plt.subplots(figsize=(2,2))
# sns.kdeplot(data=df[df["skew"]],cut=0,lw=4)
# sns.kdeplot(data=df["skew"],cut=0,lw=4)
# plt.xlim([-0.5,0.5])
# plt.xticks([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,.2,.3,.4,0.5])
### plot rpkm vs skew

# f,ax=plt.subplots(figsize=(10,4))
# plt.scatter(df[df["chrom"]!="X"]["rpkm"],abs(df[df["chrom"]!="X"]["skew"]),s=10,lw=0.2,edgecolor="black")
# plt.xlim([0,250])
# plt.show()



##############
# try violin plots
long_df = pd.concat(li,axis=0)
# long_df = long_df[(long_df['sample']=="gm12878.4x.hg19") |
# 				 (long_df['sample']=="gm12878.5x.hg19") | 
# 				 (long_df['sample']=="bouha.trim.2Ali") |
# 				 (long_df['sample']=="bouha.trim.10Al")]
print(long_df)
# sns.violinplot(x=long_df["sample"],y=long_df["skew"],color="skyblue")
# plt.ylim([-0.5,.5])
# plt.show()
# plt.close()
#############
# try violin plots
# long_df_skewed = long_df[long_df["reject"]=="True"]
# sns.violinplot(x=long_df_skewed["sample"],y=long_df_skewed["l1_fraction"],color="skyblue")
# plt.ylim([0,1])
# plt.show()

# sns.violinplot(x=long_df_skewed["sample"],y=np.log2(long_df_skewed["rpkm"]),color="skyblue")
# # plt.ylim([])
# plt.show()


# sns.violinplot(x=long_df_skewed["sample"],y=np.log10(long_df_skewed["stop"] - long_df_skewed["start"]),color="skyblue")
# # plt.ylim()
# plt.show()

sns.kdeplot(data=l1_df["l1_fraction"],cut=0,lw=4)
# sns.kdeplot(data=l1_x["l1_fraction"],cut=0,lw=4)
sns.kdeplot(data=long_df[(long_df["reject"]=="True") & (long_df["pval"]<=0.00000001) & (long_df["chrom"]!="X")]["l1_fraction"],cut=0,lw=4)
sns.kdeplot(data=long_df[(long_df["reject"]=="False") & (long_df["chrom"]!="X")]["l1_fraction"],cut=0,lw=4)
plt.xlim([0,1])
plt.xticks([0,.05,.1,.15,.2,.25,.5,.8,1])
plt.legend()
plt.show()
plt.close()

# plt.show()
f,ax = plt.subplots(figsize=(2.5,2))
sns.kdeplot(data=l1_df["l1_fraction"],cut=0,lw=4,c="mediumblue",ax=ax)
# sns.kdeplot(data=l1_x["l1_fraction"],cut=0,lw=4)
sns.kdeplot(data=long_df["l1_fraction"],cut=0,lw=4,c="red",ax=ax)
# sns.kdeplot(data=long_df[long_df["reject"]=="False"]["l1_fraction"],cut=0,lw=4,c="orange")
# sns.kdeplot(data=long_df[long_df["pval"]<=0.00001]["l1_fraction"],cut=0,lw=4,c="orange")
plt.xlim([0,0.5])
plt.xticks(list(np.arange(0,0.6,0.1)),fontsize=15)
plt.savefig("l1_hist.png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
###############
long_df["log_pval"] = -np.log10(long_df["pval"])
long_df["abs_skew"] = abs(long_df["skew"])
f,ax = plt.subplots(figsize=(4,4))
ax.scatter(abs(long_df[(long_df["reject"]=="True") & (long_df["pval"]<=0.001)]["skew"]), 
	long_df[(long_df["reject"]=="True")& (long_df["pval"]<=0.001)]["l1_fraction"],s=30,lw=0.2,edgecolor="black")
# sns.clustermap(long_df.loc[:,["rpkm","l1_fraction","abs_skew","log_pval","total_reads"]].corr(method="pearson"),vmin=0,vmax=1)
plt.show()




plt.subplots(figsize=(2,2))
for i in range(len(li)):
	tmp=li[i]

	sns.kdeplot(data=tmp[tmp["reject"]=="True"]["logr_diff"],c="darkblue",cut=0,label=filenames[i])
plt.xlim([0,6])
plt.legend()
plt.show()
plt.close()


plt.subplots(figsize=(2,2))
for i in range(len(li)):
	tmp=li[i]
	sns.kdeplot(data=tmp[tmp["reject"]=="True"]["l1_fraction"],cut=0,label=filenames[i])
	sns.kdeplot(data=tmp[tmp["reject"]=="False"]["l1_fraction"],cut=0,label=filenames[i])
plt.xlim([0,1])
plt.legend()
plt.show()
plt.close()
# ################
# plt.subplots(figsize=(12,12))
# for i in range(len(li)):
# 	tmp=li[i]

# 	sns.kdeplot(data=tmp[tmp["reject"]=="True"]["skew"],cut=0,label=filenames[i])
# plt.xlim([0,1])
# plt.legend()
# plt.show()
# plt.close()
# ########

# plt.subplots(figsize=(12,12))
# for i in range(len(li)):
# 	tmp=li[i]

# 	sns.kdeplot(data=tmp[tmp["reject"]=="True"]["rpkm"],cut=0,label=filenames[i])
# plt.xlim([0,300])
# plt.legend()
# plt.show()
# plt.close()
# ######plt.subplots(figsize=(12,12))
# plt.subplots(figsize=(12,12))
# for i in range(len(li)):
# 	tmp=li[i]

# 	sns.kdeplot(data=np.log10(tmp[tmp["reject"]=="True"]["stop"] - tmp[tmp["reject"]=="True"]["start"]),cut=0,label=filenames[i])
# plt.legend()
# plt.show()
# plt.close()



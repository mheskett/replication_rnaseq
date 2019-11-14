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
from sklearn import metrics
from sklearn.cluster import DBSCAN
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn import manifold, datasets
###

###


all_vlincs = pybedtools.BedTool("gm12878.rep2.hg19Aligned.outgm12878.rep2.expressed.vlincs.all.bed").sort(g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
line_age = pybedtools.BedTool("/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.age.hg19.bed").sort(g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")

## now do a bedtool map with lines and -c mean
df = all_vlincs.map(line_age, c=7 ,o="mean", g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")\
				.to_dataframe(names=["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
	"binom_pval","fdrpval","reject","total_reads","skew","average_line_age"],
	dtype={"chrom":str,"start":int,"stop":int,"name":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
	"binom_pval":float,"fdrpval":float,"reject":str,"total_reads":int,"skew":float,"average_line_age":str} )

df=df[df["average_line_age"]!="."]
df.loc[:,"average_line_age"] = df["average_line_age"].astype(float)
print(df)

# df = pd.read_csv("gm12878.rep1.hg19Aligned.outgm12878.rep1.expressed.vlincs.all.bed",
# 	sep="\t",header=None,
# 	names=["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
# 	"binom_pval","fdrpval","reject","total_reads","skew"],
# 	dtype={"chrom":str,"start":int,"stop":int,"name":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
# 	"binom_pval":float,"fdrpval":float,"reject":str,"total_reads":int,"skew":float})


df.loc[:,"length"] = df["stop"] - df["start"]
df.loc[:,"length_zscore"] = scipy.stats.mstats.zscore(df["length"])
df.loc[:,"reads/kb_zscore"] = scipy.stats.mstats.zscore(df["rpkm"])
df.loc[:,"fraction_l1_zscore"] = scipy.stats.mstats.zscore(df["l1_fraction"])
df.loc[:,"skew_zscore"] = scipy.stats.mstats.zscore(abs(df["skew"]))
df.loc[:,"average_line_age_zscore"] = scipy.stats.mstats.zscore(df["average_line_age"])

xist = df[df["name"]=="3438_gm12878.rep1.hg19Aligned.out.samtool.rmdup"].reset_index(drop=True)
df = df[df["chrom"]!="X"]
df_notskewed = df[(df["binom_pval"]>10**-6) | (abs(df["skew"])>0.4)].reset_index(drop=True)
df_skewed = df[(df["binom_pval"]<=10**-6) & (abs(df["skew"])>=0.4)].reset_index(drop=True)

# sns.set(font_scale=1.5)  # crazy big
# fig, ax = plt.subplots()
# sns.kdeplot(df_notskewed["l1_fraction"],shade=True,label="biallelic_vlincs")
# sns.kdeplot(df_skewed["l1_fraction"],shade=True,label="skewed_vlincs")
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 5])
# plt.show()
# plt.close()


# fig, ax = plt.subplots()
# sns.set(font_scale=1.5)  # crazy big

# sns.kdeplot(df_notskewed["length"],shade=True,label="biallelic_vlincs")
# sns.kdeplot(df_skewed["length"],shade=True,label="skewed_vlincs")
# ax.set_xlim([0, 1500000])
# plt.show()
# plt.close()
# sns.set(font_scale=1.5)  # sns before fig,ax gibves you seaborn colors. fig ax before sns give syou matplotlib colors
# fig, ax = plt.subplots()
# sns.kdeplot(df_notskewed["rpkm"],shade=True,label="biallelic_vlincs")
# sns.kdeplot(df_skewed["rpkm"],shade=True,label="skewed_vlincs")
# ax.set_xlim([0, 300])

# plt.show()
# plt.close()


### this is doing all of them

X = df.loc[:,["reads/kb_zscore","length_zscore","fraction_l1_zscore","skew_zscore","average_line_age_zscore"]]
reduced_data = PCA(n_components=2).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(151)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df["rpkm"],cmap="Reds",vmin=1,vmax=150)

ax = fig.add_subplot(152)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df["length"],cmap="Reds",vmin=50*10**3,vmax=500*10**3)

ax = fig.add_subplot(153)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df["l1_fraction"],cmap="Reds",vmin=0.05,vmax=0.4)

ax = fig.add_subplot(154)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df["skew"],cmap="Reds",vmin=0,vmax=1)

ax = fig.add_subplot(155)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df["average_line_age"],cmap="Reds",vmin=15,vmax=200)

plt.show()
plt.close()

print("computing tsne")
tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(151)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black", c=df["rpkm"], cmap="Reds", vmin=1, vmax=150)

ax = fig.add_subplot(152) 
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black", c=df["length"], cmap="Reds", vmin=50*10**3, vmax=500*10**3)

ax = fig.add_subplot(153)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black", c=df["l1_fraction"], cmap="Reds", vmin=0.05, vmax=0.4

ax = fig.add_subplot(154)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black", c=df["skew"],cmap="Reds", vmin=0,vmax=1)

ax = fig.add_subplot(155)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black", c=df["average_line_age"], cmap="Reds",vmin=15,vmax=200)

plt.show()
plt.close()

#### now do this for just the skewed ones and add XIST
df_skewed = df_skewed.append(xist)
print(df_skewed)
X = df_skewed.loc[:,["reads/kb_zscore","length_zscore","fraction_l1_zscore","skew_zscore","average_line_age_zscore"]]
reduced_data = PCA(n_components=2).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(151)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df_skewed["rpkm"],cmap="Reds",vmin=-2,vmax=2)

ax = fig.add_subplot(152)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df_skewed["length"],cmap="Reds",vmin=-2,vmax=2)

ax = fig.add_subplot(153)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df_skewed["l1_fraction"],cmap="Reds",vmin=-2,vmax=2)

ax = fig.add_subplot(154)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df_skewed["skew"],cmap="Reds",vmin=-2,vmax=2)

ax = fig.add_subplot(155)
ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black",c=df_skewed["average_line_age"],cmap="Reds",vmin=-2,vmax=2)


plt.show()
plt.close()

print("computing tsne")
tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(151)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df_skewed["rpkm"],cmap="Reds",vmin=-2,vmax=2)
ax.scatter(x=tsne[-1,0],y=tsne[-1,1],lw=0.1,edgecolor="black",c=xist["rpkm"],cmap="Reds",s=200,vmin=-2,vmax=2)


ax = fig.add_subplot(152)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df_skewed["l1_fraction"],cmap="Reds",vmin=-2,vmax=2)
ax.scatter(x=tsne[-1,0],y=tsne[-1,1],lw=0.1,edgecolor="black",c=xist["l1_fraction"],cmap="Reds",vmin=-2,vmax=2,s=200)


ax = fig.add_subplot(153)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df_skewed["length"],cmap="Reds",vmin=-2,vmax=2)
ax.scatter(x=tsne[-1,0],y=tsne[-1,1],lw=0.1,edgecolor="black",c=xist["length"],cmap="Reds",vmin=-2,vmax=2,s=200)


ax = fig.add_subplot(154)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df_skewed["skew"],cmap="Reds",vmin=-2,vmax=2)
ax.scatter(x=tsne[-1,0],y=tsne[-1,1],lw=0.1,edgecolor="black",c=xist["skew"],cmap="Reds",vmin=-2,vmax=2,s=200)

ax = fig.add_subplot(155)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df_skewed["average_line_age"],cmap="Reds",vmin=-2,vmax=2)
ax.scatter(x=tsne[-1,0],y=tsne[-1,1],lw=0.1,edgecolor="black",c=xist["average_line_age"],cmap="Reds",vmin=-2,vmax=2,s=200)


plt.show()
plt.close()



###make plot with x is length zscore, y is reads/kb z score, color is red for high L1, blue for low L1
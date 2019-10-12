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

df_plus = pd.read_csv("/Users/heskett/replication_rnaseq/data/gm12878.4x.hg19Aligned.out.samtool.rmdup.plus.1000.10000.50000.vlinc.discovery.bed",
				header=None,
				index_col=None,
				names=["chrom","start","stop","name","reads/kb","strand","fraction_l1"],
				sep="\t")

df_minus = pd.read_csv("/Users/heskett/replication_rnaseq/data/gm12878.4x.hg19Aligned.out.samtool.rmdup.minus.1000.10000.50000.vlinc.discovery.bed",
				header=None,
				index_col=None,
				names=["chrom","start","stop","name","reads/kb","strand","fraction_l1"],
				sep="\t")

df = pd.concat([df_plus,df_minus])
df.loc[:,"length"] = df["stop"] - df["start"]
df.loc[:,"length_zscore"] = scipy.stats.mstats.zscore(df["length"])
df.loc[:,"reads/kb_zscore"] = scipy.stats.mstats.zscore(df["reads/kb"])
df.loc[:,"fraction_l1_zscore"] = scipy.stats.mstats.zscore(df["fraction_l1"])

print(df["reads/kb"].describe())
print(df["length"].describe())
print(df["fraction_l1"].describe())

plt.hist(df["fraction_l1"])
plt.show()
plt.close()
plt.hist(df["length"])
plt.show()
plt.close()
plt.hist(df["reads/kb"])
plt.show()
plt.close()


print(df)

print(df["reads/kb_zscore"].describe())
print(df["length_zscore"].describe())
print(df["fraction_l1_zscore"].describe())
X = df.loc[:,["reads/kb_zscore","length_zscore","fraction_l1_zscore"]]
reduced_data = PCA(n_components=2).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(x=reduced_data[:,0],y=reduced_data[:,1],lw=0.1,edgecolor="black")
# plt.show()
plt.close()

print("computing tsne")
tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0).fit_transform(X.values)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df["reads/kb"],cmap="Reds",vmin=-2,vmax=4)
plt.show()
plt.close()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df["fraction_l1"],cmap="Reds",vmin=-2,vmax=4)
plt.show()
plt.close()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x=tsne[:,0],y=tsne[:,1],lw=0.1,edgecolor="black",c=df["length"],cmap="Reds",vmin=-2,vmax=4)
plt.show()
plt.close()


###make plot with x is length zscore, y is reads/kb z score, color is red for high L1, blue for low L1
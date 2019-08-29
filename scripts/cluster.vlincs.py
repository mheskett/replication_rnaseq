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
from sklearn.cluster import OPTICS, cluster_optics_dbscan

df = pd.read_table("/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed",header=None,index_col=None,
	names=["chrom","start","stop","name","score","strand"])

chromosomes = list(df["chrom"].unique())
colors = ['g.', 'r.', 'b.', 'y.', 'c.']

# for i in range(len(chromosomes)):
# 	clust = OPTICS(min_samples=5, xi=.05, min_cluster_size=2)
# 	X= df[df["chrom"]==chromosomes[i]].loc[:,["start","stop"]].values
# 	clust.fit(X)
# 	# OPTICS
# 	plt.figure(figsize=(10, 7))
# 	colors = ['g.', 'r.', 'b.', 'y.', 'c.']
# 	for klass, color in zip(range(0, 5), colors):
# 	    Xk = X[clust.labels_ == klass]
# 	    print(clust.labels_)
# 	    plt.plot(Xk[:, 0], Xk[:, 1], color, alpha=0.3)
# 	plt.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], 'k+', alpha=0.1)
# 	plt.show()
# 	break



###  try DBSCAN
for i in range(len(chromosomes)):
	X= df[df["chrom"]==chromosomes[i]].loc[:,["start","stop"]].values
	db = DBSCAN(eps=2*10**6, min_samples=5,metric="chebyshev").fit(X)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_
	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	n_noise_ = list(labels).count(-1)
	plt.figure(figsize=(10, 7))
	print('Estimated number of clusters: %d' % n_clusters_)
	print('Estimated number of noise points: %d' % n_noise_)

	# Black removed and is used for noise instead.
	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
	          for each in np.linspace(0, 1, len(unique_labels))]
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]

	    class_member_mask = (labels == k)

	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=14)

	    xy = X[class_member_mask & ~core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)

	plt.title('Estimated number of clusters: %d' % n_clusters_)
	plt.show()

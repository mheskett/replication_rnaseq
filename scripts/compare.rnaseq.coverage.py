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

path = r'' # use your path
all_files = glob.glob("/Users/mike/replication_rnaseq" + "/*10kb.counts.bed.filtered.bed")
filenames=[os.path.basename(x)[0:15] for x in all_files]
print(filenames)
li = []

for filename in all_files:
    df = pd.read_csv(filename, names=["chrom","start","stop","counts"+"_"+os.path.basename(filename)[0:15]],sep="\t",
						dtype = {"chrom":str,"start":int,"stop":int,"counts"+"_"+os.path.basename(filename)[0:15]:int}).set_index(["chrom","start","stop"])
    li.append(df)
df=pd.concat(li,axis=1)

sns.clustermap(df,row_cluster=False,vmin=0,vmax=1000,cmap="Reds")
plt.show()
plt.close()
sns.clustermap(df.corr(method="pearson"),cmap="Reds")
plt.show()
##########################
## stats on all vlincs called per sample

all_files = glob.glob("/Users/mike/replication_rnaseq/scripts" + "/*vlinc.discovery.all.bed")
filenames=[os.path.basename(x)[0:15] for x in all_files]

lengths = []
number = []
for filename in all_files:

	df = pd.read_csv(filename,sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
	lengths += [np.mean(df["stop"] - df["start"])]
	number += [len(df)]

print("median average length: ",np.median(lengths))
print("median number", np.median(number) )

###### dont use X for these

all_files = glob.glob("/Users/mike/replication_rnaseq/scripts" + "/*vlinc.discovery.skewed.bed")
filenames=[os.path.basename(x)[0:15] for x in all_files]

lengths = []
number = []
for filename in all_files:

	df = pd.read_csv(filename,sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
	df = df[df["chrom"]!="X"]
	lengths += [np.mean(df["stop"] - df["start"])]
	number += [len(df)]

print("median average length: ",np.median(lengths))
print("median number", np.median(number) )



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
import pickle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Shadow
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from sklearn.cluster import KMeans
from scipy.stats import ttest_ind
import glob


all_files_repli = glob.glob("GSM*")
cast = [x for x in all_files_repli if "CAST" in x]
c57 = [x for x in all_files_repli if "C57" in x]
pure = [x for x in all_files_repli if "pure" in x]
repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",header=0,
                        names= ["chrom","start","stop",all_files_repli[i]],
                        dtype = {"chrom":str,"start":int,"stop":int, all_files_repli[i]:float} )
    repli_li += [df_repli.set_index(["chrom","start","stop"])]

repli_df = pd.concat(repli_li,axis=1).reset_index()

repli_df["std_dev"] = repli_df.set_index(["chrom","start","stop"]).std()
mean_std_dev = repli_df["std_dev"].mean()
std_dev_dev = repli_df["std_dev"].std()
threshold = mean_std_dev + 2.5*std_dev_dev
print(threshold)
f,ax = plt.subplots(figsize=(12,2))
for i in range(len(cast)):
     ax.plot(repli_df[(repli_df["chrom"]=="chr1")]["start"],
                repli_df[(repli_df["chrom"]=="chr1")][cast[i]],lw=0.5,c="red")
for i in range(len(c57)):

    ax.plot(repli_df[(repli_df["chrom"]=="chr1")]["start"],
                repli_df[(repli_df["chrom"]=="chr1")][c57[i]],lw=0.5,c="blue")
for i in range(len(pure)):

    ax.plot(repli_df[(repli_df["chrom"]=="chr1")]["start"],
                repli_df[(repli_df["chrom"]=="chr1")][pure[i]],lw=0.5,c="green")
for index3,row3 in repli_df[(repli_df["chrom"]=="chr1") & (repli_df["std_dev"]>=threshold)].iterrows():
    rect=Rectangle((row3["start"]-250000, -5), width=row3["stop"]-row3["start"]+500000, height=10,
             facecolor="gray",alpha=1,fill=True) ## red if hap1 early, blue if hap2 early
    ax.add_patch(rect)
plt.show()
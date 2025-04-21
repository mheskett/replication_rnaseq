import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

# header line
# "source","term_name","term_id","highlighted","adjusted_p_value","negative_log10_of_adjusted_p_value","term_size","query_size","intersection_size","effective_domain_size","intersections"


files = glob.glob("gprofiler*csv")
print(files)

dfs=[]
for file in files:
	df = pd.read_csv(file,header=0,index_col=None)
	df["sample"] = file
	df["term_name_id"] = df["term_name"] + "_" + df["term_id"]
	dfs += [df]
df = pd.concat(dfs)

print(df)
# print(df[df["sample"]=="gprofiler.test2.txt"]["term_name"].value_counts().sort_values(ascending=False))
# print(df.loc[:,["sample","term_name"]].unique())
# df=df[~df["source"].isin(["TF","HP"])]
df=df[~df["sample"].isin(["gprofiler_mouse.csv"])]
dat = df.pivot(columns="sample",index="term_name_id",values="negative_log10_of_adjusted_p_value")#.reset_index()
# dat = dat.fillna(0)
dat['non_missing_count'] = dat.notna().sum(axis=1)
dat_sorted = dat.sort_values(by='non_missing_count', ascending=False)
dat_sorted = dat_sorted[dat_sorted["non_missing_count"]>=2]
dat = dat_sorted.drop(columns='non_missing_count')
print(dat)
plt.figure(figsize=(4,20))
# sns.set(font_scale=1.4)
sns.heatmap(dat,yticklabels=1,cmap="binary",vmin=1,vmax=1.01,linewidth=0.5,cbar=False)
plt.savefig("gprofiler.heatmap.png",dpi=300, bbox_inches="tight")



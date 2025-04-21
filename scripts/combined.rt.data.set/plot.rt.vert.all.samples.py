import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle



df = pd.read_csv("rt.vert.all.samples.txt",sep="\t",header=0)

#names=["chrom","start","end","num",
#	"list","acp6","acp7","gm12878","eb3_2","c57/b6"])
df = df[df["num"]>1]
df = df.set_index(["chrom","start","end"])
print(df.loc[:,["acp6","acp7","gm12878","eb3_2","c57/b6"]])

plt.figure()
sns.clustermap(df.loc[:,["acp6","acp7","gm12878","eb3_2","c57/b6"]],
	figsize=(6,12),dendrogram_ratio=0.02,yticklabels=1,cmap="binary",cbar=False,vmin=0,vmax=1)
plt.savefig("shared.verts.png",dpi=300, bbox_inches="tight")

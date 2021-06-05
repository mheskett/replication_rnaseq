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
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]

colors =["#208eb7", "#15dec5", "#214d4e", "#9be0e6", "#30408d", "#ec9fe7", "#e33ab9", "#e3cede", "#830c6f", "#a2e67c", "#519169", "#1be46d", "#65a10e", "#754819", "#bb8553", "#af3007", "#f1d438", "#fb7810", "#fd1e6e", "#f87574", "#432ab7", "#5979fe", "#7d0af6"]
color_dict={chromosomes[i]: colors[i] for i in range(len(chromosomes))}

#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
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

df = pd.read_table("/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/all.samples.overlap.sorted.txt",sep=" ",header=None,index_col=None,
	names= ["chrom","start","stop","1","2","3","4","5","6","7","8"])
df["color"] = df.apply(lambda x: color_dict[x["chrom"]],axis=1)

sns.clustermap(df.loc[:,["1","2","3","4","5","6","7","8"]].transpose(),col_colors=df["color"],col_cluster=True,cmap="binary")

# sns.heatmap(df.loc[:,["1","2","3","4","5","6","7","8"]].transpose())
plt.show()
print(df)
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


all_files = list(glob.glob("/Users/mike/replication_rnaseq/all.final.data/vlinc.calls" + "/*vlinc.discovery.all.bed"))
# all_files = [x  for x in all_files if "rep" not in x]
all_files.sort()
filenames=[os.path.basename(x)[0:15] for x in all_files]
li=[]
number=[]
skew=[]
for i in range(len(all_files)):
	df = pd.read_csv(all_files[i],sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
							"hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
							"l1_fraction":float,"hap1_counts":int,"hap2_counts":int,"reject":str})
	li.append(df)

print(li)

samples = filenames
barWidth = 0.4
num_total = [len(x) for x in li]
number_skewed = [len(x[x["reject"]=="True"]) for x in li]
print(num_total)
print(number_skewed)

r1 = np.arange(len(num_total))
r2 = [x + barWidth for x in r1]

f,ax=plt.subplots(figsize=(12,4))
# Create bars
plt.bar(r1, num_total,color='mediumblue', width=barWidth, edgecolor='white', label='Total vlincs')
plt.bar(r2, number_skewed,color='Black', width=barWidth, edgecolor='white', label='Skewed vlincs')
# Create names on the x-axis
plt.xticks([r + barWidth/2 for r in range(len(num_total))], samples)
plt.xticks(rotation = 315) # Rotates X-Axis Ticks by 45-degrees
plt.legend()

# Show graphic
plt.show()




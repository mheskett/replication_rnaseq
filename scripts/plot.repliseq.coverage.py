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
from sys import argv
import glob
def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)
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

#######################
all_files = glob.glob("/Users/mike/replication_rnaseq" + "/*s100.coverage.bed")
all_files.sort()
filenames=[os.path.basename(x)[0:15] for x in all_files]
li = {}

for i in range(len(all_files)):
	li[filenames[i]] = pd.read_table(all_files[i],
		sep="\t", names=["chrom","start","stop","reads_"+filenames[i]],
		dtype={"chrom":str,"start":int,'stop':int,"reads":int},
		header=None,index_col=None).set_index(['chrom','start','stop'])

df = pd.concat(li.values(),axis=1).reset_index()
df= df[~df.chrom.str.contains("GL")]
df= df[~df.chrom.str.contains("Y")]
df= df[~df.chrom.str.contains("MT")]

df = df.set_index(["chrom","start","stop"])
df_for_coverage = df.copy()
print(df_for_coverage.describe(percentiles=[0.025,0.25,0.5,0.75,0.975]))
for name in df.columns:
	if "bouha" in name:
		plt.hist(df_for_coverage[name])
plt.xlim([0,30000])
plt.show()
plt.close()
for name in df.columns:
	if "gm" in name:
		plt.hist(df_for_coverage[name])
# plt.xlim([0,30000])
plt.show()
plt.close()

df = df.apply(lambda x: x/x.sum() * 10**6, axis=0 )### library size normalization
df_logr = []
sample_names = []
for i in range(0,len(filenames),2):
	sample_names += [filenames[i]+"_log2r"]
	# df_logr += [pd.DataFrame(np.log2( (df["reads_"+filenames[i]]+1) / (df["reads_"+filenames[i+1]]+1) ).rename("log2r").reset_index())] # rename("log2r").reset_index()
	df_logr += [pd.DataFrame(np.log2( (df["reads_"+filenames[i]]+1) / (df["reads_"+filenames[i+1]]+1) ).rename("log2r_"+str(i)))] # rename("log2r").reset_index()

qnorm_on_ratios = quantile_normalize(pd.concat(df_logr,axis=1)).reset_index()

for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))
	for j in range(0,len(filenames),2):
		plt.plot(qnorm_on_ratios[qnorm_on_ratios["chrom"]==chromosomes[i]]["start"],
		 qnorm_on_ratios[qnorm_on_ratios["chrom"]==chromosomes[i]]["log2r_"+str(j)],label=filenames[j])
	plt.legend()
	plt.show()
plt.close()

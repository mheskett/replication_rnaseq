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
import statsmodels.api as sm
from sys import argv

# ## list of imprinted genes
# df_imprinted = pd.read_table("/Users/heskett/replication_rnaseq/scripts/imprinted.genes.fixed.bed",
# 	sep="\t", names=["chrom","start","stop","gene_name"],dtype={"chrom":str},
# 	header=None,index_col=None)
# print(df_imprinted)

## list of chromosome names
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
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

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

arm_dict = get_arms(cytoband)


## example file df = "/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed",sep="\t",


df = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/bouha.expression/bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed",sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
df_auto = df[df["chrom"]!="X"]
df_x = df[df["chrom"]=="X"]
####
# now make a plot with skew, rpkm, strand, size, qval?
## triangle if X linked?
f,ax = plt.subplots(figsize=(8,8))

for table in [df_auto,df_x]:
	color_vector = []
	for index,row in table.iterrows():
		if row["strand"] == "+" and row["reject"]==True:
			color_vector +=[ (1,0,0,1) ]
		if row["strand"] == "+" and row["reject"]==False:
			color_vector +=[ (1,0,0,0.2) ]
		if row["strand"] == "-" and row["reject"]==True:
			color_vector +=[ (0,1,0,1) ]
		if row["strand"] == "-" and row["reject"]==False:
			color_vector +=[ (0,1,0,0.2) ]

	size_vector = [10 if x > 0.05 else 20 + 2*(-np.log10(x)) for x in table["qval"]  ]
	alpha_vector = [0.1 if x=="False" else None for x in table["reject"]]
	if "X" in table["chrom"].values:
		marker = '^'
	else:
		marker = 'o'
	plt.scatter(table["skew"],
				table["stop"] - table["start"],
				lw=0.2,
				edgecolor="black",
				c = color_vector,
				s = size_vector,
				marker = marker
				)
	plt.ylim([0,10**6])
	plt.xlim([-0.5,0.5])
	plt.xticks([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5],fontsize=16)
	plt.yticks([0,250000,500000,750000,1000000],fontsize=16)

plt.savefig(os.path.basename(argv[1].replace(".bed",".jpg")), dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)

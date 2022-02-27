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


# example file ../bouhassira_data/bouha.expression/bouha.2.all.vlincs.and.repliseq.bed

df = pd.read_csv(argv[1],sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
					"hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew","hap1_early","hap2_early","hap1_late","hap2_late"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
					"l1_fraction":float,"hap1_counts":int,"hap2_counts":int,"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})

## annoying data type manipulation
tmp = df.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
tmp = tmp.astype(int)
df.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
###########
df.loc[:,"logr_hap1"] = df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
df.loc[:,"logr_hap2"] = df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
df.loc[:,"logr_diff"] = abs(df["logr_hap1"] - df["logr_hap2"])
##########

df_auto = df[df["chrom"]!="X"]
df_x = df[df["chrom"]=="X"]
####

# f,ax = plt.subplots(figsize=(8,8))
# sns.clustermap(df.loc[:,["skew","logr_diff"]],cmap="RdBu_r")
# plt.show()
# plt.close()
f,ax = plt.subplots(figsize=(8,8))
ax.set_xlim([-.5,.5])
sns.kdeplot(data=df,x="skew",cut=0,linewidth=4)
plt.savefig(os.path.basename(argv[1].replace(".bed",".skew.hist.jpg")), dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()

f,ax = plt.subplots(figsize=(8,8))
# ax.set_xlim([-.5,.5])
sns.kdeplot(data=df,x="logr_diff",cut=0,linewidth=4)
plt.savefig(os.path.basename(argv[1].replace(".bed",".logr.diff.hist.jpg")), dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
# now make a plot with skew, rpkm, strand, size, qval?
## triangle if X linked?
f,ax = plt.subplots(figsize=(8,8))

for table in [df_auto,df_x]:
	color_vector = []
	for index,row in table.iterrows():
		if row["chrom"] == "X" and row["reject"]==True:
			color_vector +=[ (1,0,0,1) ]
		if row["chrom"] == "X" and row["reject"]==False:
			color_vector +=[ (1,0,0,0.1) ]
		if row["chrom"] != "X" and row["reject"]==True:
			color_vector +=[ (0,1,0,1) ]
		if row["chrom"] != "X" and row["reject"]==False:
			color_vector +=[ (0,1,0,0.1) ]
	size_vector = [10+x for x in table["rpkm"]  ]
	alpha_vector = [0.1 if x=="False" else None for x in table["reject"]]
	# if "X" in table["chrom"].values:
	# 	marker = '^'
	# else:
	# 	marker = 'o'
	plt.scatter(table["skew"],
				np.log10(table["stop"] - table["start"]),
				lw=0.2,
				edgecolor="black",
				c = color_vector,
				s = size_vector,
				marker = "o"
				)
	# plt.ylim([0,10**6])
	plt.xlim([-0.5,0.5])
	plt.xticks([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5],fontsize=16)
	# plt.yticks([0,250000,500000,750000,1000000],fontsize=16)

plt.savefig(os.path.basename(argv[1].replace(".bed",".jpg")), dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
##########
f,ax = plt.subplots(figsize=(8,8))


validated_list = ["235_plus_bouha.trim.2Aligned.out.samtool.rmdup"]

for table in [df_x,df_auto]:#,df_x
	color_vector = []
	for index,row in table.iterrows():
		if row["chrom"] == "X" and row["reject"]==True:
			color_vector +=[ (1,0,0,1) ]
		if row["chrom"] == "X" and row["reject"]==False:
			color_vector +=[ (1,0,0,0.1) ]
		if row["chrom"] != "X" and row["reject"]==True:
			color_vector +=[ (0,1,0,1) ]
		if row["chrom"] != "X" and row["reject"]==False:
			color_vector +=[ (0,1,0,0.1) ]

	# size_vector = [10+x for x in table["rpkm"]  ]
	alpha_vector = [0.1 if x=="False" else None for x in table["reject"]]
	# if "X" in table["chrom"].values:
	# 	marker = '^'
	# else:
	# 	marker = 'o'
	plt.scatter(abs(table["skew"]),
				table["logr_diff"],
				lw=0.2,
				edgecolor="black",
				c = color_vector,
				s = 30,
				marker = "o"
				)
	# plt.text(df[df["name"]=="235_plus_bouha.trim.2Aligned.out.samtool.rmdup"]["skew"],
	# 	df[df["name"]=="235_plus_bouha.trim.2Aligned.out.samtool.rmdup"]["logr_diff"],
	# 	s="vlinc187",size=16)
	# plt.ylim([0,2])
	plt.xlim([-0.5,0.5])
	plt.xticks([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5],fontsize=16)
	#plt.yticks(,fontsize=16)
plt.show()
plt.savefig(os.path.basename(argv[1].replace(".bed","repli.jpg")), dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()
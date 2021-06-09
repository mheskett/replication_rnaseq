import os
import re
import csv
import numpy as np
import pandas as pd
import re
import seaborn as sns
import scipy.stats
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
import pickle
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

def add_binom_pval(df):
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
							row["hap1_counts"]+row["hap2_counts"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)
	switcherss = mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")
	df["fdr_pval"] = switcherss[1]
	df["fdr_reject"] = switcherss[0]

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
def helper_func(x):
	if x["total_reads"]==0: # try this for filtering
		return 0
	elif x["hap1_counts"] >= x["hap2_counts"]:
		return x["hap1_counts"]  / x["total_reads"] - 0.5
	else:
		return -x["hap2_counts"]  / x["total_reads"] + 0.5
	return

arm_dict = get_arms(cytoband)

coding_files=["bouha2.protein.coding.all.counts.bed",
"bouha3.protein.coding.all.counts.bed",
"bouha4.protein.coding.all.counts.bed",
"bouha10.protein.coding.all.counts.bed",
"bouha13.protein.coding.all.counts.bed",
"bouha15.protein.coding.all.counts.bed"]
coding_dfs = []
for i in range(len(coding_files)):
	coding_df = pd.read_csv(coding_files[i],sep="\t",
							names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	coding_df["total_reads"] = coding_df["hap1_counts"] + coding_df["hap2_counts"]
	coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
	coding_df["sample"] = coding_files[i][0:7]
	add_binom_pval(coding_df)
	coding_dfs += [coding_df]
df_coding = pd.concat(coding_dfs)



######
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))

df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
	axis=1)
print(df_coding)
#######
df_coding = df_coding[df_coding["total_reads"]>=10]
unique_genes = list(df_coding["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
# df_significant_rows = df_coding[df_coding["binom_pval"]<=0.001]
df_significant_rows = df_coding[df_coding["significant_deviation"]==True]
# df_nonsignificant_rows = df_coding[df_coding["binom_pval"] >=0.001]
df_nonsignificant_rows = df_coding[df_coding["significant_deviation"]==False]
biallelic = df_coding[df_coding["binom_pval"]>=0.001]
### switchers algorithm
for i in range(len(unique_genes)):
	samples = df_significant_rows[df_significant_rows["name"]==unique_genes[i]]
	# samples = samples[(samples["binom_pval_plus"]<=0.05) | (samples["binom_pval_minus"] <=0.05)]
	if len(samples)<=1:
		continue
	# samples = samples.reset_index(drop=True)
	# print(samples)
	hap1_skew,hap2_skew= False,False
	for index,row in samples.iterrows():
		if (row["skew"]>=0.1):
			hap1_skew = True
		if (row["skew"]<=-0.1):
			hap2_skew = True

	if hap1_skew and hap2_skew:
		switchers += [samples]
	elif hap1_skew ^ hap2_skew:
		nonswitchers += [samples]

switchers = pd.concat(switchers)
nonswitchers = pd.concat(nonswitchers)
######
## get the genes list
genes = list(switchers["name"].drop_duplicates())
genes = [x.split(",")[0] for x in genes ]
# print("\n".join(genes))
###

# nonswitching gene list
genes = list(nonswitchers["name"].drop_duplicates())
genes = [x.split(",")[0] for x in genes ]
# print("\n".join(genes))
# print(len(genes))
#plotting shit

color_dict = {"bouha4.":"r","bouha15":"c","bouha10":"y","bouha3.":"g",
"bouha2.":"b","bouha13":"m"}
switchers["color"]= [color_dict[x] for x in switchers["sample"]]
nonswitchers["color"]= [color_dict[x] for x in nonswitchers["sample"]]
biallelic["color"]= [color_dict[x] for x in biallelic["sample"]]

legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha13',markerfacecolor='r', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha15',markerfacecolor='c', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha3',markerfacecolor='g', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha4',markerfacecolor='m', markersize=10)]
##### plot coding switchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = switchers[switchers["chrom"]==chromosomes[i]]
	for index,row in tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
		# rectvar = Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.16,
	 #                 facecolor="gray", edgecolor="gray",hatch="/",fill=False) ## plot vlincs as rectangles
		ax.add_patch(rect)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_ylim([-0.6,0.6])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	plt.savefig("bouha.coding.switchers.region."+str(chromosomes[i])+ ".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
switchers.to_csv("bouha.coding.switchers.bed",sep="\t",index=False,header=False)
###### nonswitchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	for index,row in tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
		ax.add_patch(rect)
	ax.set_ylim([-0.6,0.6])
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	plt.savefig("bouha.coding.nonswitchers.region."+str(chromosomes[i])+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
nonswitchers.to_csv("bouha.coding.nonswitchers.bed",sep="\t",index=False,header=False)
##############
# plot BIALLELCI
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = biallelic[biallelic["chrom"]==chromosomes[i]]
	for index,row in tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
		# rectvar = Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.16,
	 #                 facecolor="gray", edgecolor="gray",hatch="/",fill=False) ## plot vlincs as rectangles
		ax.add_patch(rect)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_ylim([-0.6,0.6])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	plt.savefig("bouha.coding.biallelic.region."+str(chromosomes[i])+ ".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
biallelic.to_csv("bouha.coding.biallelic.bed",sep="\t",index=False,header=False)
##############
##############
## plot sorted list of switchers
switching_loci = switchers.loc[:,["chrom","start","stop"]].apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)

print("variance of the skew of all switchers", switchers["skew"].var())

test = df_coding.pivot(index=["chrom","start","stop"],columns="sample",values="skew").reset_index() # yay this works 
test["position"] = test.apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)
test = test[test["position"].isin(switching_loci)].drop("position",axis=1)#.set_index(["chrom","start","stop"])
# sns.clustermap(test[test["chrom"]=="1"].drop(["chrom","start","stop"],axis=1).fillna(0),cmap="RdBu_r")
# plt.show()
chrom_color_dict = {"1":"red","2":"blue","3":"red","4":"blue","5":"red","6":"blue","7":"red","8":"blue","9":"red","10":"blue","11":"red","12":"blue",
				"13":"red","14":"blue","15":"red","16":"blue","17":"red","18":"blue","19":"red","20":"blue","21":"red","22":"blue","X":"red"}
colors = test.sort_values(["chrom","start"]).apply(lambda x: chrom_color_dict[x["chrom"]],axis=1)
print(test.sort_values(["chrom","start"]))
sns.clustermap(test.sort_values(["chrom","start"]).drop(["chrom","start","stop"],axis=1).transpose().fillna(0),cmap="RdBu_r",
	col_cluster=False,row_cluster=False,figsize=(12,3),vmin=-.25,vmax=.25,col_colors = colors)
plt.show()
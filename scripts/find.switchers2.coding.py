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
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))

coding_files=["bouha2.protein.coding.all.counts.bed",
                "bouha3.protein.coding.all.counts.bed",
                "bouha4.protein.coding.all.counts.bed",
                "bouha10.protein.coding.all.counts.bed",
                "bouha13.protein.coding.all.counts.bed",
                "bouha15.protein.coding.all.counts.bed"]

## GM FILES
# coding_files=["gm12878.4.protein.coding.all.counts.bed",
#                 "gm12878.5.protein.coding.all.counts.bed"]
coding_dfs = []
for i in range(len(coding_files)):
    coding_df = pd.read_csv(coding_files[i],sep="\t",
                            names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    coding_df["total_reads"] = coding_df["hap1_counts"] + coding_df["hap2_counts"]
    coding_df["skew"] = coding_df.apply(helper_func, axis = 1)
    coding_df["sample"] = coding_files[i][0:9]
    add_binom_pval(coding_df)
    coding_dfs += [coding_df]
df_coding = pd.concat(coding_dfs)
df_coding = df_coding[df_coding["total_reads"]>=20]
df_coding = df_coding[df_coding["chrom"]!="X"]
df_coding["reads_per_kb"] = df_coding["total_reads"] / ((df_coding["stop"] - df_coding["start"]) / 1000 )
df_coding  = df_coding[df_coding["reads_per_kb"]>=1]
df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]])\
    .reshape(1,-1))*2.5 else False,
    axis=1)

# df["reads_per_kb"] = df["total_reads"] / ((df["stop"]-df["start"])/1000)
# sns.kdeplot(df["reads_per_kb"],cut=0,clip=(0,10))
# plt.show()
df_coding=df_coding[df_coding["total_reads"]>=20]
print("number of Bouhassira coding genes that pass minimum read filter ")
print(len(df_coding.loc[:,["chrom","start","stop"]].drop_duplicates()))
######
#######
unique_genes = list(df_coding["name"].drop_duplicates())
switchers = [] # list of rows that are switchers
nonswitchers=[]
# df_significant_rows = df[df["binom_pval"]<=0.001]
# df_nonsignificant_rows = df[df["binom_pval"] >=0.001]
df_coding["significant_deviation"] = df_coding.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
	axis=1)
df_significant_rows = df_coding[df_coding["significant_deviation"]==True]
df_nonsignificant_rows = df_coding[df_coding["significant_deviation"]==False]

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

#####


color_dict = {"bouha4.pr":"green",
"bouha15.p":"royalblue",
"bouha10.p":"red",
"bouha3.pr":"yellow",
"bouha2.pr":"cyan",
"bouha13.p":"plum",
"gm12878.4":"red",
"gm12878.5":"blue"}

hatch_ditct = {"bouha4.pr":"//","bouha15.p":"||","bouha10.p":"\\\\","bouha3.pr":"--",
"bouha2.pr":"oo","bouha13.p":".."}
switchers["color"]= [color_dict[x] for x in switchers["sample"]]
nonswitchers["color"]= [color_dict[x] for x in nonswitchers["sample"]]
switchers["hatch"] = [hatch_ditct[x] for x in switchers["sample"]]
df_coding["color"]= [color_dict[x] for x in df_coding["sample"]]


#### make scatter plot dots including faded dots
tmp = df_coding[df_coding["name"].isin(switchers["name"])]
tmp["color"]= [color_dict[x] for x in tmp["sample"]]
tmp["alpha"] = [1 if x==True else 0.1 for x in tmp["significant_deviation"]]
tmp["unique_pos"] = [row["chrom"]+":"+str(row["start"]) for index,row in tmp.iterrows()]
f,ax=plt.subplots(1,1,figsize=(10,2))
print(tmp)
ax.scatter(tmp["unique_pos"],tmp["skew"],c=tmp["color"],s=25,edgecolor="black",lw=0.1,zorder=3,alpha=tmp["alpha"])
for index,row in tmp.drop_duplicates(["unique_pos"]).iterrows():
	ax.axvline(x=row["unique_pos"],linestyle="--",lw=0.4,c="black")
plt.xticks(rotation = 315,fontsize=5)
ax.set_xticklabels(tmp.drop_duplicates(["unique_pos"])["name"])
ax.margins(x=.015,y=0)
ax.set_ylim([-0.52,.52])
ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
ax.set_yticks([-0.5,-.25,0,.25,.5])
plt.savefig("bouha.coding.switchers.all.test.alpha.png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()

print("number of switching loci ")
print(len(tmp["unique_pos"].drop_duplicates()))
#exit()######## EXIT EXIT EXIT EXIT EXIT EXIT EXIT EXIT
## biallelic
result=[]
for index,row in df_coding.iterrows():
	if sum(list(df_coding[df_coding["name"]==row["name"]]["significant_deviation"]))==0:
		result+=[False]
	else:
		result+=[True]
df_coding["any_significant_deviation"] = result

tmp = df_coding[df_coding["any_significant_deviation"]==False]
tmp["color"]= [color_dict[x] for x in tmp["sample"]]
tmp["unique_pos"] = [row["chrom"]+":"+str(row["start"]) for index,row in tmp.iterrows()]
f,ax=plt.subplots(1,1,figsize=(10,2))
ax.scatter(tmp["unique_pos"],tmp["skew"],c=tmp["color"],s=25,edgecolor="black",lw=0.1,zorder=3)
for index,row in tmp.drop_duplicates(["unique_pos"]).iterrows():
	ax.axvline(x=row["unique_pos"],linestyle="--",lw=0.4,c="black")
plt.xticks(rotation = 315,fontsize=5)
ax.margins(x=.015,y=0)
ax.set_ylim([-0.52,.52])
ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
ax.set_yticks([-0.5,-.25,0,.25,.5])
plt.savefig("bouha.coding.biallelic.all.test.alpha.png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
plt.close()


# legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha13',markerfacecolor='r', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha15',markerfacecolor='c', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha3',markerfacecolor='g', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10),
# Line2D([0], [0], marker='o', color='w', label='bouha4',markerfacecolor='m', markersize=10)]

##### plot switchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = switchers[switchers["chrom"]==chromosomes[i]]
	for index,row in tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.025), width=row["stop"]-row["start"], height=0.05,
	                 facecolor=row["color"],fill=True) ## plot vlincs as rectangles
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

test = list(switchers["tmp"])
print(test)
exit()
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
##############
##############
## plot sorted list of switchers
switching_loci = switchers.loc[:,["chrom","start","stop"]].apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)

print("variance of the skew of all switchers", switchers["skew"].var())

test = df.pivot(index=["chrom","start","stop"],columns="sample",values="skew").reset_index() # yay this works 
test["position"] = test.apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)
test = test[test["position"].isin(switching_loci)].drop("position",axis=1)#.set_index(["chrom","start","stop"])
# sns.clustermap(test[test["chrom"]=="1"].drop(["chrom","start","stop"],axis=1).fillna(0),cmap="RdBu_r")
# plt.show()
# chrom_color_dict = {"1":"red","2":"blue","3":"red","4":"blue","5":"red","6":"blue","7":"red","8":"blue","9":"red","10":"blue","11":"red","12":"blue",
# 				"13":"red","14":"blue","15":"red","16":"blue","17":"red","18":"blue","19":"red","20":"blue","21":"red","22":"blue","X":"red"}
# colors = test.sort_values(["chrom","start"]).apply(lambda x: chrom_color_dict[x["chrom"]],axis=1)
# print(test.sort_values(["chrom","start"]))
# sns.clustermap(test.sort_values(["chrom","start"]).drop(["chrom","start","stop"],axis=1).transpose().fillna(0),cmap="RdBu_r",
# 	col_cluster=False,row_cluster=False,figsize=(12,3),vmin=-.25,vmax=.25,col_colors = colors)
# plt.show()
##############
regions = [
	["1",186900000,187800000],
	["1",238000000,239000000],
	["2",80400000,81100000], 
	["2",227300000,228500000],
	["2",154400000,154900000],
	["5",72000000,75000000],
	["5",96000000,98000000],
	['5',8300000,9700000],
	["5",84400000,85000000],
	["5",157500000,158700000],
	["6",21000000,23000000],
	["6",75000000,76000000],
	["6",26000000,36000000], ## this is the IG locus
	["8",62000000,64000000],
	["8",119000000,120500000],
	["14",106052774,107288051], # IGH
	["2",89156874,90274235], #IGK
	["22",22380474,23265085], #IGL
	["21",17400000,18200000],
	["21",29600000,31000000],
	]

# ### for zooming in on regions
# for i in range(len(regions)):
# 	chrom=regions[i][0]
# 	start=regions[i][1]
# 	stop=regions[i][2]
# ### switchers RNA and vlincs
# 	f,ax = plt.subplots(1,figsize=(1.5,2),sharex=False)
# 	plt.suptitle(chrom)
# 	plt.xticks(rotation=70,fontsize=6)
# 	# ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
# 	# ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
# 	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
# 	ax.set_xlim([start, stop])
# 	ax.set_ylim([-.52,.52])
# 	ax.tick_params(axis="x", labelsize=6,labelrotation=335) 
# 	ax.set_xticks(list(np.linspace(start,stop, 4)))
# 	print("plotting vlincs ")
# 	for index,row in switchers[(switchers["chrom"]==chrom) & (switchers["start"]>=start-2000000) & (switchers["stop"]<=stop+2000000)
# 					& (switchers["fdr_pval"]<=0.01)].iterrows():
# 		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
#                      facecolor=row["color"], edgecolor=row["color"],alpha=0.8,fill=False,hatch="/",lw=1)
# 		ax.add_patch(rect)
# 	plt.show()
#### for zooming in on all regions in the switchers df
unique_locations = switchers.loc[:,["chrom","start","stop"]].drop_duplicates()
print("number of non-coding switching loci",len(unique_locations))
print("base pairs of noncoding switching aDE", (unique_locations["stop"] - unique_locations["start"]).sum())

## without showing biallelics
for index,row in unique_locations.iterrows():
	chrom=row["chrom"]
	start=row["start"]
	stop=row["stop"]
### switchers RNA and vlincs
	f,ax = plt.subplots(1,figsize=(1.5,2),sharex=False)
	plt.suptitle(chrom)
	plt.xticks(np.linspace(start-100000,stop+100000,4),rotation=70,fontsize=8)
	# ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
	# ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([start-100000, stop+100000])
	ax.set_ylim([-.52,.52])
	ax.tick_params(axis="x", labelsize=6,labelrotation=335) 
	print("plotting coding ")
	### july addition: make the rectangles thinner and add jitter?
	for index,row in switchers[(switchers["chrom"]==chrom) & (switchers["start"]>=start-2000000) & (switchers["stop"]<=stop+2000000)
					& (switchers["fdr_pval"]<=0.01)].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.025), width=row["stop"]-row["start"], height=0.05,
                     facecolor=row["color"], edgecolor="black",alpha=1,fill=True,lw=0.5) #
		ax.add_patch(rect)
	plt.savefig("bouha.coding.switchers.region."+str(chrom)+"."+str(start)+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
###
###
### with showing biallelics
for index,row in unique_locations.iterrows():
	chrom=row["chrom"]
	start=row["start"]
	stop=row["stop"]
### switchers RNA and vlincs
	f,ax = plt.subplots(1,figsize=(1.5,2),sharex=False)
	plt.suptitle(chrom)
	plt.xticks(np.linspace(start-100000,stop+100000,4),rotation=70,fontsize=8)
	# ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
	# ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([start-100000, stop+100000])
	ax.set_ylim([-.52,.52])
	ax.tick_params(axis="x", labelsize=6,labelrotation=335) 
	for index,row in df[(df["chrom"]==chrom) & (df["start"]>=start-2000000) & (df["stop"]<=stop+2000000)].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.025), width=row["stop"]-row["start"], height=0.05,
                     facecolor=row["color"], edgecolor="black",alpha=1 if row["significant_deviation"]==True else 0.1,fill=True,lw=0.5) #
		ax.add_patch(rect)
	plt.savefig("bouha.coding.switchers.region.test.alpha."+str(chrom)+"."+str(start)+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()

exit()

###############

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed"]
filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
	df_repli = pd.read_csv(all_files_repli[i],sep="\t",
						names= ["chrom","start","stop","hap1_early","hap2_early","hap1_late","hap2_late"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str,
						"hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})

	tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
	tmp = tmp.astype(int)
	df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
	df_repli = df_repli.set_index(["chrom","start","stop"])
	df_repli = df_repli[df_repli.sum(axis="columns")!=0]
	df_repli = df_repli.reset_index()
	df_repli["sample"] = filenames_repli[i]
	repli_li.append(df_repli)
repli_df = pd.concat(repli_li)
repli_df.loc[:,"logr_hap1"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
repli_df.loc[:,"logr_hap2"] = repli_df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
repli_df.loc[:,"logr_diff"] = abs(repli_df["logr_hap1"] - repli_df["logr_hap2"]) ## 
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
sig_repli = np.percentile(a = repli_df[repli_df["chrom"]!="X"]["logr_diff"], q = 95)
repli_df["sig_repli"]=["True" if x > sig_repli else "False" for x in repli_df["logr_diff"]]
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector

# coding switchers with asynchrony....
##### plot coding switchers
for i in range(len(chromosomes)):
	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = switchers[switchers["chrom"]==chromosomes[i]]
	for index,row in tmp.iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
	                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
		# rectvar = Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.16,
	 #                 facecolor="gray", edgecolor="gray",hatch="/",fill=False) ## plot vlincs as rectangles
		ax[0].add_patch(rect)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_ylim([-0.6,0.6])
	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j+1].set_ylim([0,1.6])
		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]

		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df[tmp_df["arm"]=="p"])
		smoothed_logrdiff_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
			return_sorted=False, frac = frac1 )

		smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
			return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_logrdiff_q,c="black",zorder=1,lw=0.8)
		ax[j+1].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j+1].add_patch(rect)
	plt.savefig("bouha.coding.switchers.repli."+str(chromosomes[i])+ ".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.show()
### coding non switchers with asynchrony
for i in range(len(chromosomes)):
	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
	plt.suptitle(chromosomes[i])
	tmp = nonswitchers[nonswitchers["chrom"]==chromosomes[i]]
	ax[0].scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j+1].set_ylim([0,1.6])

		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]

		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df[tmp_df["arm"]=="p"])
		smoothed_logrdiff_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
			return_sorted=False, frac = frac1 )

		smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
			return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_logrdiff_q,c="black",zorder=1,lw=0.8)
		ax[j+1].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j+1].add_patch(rect)
	plt.show()

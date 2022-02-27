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
import glob
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as mt
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
	results = mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")
	df["fdr_pval"] = results[1]
	df["fdr_reject"] = results[0]

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

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.500kb.bed"]
filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
	df_repli = pd.read_csv(all_files_repli[i],sep="\t",
						names= ["chrom","start","stop","hap1_early","hap2_early","hap1_late","hap2_late"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str,
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
repli_df.loc[:,"logr_diff_raw"] = repli_df["logr_hap1"] - repli_df["logr_hap2"] # positive if hap1 early, negative if hap2 early

repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector
#####
logr_diff_bouha10= repli_df[repli_df["sample"]=="bouha.10.repli."].loc[:,["chrom","start","stop","logr_diff_raw"]].set_index(["chrom","start","stop"])
logr_diff_bouha2= repli_df[repli_df["sample"]=="bouha.2.repli.5"].loc[:,["chrom","start","stop","logr_diff_raw"]].set_index(["chrom","start","stop"])
logr_diff_bouha10 =  logr_diff_bouha10.rename(columns={"logr_diff_raw":"logr_diff_raw_bouha10"})#.reset_index()
logr_diff_bouha2= logr_diff_bouha2.rename(columns={"logr_diff_raw":"logr_diff_raw_bouha2"})#.reset_index()
#############
colors =["Red", "Blue", "#214d4e", "#9be0e6", "#30408d", "#ec9fe7", "#e33ab9", "#e3cede", "#830c6f", "#a2e67c", "#519169", "#1be46d", "#65a10e", "#754819", "#bb8553", "#af3007", "#f1d438", "#fb7810", "#fd1e6e", "#f87574", "#432ab7", "#5979fe", "#7d0af6"]
###########
############
#####
######
######
## make example plots 
test = repli_df[(repli_df["sample"]=="bouha.10.repli.") & (repli_df["chrom"]=="1")]
if len(test[test["arm"]=="p"]) <= 10:
	frac1 = 1
else:
	frac1 = 6 / len(test[test["arm"]=="p"])
hap1_p = sm.nonparametric.lowess(endog=test[test["arm"]=="p"]["logr_hap1"], 
	exog=test[test["arm"]=="p"]["start"], 
	return_sorted=False, frac = frac1 )
hap1_q = sm.nonparametric.lowess(endog=test[test["arm"]=="q"]["logr_hap1"], 
	exog=test[test["arm"]=="q"]["start"], 
	return_sorted=False, frac = 6/len(test[test["arm"]=="q"]["start"].index) )
hap2_p = sm.nonparametric.lowess(endog=test[test["arm"]=="p"]["logr_hap2"], 
	exog=test[test["arm"]=="p"]["start"], 
	return_sorted=False, frac = frac1 )
hap2_q = sm.nonparametric.lowess(endog=test[test["arm"]=="q"]["logr_hap2"], 
	exog=test[test["arm"]=="q"]["start"], 
	return_sorted=False, frac = 6/len(test[test["arm"]=="q"]["start"].index) )

plt.subplots(figsize=(6,2))
plt.xlim([0, 240000000])
plt.ylim([-3.5,3.5])
plt.xticks(np.linspace(0,240000000,20))
plt.plot(test[test["arm"]=="p"]["start"],hap1_p,c="Red")
plt.plot(test[test["arm"]=="q"]["start"],hap1_q,c="Red")
plt.plot(test[test["arm"]=="p"]["start"],hap2_p,c="Red")
plt.plot(test[test["arm"]=="q"]["start"],hap2_q,c="Red")
test = repli_df[(repli_df["sample"]=="bouha.2.repli.5") & (repli_df["chrom"]=="1")]
print(test)
if len(test[test["arm"]=="p"]) <= 10:
	frac1 = 1
else:
	frac1 = 6 / len(test[test["arm"]=="p"])
hap1_p = sm.nonparametric.lowess(endog=test[test["arm"]=="p"]["logr_hap1"], 
	exog=test[test["arm"]=="p"]["start"], 
	return_sorted=False, frac = frac1 )
hap1_q = sm.nonparametric.lowess(endog=test[test["arm"]=="q"]["logr_hap1"], 
	exog=test[test["arm"]=="q"]["start"], 
	return_sorted=False, frac = 6/len(test[test["arm"]=="q"]["start"].index) )
hap2_p = sm.nonparametric.lowess(endog=test[test["arm"]=="p"]["logr_hap2"], 
	exog=test[test["arm"]=="p"]["start"], 
	return_sorted=False, frac = frac1 )
hap2_q = sm.nonparametric.lowess(endog=test[test["arm"]=="q"]["logr_hap2"], 
	exog=test[test["arm"]=="q"]["start"], 
	return_sorted=False, frac = 6/len(test[test["arm"]=="q"]["start"].index) )
plt.xlim([0, 240000000])
plt.ylim([-3.5,3.5])
plt.xticks(np.linspace(0,240000000,20))
plt.plot(test[test["arm"]=="p"]["start"],hap1_p,c="Blue")
plt.plot(test[test["arm"]=="q"]["start"],hap1_q,c="Blue")
plt.plot(test[test["arm"]=="p"]["start"],hap2_p,c="Blue")
plt.plot(test[test["arm"]=="q"]["start"],hap2_q,c="Blue")



plt.show()
plt.close()

plt.xlim(0,240000000)
plt.plot(test[test["arm"]=="p"]["start"],abs(hap1_p - hap2_p))
plt.show()
#############
############
############
###########
clustermap_df = pd.concat([logr_diff_bouha10,logr_diff_bouha2],axis=1).reset_index().dropna(how="any")
color_vector=[colors[0] if x=="X" else colors[1] for x in clustermap_df["chrom"]]
print(clustermap_df[clustermap_df.isna().any(axis=1)])
clustermap_df["logr_diff_diff"] = abs(clustermap_df["logr_diff_raw_bouha2"] - clustermap_df["logr_diff_raw_bouha10"])
logr_diff_diff_min = np.percentile(a = clustermap_df[clustermap_df["chrom"]!="X"]["logr_diff_diff"], q = 90)
clustermap_df["color"] = color_vector
print(clustermap_df)
f,ax=plt.subplots(1,1,figsize=(3,3))
print(np.log2(abs(clustermap_df[clustermap_df["chrom"]!="X"]["logr_diff_raw_bouha2"])))
# plt.xlim([0,2])
# sns.kdeplot(abs(clustermap_df[clustermap_df["chrom"]!="X"]["logr_diff_raw_bouha2"]),cut=0,lw=3)
# sns.kdeplot(abs(clustermap_df[clustermap_df["chrom"]!="X"]["logr_diff_raw_bouha10"]),cut=0,lw=3)
# plt.hist(np.log2(abs(clustermap_df[(clustermap_df["chrom"]!="X")]["logr_diff_raw_bouha2"]+1)),lw=3)
# sns.kdeplot(np.log2(abs(clustermap_df[clustermap_df["chrom"]!="X"]["logr_diff_raw_bouha10"])),lw=3)
# plt.show()
# plt.savefig("bouha.RT.dist.png",
# 		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.close()

#########
print("distribution of AS_RT switching among the most AS_RT regions of the genome")
sig_repli_bouha2 = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]=="bouha.2.repli.5")]["logr_diff"], q = 90)
sig_repli_bouha10 = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]=="bouha.10.repli.")]["logr_diff"], q = 90)
f,ax=plt.subplots(figsize=(3,3))
plt.xlim([0,4])
sns.kdeplot(clustermap_df[(clustermap_df["chrom"]!="X") & (abs(clustermap_df["logr_diff_raw_bouha2"])>=sig_repli_bouha2) & 
			(abs(clustermap_df["logr_diff_raw_bouha10"])>=sig_repli_bouha10)]["logr_diff_diff"],lw=3)
plt.show()
sns.kdeplot(np.log2(clustermap_df[(clustermap_df["chrom"]!="X") & (abs(clustermap_df["logr_diff_raw_bouha2"])>=sig_repli_bouha2) & 
			(abs(clustermap_df["logr_diff_raw_bouha10"])>=sig_repli_bouha10)]["logr_diff_diff"]),lw=3)########
plt.show()
########
########
#######
print(repli_df)
# coding switchers with asynchrony....
for i in range(len(chromosomes)):
	f,ax = plt.subplots(4,1,figsize=(10,8),sharex=True)
	# plt.suptitle(chromosomes[i])
	# tmp = switchers[switchers["chrom"]==chromosomes[i]]
	# ax[0].scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=30)
	# ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	# ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	# ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 16))
	for j in range(len(filenames_repli)):
		ax[j].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j].set_ylim([0,1.6])
		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]
		print("repli 90%",filenames_repli[j],sig_repli)
		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df[tmp_df["arm"]=="p"])
		smoothed_logrdiff_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
			return_sorted=False, frac = frac1 )

		smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
			return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
			ax[j].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_logrdiff_q,c="black",zorder=1,lw=0.8)
		ax[j].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j].add_patch(rect)
	#plt.show()
	plt.close()

########## specific regions

	regions=[
	        ["1",186000000,188000000],
	        ["8",2400000,2900000],
	        ["1",13500000,16000000],
	        ["1",56000000,57500000],
	        ["10",18500000,21000000],
	        ["10",75000000,81000000],
	        ["10",108000000,112000000],
	        ["1",187047872,187088572], ## regions that matt validated with fish
			["1",187194601,187239328],
 			["1",187416000,187459201],
 			["3",174933729,174972546],
 			["6",40074887,40115741],
 			["6",73394774,73431597],
 			["6",77914535,77951564],
 			["6",120454710,120497591],
 			["6",130585933,130623577],
 			["6",141108314,141151646],
 			["6",141184172,141225657],
 			["8",2586332,2625241],
 			["8",22568321,22607943],
 			["9",22765327,22807785],
 			["9",23643340,23683461],
 			["9",24038517,24077283],
 			["9",29802888,29843986],
 			["15",47102705,47144857],
 			["15",54386886,54427639],
 			["15",54629690,54670049],
 			["15",92002424,92044668],
 			["15",96464802,96501683],
 			["15",97096204,97134188],

	         ]
####

for i in range(len(regions)):

	f,ax = plt.subplots(4,1,figsize=(2,8),sharex=False)
	chrom = regions[i][0]
	start = regions[i][1]
	stop = regions[i][2]
	plt.suptitle(chrom)
#### repli
	for j in range(len(filenames_repli)):
		ax[j].set_xlim([start - 500000, stop + 500000])
		ax[j].set_ylim([0,1.5])

		tmp_df = repli_df[(repli_df["chrom"]==chrom) & (repli_df["sample"]==filenames_repli[j]) & (repli_df["start"]>=start-5000000 ) & (repli_df["stop"]<=stop+5000000  )] # chromosomes and sample specific now.
		sig_repli = np.percentile(a = repli_df[(repli_df["chrom"]!="X") & (repli_df["sample"]==filenames_repli[j])]["logr_diff"], q = 95)
		tmp_df["sig_repli"]=["True" if x > sig_repli else "False" for x in tmp_df["logr_diff"]]
		print(tmp_df)
		if len(tmp_df) <= 10:
			frac1 = 1
		else:
			frac1 = 4 / len(tmp_df)
		smoothed_logrdiff = sm.nonparametric.lowess(endog=tmp_df["logr_diff"], exog=tmp_df["start"], 
			return_sorted=False, frac = frac1 )

		# smoothed_logrdiff_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_diff"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
		# 	return_sorted=False, frac = 4/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

		# if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
		# 	ax[j].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_logrdiff_p,c="black",zorder=1,lw=0.8)
		ax[j].plot(tmp_df["start"],smoothed_logrdiff,c="black",zorder=1,lw=4)
		ax[j].set_xticks(np.linspace(start - 500000, stop + 500000, 4))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j].add_patch(rect)
	plt.show()

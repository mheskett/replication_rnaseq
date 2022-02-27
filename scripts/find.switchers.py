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
import glob
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
from matplotlib.lines import Line2D

### might not need stranded analysis at all for switchers analysis. 
def helper_func(x):
	if x["total_reads"]==0: # try this for filtering
		return 0
	elif x["hap1_counts"] >= x["hap2_counts"]:
		return x["hap1_counts"]  / x["total_reads"] - 0.5
	else:
		return -x["hap2_counts"]  / x["total_reads"] + 0.5
	return
def helper_func_plus(x):
	if x["total_plus"]<10:
		return 0
	elif x["hap1_counts_plus"] >= x["hap2_counts_plus"]:
		return x["hap1_counts_plus"]  / x["total_plus"] - 0.5
	else:
		return -x["hap2_counts_plus"]  / x["total_plus"] + 0.5
	return

def helper_func_minus(x):
	if x["total_minus"]<10:
		return 0
	elif x["hap1_counts_minus"] >= x["hap2_counts_minus"]:
		return x["hap1_counts_minus"]  / x["total_minus"] - 0.5
	else:
		return -x["hap2_counts_minus"]  / x["total_minus"] + 0.5
	return

def add_binom_pval(df):
	df["binom_pval_plus"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts_plus"],
							row["hap1_counts_plus"]+row["hap2_counts_plus"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)
	df["binom_pval_minus"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts_minus"], ## making minus strand pval negative for ease of analysis.
							row["hap1_counts_minus"]+row["hap2_counts_minus"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)
	return

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
#########################
arm_dict = get_arms(cytoband)
#### coding genes just for plotting
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
	coding_df["sample"] = coding_files[i][0:15]
	add_binom_pval(coding_df)
	dfs += [coding_df]
df_coding = pd.concat(dfs)

### repli
# all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.4.repli.250kb.bed",
# "/Users/mike/replication_rnaseq/all.final.data/gm12878.5.repli.250kb.bed"]

all_files_repli = ["/Users/mike/replication_rnaseq/all.final.data/bouha.10.repli.500kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/bouha.2.repli.500kb.bed"]
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
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

#### loading repli

# all_files = ["/Users/mike/replication_rnaseq/all.final.data/gm12878.4.rna.50kb.bed",
# 			"/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rna.50kb.bed"]
all_files = ["/Users/mike/replication_rnaseq/all.final.data/bouha10.rna.50kb.bed",
			"/Users/mike/replication_rnaseq/all.final.data/bouha2.rna.50kb.bed",
		"/Users/mike/replication_rnaseq/all.final.data/bouha3.rna.50kb.bed",
		"/Users/mike/replication_rnaseq/all.final.data/bouha4.rna.50kb.bed",
		"/Users/mike/replication_rnaseq/all.final.data/bouha13.rna.50kb.bed",
		"/Users/mike/replication_rnaseq/all.final.data/bouha15.rna.50kb.bed"]

all_vlincs = glob.glob("/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/bouha*all.bed")
# all_vlincs = glob.glob("/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.*all.bed")

# ["bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.skewed.bed",
# "bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.10Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.skewed.bed"]
filenames_vlincs=[os.path.basename(x)[0:15] for x in all_vlincs]
vlinc_li = []
for i in range(len(all_vlincs)):
	vlinc_df = pd.read_csv(all_vlincs[i],sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
						"hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
						"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
	vlinc_df["sample"] = filenames_vlincs[i]
	vlinc_li.append(vlinc_df)
vlincs_df = pd.concat(vlinc_li)
# df = pd.read_csv("/Users/mike/replication_rnaseq/all.final.data/gm12878.4.rna.50kb.bed",sep="\t",
# 						names= ["chrom","start","stop","hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"],
# 						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str})
filenames=[os.path.basename(x)[0:15] for x in all_files]
li = []
wide_li = []
for i in range(len(all_files)):
	df = pd.read_csv(all_files[i],sep="\t",
							names= ["chrom","start","stop","hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str})
	df["sample"] = filenames[i]
	tmp = df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]].replace(".",0)
	tmp = tmp.astype(int)
	df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]] = tmp
	df = df.set_index(["chrom","start","stop"])
	# df = df[df.sum(axis="columns")!=0] 
	df = df[(df.loc[:,["hap1_counts_plus","hap2_counts_plus"]].sum(axis=1)>=10) | (df.loc[:,["hap1_counts_minus","hap2_counts_minus"]].sum(axis=1)>=10)]
	df = df.reset_index()
	df["total_plus"] = df["hap1_counts_plus"] + df["hap2_counts_plus"]
	df["total_minus"] = df["hap1_counts_minus"] + df["hap2_counts_minus"]
	df["skew_plus"] = df.apply(helper_func_plus, axis = 1)
	df["skew_minus"] = df.apply(helper_func_minus, axis = 1)
	li.append(df)
	wide_li.append(df.set_index(["chrom","start","stop"]).loc[:,["skew_plus","skew_minus"]].add_suffix("_"+filenames[i]))
df =pd.concat(li)
df_wide = pd.concat(wide_li,axis=1)
df_wide = df_wide.dropna(how="any",axis=0).reset_index()
print("adding pvals") # I guess should have had a separate script to do this once......
add_binom_pval(df)
print("done adding pvals")
unique_rows = df.loc[:,["chrom","start","stop"]].drop_duplicates()
switchers = [] # list of rows that are switchers
nonswitchers=[]
df_significant_rows = df[(df["binom_pval_plus"]<=0.001) | (df["binom_pval_minus"] <=0.001)]
df_nonsignificant_rows = df[(df["binom_pval_plus"]>=0.001) | (df["binom_pval_minus"] >=0.001)]
significant_unique_rows = df_significant_rows.loc[:,["chrom","start","stop"]].drop_duplicates()


### switchers algorithm
for index,row in unique_rows.iterrows():
	samples = df_significant_rows[(df_significant_rows["chrom"]==row["chrom"]) & (df_significant_rows["start"]==row["start"]) & (df_significant_rows["stop"]==row["stop"]) ]
	# samples = samples[(samples["binom_pval_plus"]<=0.05) | (samples["binom_pval_minus"] <=0.05)]
	if len(samples)<=1:
		continue
	samples = samples.reset_index(drop=True)
	# print(samples)
	hap1_plus_skew,hap2_plus_skew,hap1_minus_skew,hap2_minus_skew = False,False,False,False
	for index,row in samples.iterrows():
		if (row["skew_plus"]>=0.1 and row["binom_pval_plus"]<=0.001):
			hap1_plus_skew = True
		if (row["skew_plus"]<=-0.1 and row["binom_pval_plus"]<=0.001):
			hap2_plus_skew = True
		if (row["skew_minus"]>=0.1 and row["binom_pval_minus"]<=0.001):
			hap1_minus_skew = True
		if (row["skew_minus"]<=-0.1 and row["binom_pval_minus"]<=0.001):
			hap2_minus_skew = True
	if (hap1_plus_skew == True or hap1_minus_skew == True) and (hap2_plus_skew ==True or hap2_minus_skew == True):
		switchers += [samples]
	else:
		nonswitchers +=[samples]
result = pd.concat(switchers)
nonswitch_result = pd.concat(nonswitchers)
print(result)
##############
#### heirarchical plot of all samples by skew plus and skew minus?
sig_index = result.loc[:,["chrom","start","stop"]].apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)
print("sig index length",len(sig_index))
df_wide["position"] = df_wide.apply(lambda x: x["chrom"]+":"+str(x["start"])+"-"+str(x["stop"]),axis=1)
df_sig = df_wide[df_wide["position"].isin(sig_index)].drop("position",axis=1).set_index(["chrom","start","stop"])

dfplus = df_sig.filter(like='skew_plus',axis=1)
dfplus.columns = np.arange(0,len(dfplus.columns))
dfplus = dfplus[abs(dfplus.sum(axis=1))>0]
dfmin = df_sig.filter(like="skew_minus",axis=1)
dfmin.columns = np.arange(0,len(dfmin.columns))
dfmin = dfmin[abs(dfmin.sum(axis=1))>0]
test =pd.concat([dfplus ,dfmin])
test=test[test.var(axis=1)>=0.077612] # get top 100 ish?
#f,x=plt.subplots(figsize=(3,10))
############ seaborn heatmap
# sns.set(font_scale=.5)
# sns.clustermap(test,cmap="RdBu_r",vmin=-0.3,vmax=0.3,figsize=(5,12))
# # find rows with the greatest variance?
# #print(test.var(axis=1).sort_values(ascending=False)[0:100])
# plt.savefig("switchers.heatmap.png",
# 		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
# plt.close()
#################

##########
###########

color_dict = {"gm12878.4.rna.5":"plum","gm12878.5.rna.5":"olivedrab","bouha13.rna.50k":"r","bouha15.rna.50k":"c","bouha10.rna.50k":"y","bouha3.rna.50kb":"g",
"bouha2.rna.50kb":"b","bouha4.rna.50kb":"m"}
color_dict_vlincs = {"bouha.trim.13Al":"r","bouha.trim.15Al":"c","bouha.trim.10Al":"y","bouha.trim.3Ali":"g",
"bouha.trim.2Ali":"b","bouha.trim.4Ali":"m","gm12878.4x.hg19":"plum","gm12878.5x.hg19":"olivedrab",
"gm12878.4.rep1.":"plum","gm12878.5.rep1.":"olivedrab",
}
result["color"]= [color_dict[x] for x in result["sample"]]
vlincs_df["color"] = [color_dict_vlincs[x] for x in vlincs_df["sample"]]
sig_repli = np.percentile(a = repli_df[repli_df["chrom"]!="X"]["logr_diff"], q = 95)
repli_df["sig_repli"]=["True" if x > sig_repli else "False" for x in repli_df["logr_diff"]]
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector
legend = [Line2D([0], [0], marker='o', color='w', label='gm12878.4',markerfacecolor='plum', markersize=10),
Line2D([0], [0], marker='o', color='w', label='gm12878.5',markerfacecolor='olivedrab', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha13',markerfacecolor='r', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha15',markerfacecolor='c', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha10',markerfacecolor='y', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha3',markerfacecolor='g', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha2',markerfacecolor='b', markersize=10),
Line2D([0], [0], marker='o', color='w', label='bouha4',markerfacecolor='m', markersize=10)]
#####################################################################
## Plot all biallelic vlincs only
for i in range(len(chromosomes)):
	f,ax = plt.subplots(1,1,figsize=(10,2),sharex=False)
	# tmp_plus = result[(result["chrom"]==chromosomes[i]) & (result["skew_plus"]!=0)]
	# tmp_minus = result[(result["chrom"]==chromosomes[i]) & (result["skew_minus"]!=0)]
	plt.suptitle(chromosomes[i])
	ax.set_ylim([-.52,.52])

	# ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
	# ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
	ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax.set_xlim([0, chromosome_length[chromosomes[i]]])
	ax.set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
	for index,row in vlincs_df[(vlincs_df["chrom"]==chromosomes[i]) & (vlincs_df["qval"]>=0.01) & (abs(vlincs_df["skew"])<=0.1)].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"], alpha=0.8,fill=False,hatch="/",edgecolor=row["color"])
		ax.add_patch(rect)
	# plt.legend(handles=legend,loc=(1.03,0))
	plt.savefig("bouha.biallelic.vlincs."+str(chromosomes[i])+ ".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
########################
######################
# Make some zooms

# regions = [ ## microsatellites in rare disease
# ["12",	6033625,	8053815],
# ["12",	110890017,	113037480],
# ["13",	69681344,	71713885],
# ["14",	91524895,	93572965],
# ["16",	65503746,	67584315],
# ["16",	86635440,	88731761],
# ["18",	51889561,	54303188],
# ["19",	12317255,	14617274],
# ["19",	45272975,	47285815],
# ["20",	1633177,	3639039],
# ["22",	45067677,	47241187],
# ["3",62850232,	64989136],
# ["3",127886657,	129902810],
# ["4",2064972,	4076241],
# ["4",2076407,	4245687],
# ["5",144969066,	147461083],
# ["6",15299342,	17761721],
# ["6",169863420,	171881958],
# ["9",26546542,	28573864],
# ["9",70650478,	72715094],
# ["X",65763873,	67950461],
# ["X",145990948,	148003676],
# ["X",145993468,	149082193]
# ]
regions = [
	["1",97500000,98200000],
	["1",186900000,187800000],
	["1",227000000,227800000],
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
	["22",33000000,35000000]
	]

### for zooming in on regions
for i in range(len(regions)):
	chrom=regions[i][0]
	start=regions[i][1]
	stop=regions[i][2]
### switchers RNA and vlincs
	f,ax = plt.subplots(3,1,figsize=(1.5,6),sharex=False)
	tmp_plus = result[(result["chrom"]==chrom) & (result["skew_plus"]!=0) & (result["start"]>=start) & (result["stop"]<=stop)]
	tmp_minus = result[(result["chrom"]==chrom) & (result["skew_minus"]!=0) & (result["start"]>=start) & (result["stop"]<=stop)]
	plt.suptitle(chrom)
	plt.xticks(rotation=70,fontsize=6)
	# ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
	# ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([start, stop])
	ax[0].set_ylim([-.52,.52])
	ax[0].tick_params(axis="x", labelsize=6,labelrotation=335) 
	ax[0].set_xticks(list(np.linspace(start,stop, 4)))
	print("plotting vlincs ")
	for index,row in vlincs_df[(vlincs_df["chrom"]==chrom) & (vlincs_df["start"]>=start-2000000) & (vlincs_df["stop"]<=stop+2000000)
					& (vlincs_df["qval"]<=0.01)].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"], edgecolor=row["color"],alpha=0.8,fill=False,hatch="/")
		ax[0].add_patch(rect)

#### repli
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([start,stop])
		ax[j+1].set_ylim([-4,4])
		ax[j+1].axhline(y=0,linestyle="--",c="black",lw=0.4)

		tmp_df = repli_df[(repli_df["chrom"]==chrom) & (repli_df["sample"]==filenames_repli[j]) & (repli_df["start"]>=start-500000) & (repli_df["stop"]<=stop+500000)] # chromosomes and sample specific now.
		
		if len(tmp_df) <= 10:
			frac1 = 1
		else:
			frac1 = 10 / len(tmp_df)
		smoothed_hap1 = sm.nonparametric.lowess(endog=tmp_df["logr_hap1"], exog=tmp_df["start"], 
			return_sorted=False, frac = frac1 )
		smoothed_hap2 = sm.nonparametric.lowess(endog=tmp_df["logr_hap2"], exog=tmp_df["start"],
			return_sorted=False, frac = frac1 )
		plt.xticks(rotation=70,fontsize=6)
		ax[j+1].plot(tmp_df["start"],smoothed_hap1,c="red",zorder=1,lw=0.9)
		ax[j+1].plot(tmp_df["start"],smoothed_hap2,c="blue",zorder=1,lw=0.9) #hap2 blue
		plt.xticks(rotation=70,fontsize=6)
		ax[j+1].tick_params(axis="x", labelsize=6,labelrotation=335) 
		ax[j+1].set_xticks(list(np.linspace(start,stop,4)))
		####################### repli
		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax[j+1].add_patch(rect)
	# plt.subplots_adjust(wspace=0, hspace=0)
	# plt.legend(handles=legend,loc=(1.03,0))
	plt.savefig("bouha.switchers.region."+str(chrom) + "." + str(start)+ ".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	# plt.show()
	plt.close()
#####################
# # for plotting whole chroms
# ## use this to plot whole chroms...

# ## switchers RNA and vlincs
# for i in range(len(chromosomes)):
# 	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
# 	tmp_plus = result[(result["chrom"]==chromosomes[i]) & (result["skew_plus"]!=0)]
# 	tmp_minus = result[(result["chrom"]==chromosomes[i]) & (result["skew_minus"]!=0)]
# 	plt.suptitle(chromosomes[i])
# 	ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
# 	ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
# 	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
# 	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
# 	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))

# 	print("plotting vlincs ")
# 	for index,row in vlincs_df[vlincs_df["chrom"]==chromosomes[i]].iterrows():
# 		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
#                      facecolor=row["color"], alpha=0.8,fill=False,hatch="/")
# 		ax[0].add_patch(rect)

# #### repli
# 	for j in range(len(filenames_repli)):
# 		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
# 		ax[j+1].set_ylim([-4,4])
# 		ax[j+1].axhline(y=0,linestyle="--",c="black",lw=0.4)

# 		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		
# 		if len(tmp_df[tmp_df["arm"]=="p"]) <= 10:
# 			frac1 = 1
# 		else:
# 			frac1 = 10 / len(tmp_df[tmp_df["arm"]=="p"])
# 		smoothed_hap1_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_hap1"], exog=tmp_df[tmp_df["arm"]=="p"]["start"], 
# 			return_sorted=False, frac = frac1 )
# 		smoothed_hap2_p = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="p"]["logr_hap2"], exog=tmp_df[tmp_df["arm"]=="p"]["start"],
# 			return_sorted=False, frac = frac1 )
# 		smoothed_hap1_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_hap1"], exog=tmp_df[tmp_df["arm"]=="q"]["start"], 
# 			return_sorted=False, frac = 10/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )
# 		smoothed_hap2_q = sm.nonparametric.lowess(endog=tmp_df[tmp_df["arm"]=="q"]["logr_hap2"], exog=tmp_df[tmp_df["arm"]=="q"]["start"],
# 			return_sorted=False, frac = 10/len(tmp_df[tmp_df["arm"]=="q"]["start"].index) )

# 		if len(tmp_df[tmp_df["arm"]=="p"]) >= 10:
# 			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_hap1_p,c="red",zorder=1,lw=0.8)
# 			ax[j+1].plot(tmp_df[tmp_df["arm"]=="p"]["start"],smoothed_hap2_p,c="blue",zorder=1,lw=0.8)
# 		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_hap1_q,c="red",zorder=1,lw=0.8)
# 		ax[j+1].plot(tmp_df[tmp_df["arm"]=="q"]["start"],smoothed_hap2_q,c="blue",zorder=1,lw=0.8) #hap2 blue
# 		ax[j+1].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))
# 		####################### repli
# 		for index,row in tmp_df[(tmp_df["sig_repli"]=="True")].iterrows():
# 			rect=Rectangle((row["start"],-10),width=row["stop"]-row["start"],height=20,
# 				facecolor=row["repli_color"],alpha=0.6,fill="True")
# 			ax[j+1].add_patch(rect)
# 	# plt.subplots_adjust(wspace=0, hspace=0)
# 	plt.legend(handles=legend,loc=(1.03,0))
# 	plt.savefig("bouha.lcl.switchers.chrom."+str(chromosomes[i])+ ".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)

# 	# plt.show()
# 	plt.close()
# result.to_csv("bouha.lcl.switchers.bed",sep="\t",index=False,header=False)




#####################
# for plotting whole chroms with log2r diff instead of regular repliseq track.
for i in range(len(chromosomes)):
	f,ax = plt.subplots(3,1,figsize=(10,6),sharex=False)
	tmp_plus = result[(result["chrom"]==chromosomes[i]) & (result["skew_plus"]!=0)]
	tmp_minus = result[(result["chrom"]==chromosomes[i]) & (result["skew_minus"]!=0)]
	plt.suptitle(chromosomes[i])
	ax[0].scatter(tmp_plus["start"],tmp_plus["skew_plus"],c=tmp_plus["color"],zorder=1,lw=0.2,edgecolor="black",s=20)
	ax[0].scatter(tmp_minus["start"],tmp_minus["skew_minus"],c=tmp_minus["color"],lw=0.2,zorder=1,edgecolor="black",s=20)
	ax[0].axhline(y=0,linestyle="--",lw=0.4,c="black")
	ax[0].set_xlim([0, chromosome_length[chromosomes[i]]])
	ax[0].set_xticks(np.linspace(0, chromosome_length[chromosomes[i]], 12))

	print("plotting vlincs ") ## all vlincs with somewhat significant skew.
	for index,row in vlincs_df[(vlincs_df["chrom"]==chromosomes[i]) & (vlincs_df["qval"]<=0.01)].iterrows():
		rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
                     facecolor=row["color"], alpha=0.8,fill=False,hatch="/",edgecolor=row["color"])
		ax[0].add_patch(rect)

#### repli
	for j in range(len(filenames_repli)):
		ax[j+1].set_xlim([0, chromosome_length[chromosomes[i]]])
		ax[j+1].set_ylim([0,1.6])

		tmp_df = repli_df[(repli_df["chrom"]==chromosomes[i]) & (repli_df["sample"]==filenames_repli[j])] # chromosomes and sample specific now.
		
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
	# plt.subplots_adjust(wspace=0, hspace=0)
	# plt.legend(handles=legend,loc=(1.03,0))
	plt.savefig("bouha.lcl.switchers.chrom."+str(chromosomes[i])+ ".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)

	# plt.show()
	plt.close()
result.to_csv("bouha.lcl.switchers.bed",sep="\t",index=False,header=False)

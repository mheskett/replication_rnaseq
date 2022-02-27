import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
import seaborn as sns
import scipy.stats
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
import glob
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
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
	df["fdr_pval_plus"]=mt.multipletests(pvals=df["binom_pval_plus"], 
								alpha=0.01,
								method="fdr_bh")[1]
	df["fdr_reject_plus"] =  mt.multipletests(pvals=df["binom_pval_plus"], 
									alpha=0.01,
									method="fdr_bh")[0]
	df["fdr_pval_minus"]=mt.multipletests(pvals=df["binom_pval_minus"], 
								alpha=0.01,
								method="fdr_bh")[1]
	df["fdr_reject_minus"] =  mt.multipletests(pvals=df["binom_pval_minus"], 
									alpha=0.01,
									method="fdr_bh")[0]
	return

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict


def helper_func_plus(x):
	if x["total_plus"]==0:
		return 0
	elif x["hap1_counts_plus"] >= x["hap2_counts_plus"]:
		return x["hap1_counts_plus"]  / x["total_plus"] - 0.5
	else:
		return -x["hap2_counts_plus"]  / x["total_plus"] + 0.5
	return

def helper_func_minus(x):
	if x["total_minus"]==0: # try this for filtering
		return 0
	elif x["hap1_counts_minus"] >= x["hap2_counts_minus"]:
		return x["hap1_counts_minus"]  / x["total_minus"] - 0.5
	else:
		return -x["hap2_counts_minus"]  / x["total_minus"] + 0.5
	return

#########################

rna_files=[
["/Users/mike/replication_rnaseq/all.final.data/gm12878.4.rna.50kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.4x.hg19Aligned.outgm12878.4x.hg19Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/gm12878.5.rna.50kb.bed",
"/Users/mike/replication_rnaseq/all.final.data/vlinc.calls/gm12878.5x.hg19Aligned.outgm12878.5x.hg19Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha4.rna.50kb.bed",
"bouha.trim.4Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.4Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha2.rna.50kb.bed",
"bouha.trim.2Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.2Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha3.rna.50kb.bed",
"bouha.trim.3Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.3Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha10.rna.50kb.bed",
"bouha.trim.10Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.10Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha15.rna.50kb.bed",
"bouha.trim.15Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.15Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
["/Users/mike/replication_rnaseq/all.final.data/bouha13.rna.50kb.bed",
"bouha.trim.13Aligned.samtool.rmdup.plus.all.chrom.allele.counts.haplotype.resolved.counts.bedbouha.trim.13Aligned.out.samtool.rmdup.intergenic.1000.10000.50000.vlinc.discovery.all.bed"],
]
arm_dict = get_arms(cytoband)
for j in range(len(rna_files)):
	df = pd.read_csv(rna_files[j][0],sep="\t",
							names= ["chrom","start","stop","hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str})
	df_vlinc = pd.read_csv(rna_files[j][1],sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
						"hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
						"l1_fraction":float,"hap1_counts":int,"hap2_counts":int,"reject":str})
	##################
	tmp = df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]].replace(".",0)
	tmp = tmp.astype(int)
	df.loc[:,["hap1_counts_plus","hap2_counts_plus","hap1_counts_minus","hap2_counts_minus"]] = tmp
	df = df.set_index(["chrom","start","stop"])
	df = df[df.sum(axis="columns")!=0]
	df = df.reset_index()
	#################
	###################
	add_binom_pval(df)
	print(df)
	################
	## filter the zeros
	df["total_plus"] = df["hap1_counts_plus"] + df["hap2_counts_plus"]
	df["total_minus"] = df["hap1_counts_minus"] + df["hap2_counts_minus"]


	df["skew_plus"] = df.apply(helper_func_plus, axis = 1)
	df["skew_minus"] = df.apply(helper_func_minus, axis = 1)


	#####
	df["arm"] = df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
	df_vlinc["arm"] = df_vlinc.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
	###########
	###########
	color_vector = ["darkred" if row["strand"] == "+" else "darkblue" for index,row in df_vlinc.iterrows() ] # red if hap1 early, blue if hap2 early
	df_vlinc["color"] = color_vector

	# color_vector_rna
	color_vector_plus= []
	color_vector_minus= []
	for index,row in df.iterrows():
		if row["fdr_reject_plus"]==True:
			color_vector_plus +=[ (1,0,0,1) ]
		if row["fdr_reject_plus"]==False:
			color_vector_plus +=[ (1,0,0,0.05) ]
		if row["fdr_reject_minus"]==True:
			color_vector_minus +=[ (0,0,1,1) ]
		if row["fdr_reject_minus"]==False:
			color_vector_minus +=[ (0,0,1,0.05) ]
	df["color_plus"] = color_vector_plus
	df["color_minus"] = color_vector_minus
	## example plots
	plt.scatter(df["total_plus"], abs(df["skew_plus"]),s=8,lw=0.05,color=df["color_plus"],edgecolor="black")
	plt.scatter(df["total_minus"], abs(df["skew_minus"]),s=8,lw=0.05,color=df["color_minus"],edgecolor="black")
	plt.xlim([0,5000])
	plt.close()
	###
	regions=[
	        ["1",186000000,188000000],
	        ["8",2400000,2900000],
	        ["1",13500000,16000000],
	        ["1",56000000,57500000],
	        ["9",0,141213431],
	        ["10",0,135534747],
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

# 	regions = [["1",0,249250621], # whole chroms
# ["2",0,243199373],
# ["3",0,198022430],
# ["4",0,191154276],
# ["5",0,180915260],
# ["6",0,171115067],
# ["7",0,159138663],
# ["8",0,146364022],
# ["9",0,141213431],
# ["10",0,135534747],
# ["11",0,135006516],
# ["12",0,133851895],
# ["13",0,115169878],
# ["14",0,107349540],
# ["15",0,102531392],
# ["16",0,90354753],
# ["17",0,81195210],
# ["18",0,78077248],
# ["19",0,59128983],
# ["20",0,63025520],
# ["21",0,48129895],
# ["22",0,51304566],
# ["X",0,155270560]]

	# regions = [ ## microsatellites
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
	print("plotting")

	size_vector_plus = [float(15+4*np.log2(x+1)) for x in df["total_plus"]]
	size_vector_minus = [float(15+4*np.log2(x+1)) for x in df["total_minus"]]
	df["size_plus"] = size_vector_plus
	df["size_minus"] = size_vector_minus
	for i in range(len(regions)):
		f,ax = plt.subplots(figsize=(10,2)) ## small region is 1.5,2
		# f,ax = plt.subplots(figsize=(8,2)) ## chr arm

		vlincs_tmp = df_vlinc[(df_vlinc["chrom"]==regions[i][0]) & (df_vlinc["start"]>=regions[i][1]-2000000) & 
						(df_vlinc["stop"]<=regions[i][2]+2000000) & (df_vlinc["reject"]=="True")]
		########################## vlincs
		ax.set_ylim([-0.52,0.52])
		ax.set_xlim(regions[i][1],regions[i][2])
		ax.axhline(y=0,linestyle="--",c="black",lw=0.2)
		for index,row in vlincs_tmp.iterrows():
			rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
		                 facecolor=row["color"], edgecolor=row["color"],hatch="/",alpha=0.9,fill=False) ## plot vlincs as rectangles
			ax.add_patch(rect)
		#### RNA
		filtered_plus = df[df["total_plus"]>=15]
		filtered_minus = df[df["total_minus"]>=15]
		rna_tmp_plus = filtered_plus[(filtered_plus["chrom"]==regions[i][0]) & (filtered_plus["start"]>=regions[i][1]-500000) & (filtered_plus["stop"]<=regions[i][2]+500000)]
		rna_tmp_minus = filtered_minus[(filtered_minus["chrom"]==regions[i][0]) & (filtered_minus["start"]>=regions[i][1]-500000) & (filtered_minus["stop"]<=regions[i][2]+500000)]
		print(rna_tmp_plus)
		ax.scatter(rna_tmp_plus["start"],
					rna_tmp_plus["skew_plus"],c=rna_tmp_plus["color_plus"],lw=0.1,zorder=1,edgecolor="black",s=15) # 15+rna_tmp_plus["total_plus"]/12
		ax.scatter(rna_tmp_minus["start"],
					rna_tmp_minus["skew_minus"],c=rna_tmp_minus["color_minus"],lw=0.1,zorder=1,edgecolor="black",s=15) # 15+rna_tmp_plus["total_plus"]/12
		ax.set_xticks(np.linspace(regions[i][1],regions[i][2],12)) #### number of ticks
		ax.set_yticks([-0.5,-.4,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5])
		plt.xticks(rotation = 315) # Rotates X-Axis Ticks by 45-degrees
		ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
		plt.savefig(os.path.basename(rna_files[j][0])[0:15]+"-region"+str(regions[i][0])+"-"+str(regions[i][1])+"-"+str(regions[i][2])+".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
		plt.close()
import os
import csv
import numpy as np
import pandas as pd
import re
import seaborn as sns
import scipy.stats
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns
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
	return

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

regions=[

			# ["1",187416000,187459201],
			# ["3",174933729,174972546],
			# ["6",40074887,40115741],
			# ["6",73394774,73431597],
			# ["6",77914535,77951564],
			# ["6",120454710,120497591],
			["6",130585933,130623577], ### for paper
			["6",130265000,131160000], ### closer zoo mfor paper
			# ["6",141108314,141151646],
			# ["6",141184172,141225657],
			# ["8",2586332,2625241],
			# ["8",22568321,22607943],
			# ["9",22765327,22807785],
			# ["9",23643340,23683461],
			# ["9",24038517,24077283],
			# ["9",29802888,29843986],
			# ["15",47102705,47144857],
			# ["15",54386886,54427639],
			# ["15",54629690,54670049],
			["15",92002424,92044668], ## for paper
			["15",91900000,92450000], ## closer zoom for paper
			# ["15",96464802,96501683],
			# ["15",97096204,97134188],

         ]

rna_files=[
["/Users/mike/replication_rnaseq/scripts/gm12878.rep1.protein.coding.all.counts.bed","/Users/mike/replication_rnaseq/all.final.data/gm12878.rep1.vlincs.all.bed"]]

arm_dict = get_arms(cytoband)
for j in range(len(rna_files)):
	df = pd.read_csv(rna_files[j][0],sep="\t",
							names= ["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
	df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]

	df["skew"] = df.apply(helper_func, axis = 1)
	df=df[df["total_reads"]>=10]
	add_binom_pval(df)
#################
	df_vlinc = pd.read_csv(rna_files[j][1],sep="\t",
							names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction",
							"hap1_counts","hap2_counts"],
							dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,
							"l1_fraction":float,"hap1_counts":int,"hap2_counts":int,"reject":str})
	df_vlinc["total_reads"] = df_vlinc["hap1_counts"] + df_vlinc["hap2_counts"]
	df_vlinc["skew"] = df_vlinc.apply(helper_func, axis = 1)
	df_vlinc["sample"] = os.path.basename(rna_files[j][1])[0:15]
	df_vlinc["informative_reads_per_kb"] = df_vlinc["total_reads"] / ((df_vlinc["stop"] - df_vlinc["start"])  / 1000)
	df_vlinc = df_vlinc[df_vlinc["total_reads"]>=15]
	add_binom_pval(df_vlinc)
##############
	color_vector_df= []
	for index,row in df.iterrows():
		if ((row["fdr_reject"]==True) & (abs(row["skew"])>=.1)):
			if row["strand"] == "+":
				color_vector_df +=[ (1,0,0,1) ]
			if row["strand"] == "-":
				color_vector_df +=[ (0,0,1,1) ]	
		if row["fdr_reject"]==False or ((row["fdr_reject"]==True) and (abs(row["skew"])<.1)):
			if row["strand"] == "+":
				color_vector_df +=[ (1,0,0,.2) ]
			if row["strand"] == "-":
				color_vector_df +=[ (0,0,1,.2) ]	
	color_vector_vlinc = []
	for index,row in df_vlinc.iterrows():
		if ((row["fdr_reject"]==True) & (abs(row["skew"])>=.1)):
			if row["strand"] == "+":
				color_vector_vlinc +=[ (1,0,0,1) ]
			if row["strand"] == "-":
				color_vector_vlinc +=[ (0,0,1,1) ]	
		if row["fdr_reject"]==False or ((row["fdr_reject"]==True) and (abs(row["skew"])<.1)):
			if row["strand"] == "+":
				color_vector_vlinc +=[ (1,0,0,.2) ]
			if row["strand"] == "-":
				color_vector_vlinc +=[ (0,0,1,.2) ]		
	df["color"] = color_vector_df
	df_vlinc["color"] = color_vector_vlinc
##################
# protein coding
	for i in range(len(regions)):
		f,ax = plt.subplots(figsize=(4,1.5)) ## small region is 1.5,2
		# f,ax = plt.subplots(figsize=(8,2)) ## chr arm
		chrom = regions[i][0]
		start = regions[i][1]
		stop = regions[i][2]
		plt.suptitle(chrom)
		tmp = df[(df["chrom"]==chrom) & (df["start"]>=start-500000) & 
						(df["stop"]<=stop+500000)]
		print(tmp)
		# ax.scatter(tmp["start"],tmp["skew"],c=tmp["color"],zorder=1,lw=0.2,edgecolor="black",s=40)
		for index,row in tmp.iterrows():
			rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
		                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
			ax.add_patch(rect)
			x = range(start,stop)
			y = range(row["start"],row["stop"])
			# if range(max(x[0], y[0]), min(x[-1], y[-1])+1):
			# 	ax.text(row["start"]+10000,row["skew"]-0.1,row["name"][0:15],size=5) ## to get the gene names

		ax.axhline(y=0,linestyle="--",lw=0.4,c="black")
		ax.set_xlim([start, stop])
		ax.set_ylim([-0.6,0.6])
		ax.set_xticks(np.linspace(start,stop, 6))
		ax.set_yticks([-0.5,-.25,0,.25,.5])
		print(np.linspace(start,stop,4))
	## vlincs########################
		vlincs_tmp = df_vlinc[(df_vlinc["chrom"]==chrom) & (df_vlinc["start"]>=start-2000000) & 
						(df_vlinc["stop"]<=stop+2000000)]
		########################## vlincs
		ax.axhline(y=0,linestyle="--",c="black",lw=0.2)
		print(vlincs_tmp)
		for index,row in vlincs_tmp.iterrows():
			rect=Rectangle((row["start"], row["skew"]-.05), width=row["stop"]-row["start"], height=0.1,
		                 facecolor=row["color"], edgecolor=row["color"],hatch="/",fill=False) ## plot vlincs as rectangles
			ax.add_patch(rect)
		plt.savefig("gm12878.rep1.protein.vlinc."+str(chrom)+"."+str(start)+"."+str(stop)+".png",
		dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
		plt.close()




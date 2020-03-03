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

def get_windows(peaks, allele_counts): 
	
	a = pybedtools.BedTool(peaks) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(allele_counts) ## read counts at alleles

	a = a.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a = a , b = b , c = [6,7] , o = ["sum","sum"], g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	### invalid int for literal...
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int, "hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df

def add_binom_pval(df):
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],
							row["hap1_reads"]+row["hap2_reads"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)

	return

ratios = [249250621/249250621,
	243199373/249250621,
	198022430/249250621,
	191154276/249250621,
	180915260/249250621,
	171115067/249250621,
	159138663/249250621,
	146364022/249250621,
	141213431/249250621,
	135534747/249250621,
	135006516/249250621,
	133851895/249250621,
	115169878/249250621,
	107349540/249250621,
	102531392/249250621,
	90354753/249250621,
	81195210/249250621,
	78077248/249250621,
	59128983/249250621,
	63025520/249250621,
	48129895/249250621,
	51304566/249250621,
	155270560/249250621]

lengths = [249250621,
	243199373,
	198022430,
	191154276,
	180915260,
	171115067,
	159138663,
	146364022,
	141213431,
	135534747,
	135006516,
	133851895,
	115169878,
	107349540,
	102531392,
	90354753,
	81195210,
	78077248,
	59128983,
	63025520,
	48129895,
	51304566,
	155270560]

centromere = {"1":124535434,
				"2":95326171,
				"3":93504854,
				"4":52660117,
				"5":49405641,
				"6":61830166,
				"7":61054331,
				"8":46838887,
				"9":50367679,
				"X":61632012,
				"Y":13104553,
				"10":42254935,
				"11":54644205,
				"12":37856694,
				"13":19000000,
				"14":19000000,
				"15":20000000,
				"16":38335801,
				"17":25263006,
				"18":18460898,
				"19":27681782,
				"20":29369569,
				"21":14288129,
				"22":16000000}
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/data/cytoband.nochr.hg19.bed",
							sep="\t",
							names =["chrom","start","stop","arm","band"])

file_list = ["H2AFZ-human.significant.monoallelic.peaks.bed",
			"H3K9me3-human.significant.monoallelic.peaks.bed",
			"H3K27ac-human.significant.monoallelic.peaks.bed",  
			"H3K4me1-human.significant.monoallelic.peaks.bed",  
			"H3K79me2-human.significant.monoallelic.peaks.bed",
			"H4K20me1-human.significant.monoallelic.peaks.bed",
			"H3K27me3-human.significant.monoallelic.peaks.bed", 
			"H3K4me2-human.significant.monoallelic.peaks.bed",
			"H3K9ac-human.significant.monoallelic.peaks.bed"]

color_map = {"H2AFZ-human":"red",
			"H3K9me3-human":"orange",
			"H3K27ac-human":"yellow",
			"H3K4me1-human":"green",
			"H3K79me2-human":"blue",
			"H4K20me1-human":"purple",
			"H3K27me3-human":"black",
			"H3K4me2-human":"pink",
			"H3K9ac-human":"gold"}

all_chip_dfs = [pd.read_csv(x,header=0,sep="\t",
		names = ["chrom","start","stop","hap1_reads","hap2_reads","pval"],
		dtype = {"chrom":str,"start":int,"stop":int,"hap1_reads":int,"hap2_reads":int,"pval":float})
		for x in file_list]

for i in range(len(all_chip_dfs)):
	all_chip_dfs[i]["label"] = file_list[i]

all_chip = pd.concat(all_chip_dfs,axis=0)

all_chip["label"] = all_chip["label"].str.split(".",expand=True)[0]
df_expression = pd.read_csv("gm12878.rep1.hg19Aligned.out.tiling.all.bed",sep="\t",
	names = ["chrom","start","stop","hap1_reads","hap2_reads","strand","pval","qval","reject","total_reads","skew"],
	dtype = {"chrom":str,"start":int,"stop":int,"hap1_reads":int,"hap2_reads":int,"strand":str,"pval":float,"qval":float,"reject":str,"total_reads":int,"skew":float})

all_chip.loc[:,"hap_label"] = all_chip.apply(lambda x: "red" if x["hap1_reads"] > x["hap2_reads"] else "green", axis = 1)
all_chip.loc[:,"total_reads"] = all_chip["hap1_reads"] + all_chip["hap2_reads"]
all_chip.loc[:,"skew"] = all_chip.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
															(-x["hap2_reads"]  / x["total_reads"])	+ 0.5, axis = 1)
all_chip.sort_values(by=["chrom","start"]).to_csv("gm12878.all.chip.significant.monoallelic.peaks.bed", sep="\t", header=None, index=None)
all_chip.loc[:,"color"] = [color_map[x] for x in all_chip["label"]]
for i in range(len(chromosomes)):
	f,ax = plt.subplots(nrows=2, ncols=1,figsize=(12,2))
	ax[0].scatter(df_expression[df_expression["chrom"]==chromosomes[i]]["start"], df_expression[df_expression["chrom"]==chromosomes[i]]["skew"],
		s=10,lw=0.2,edgecolor="black")
	tmp = all_chip[(all_chip["chrom"]==chromosomes[i])]
	ax[1].scatter(tmp["start"], 
		tmp["skew"],
		color=tmp["color"],
		s=10, lw=0.2, edgecolor="black")
	ax[1].legend()
	plt.show()
	plt.close()

# 	plt.show()
# 	plt.close()
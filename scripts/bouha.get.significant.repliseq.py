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
"""
seq directions ( do this separately for early and late)

1. take allele counts bed file and sum over 100kb windows that are tiled by 25kb -- 
2. get autosome wide distribution of logR values
3. only keep top 10% most asynchronous
4. Merge the segments, requiring at least 300kb contiguous with allowable gap of 50kb
4. Merge again, allowing for a 100kb gap.
5. Compare overlap between early and late significant regions. get p-values
"""
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

ratios=[249250621/249250621,
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

def library_size_normalization(df):
	new_df = df.copy()
	total_reads = new_df["hap1_counts"].sum() + new_df["hap2_counts"].sum()
	new_df["hap1_lsm"] = (new_df["hap1_counts"] / total_reads) * 10**6
	new_df["hap2_lsm"] = (new_df["hap2_counts"] / total_reads) * 10**6

	return new_df


def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

def add_binom_pval(df):
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],
							row["hap1_reads"]+row["hap2_reads"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)

	df["fdr_pval"]=mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")[1]
	df["fdr_reject"] =  mt.multipletests(pvals=df["binom_pval"], 
									alpha=0.01,
									method="fdr_bh")[0]
	return
		
def get_tiling_windows(read_counts_file, size):
	a = pybedtools.BedTool()
	a = a.window_maker(w=size, s=size/4, g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = pybedtools.BedTool(read_counts_file)
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	c = pybedtools.BedTool()

	window_read_counts = c.map(a=a, b=b, c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df
arm_dict = get_arms(cytoband)

### point of this script will be to
### 1. plot regular repliseq for multiple files including library size normalization
### 2. figure out allele specific differences 
### 3.
### 4.

df_windows = pd.read_csv("4e.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
df2_windows = pd.read_csv("4l.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","hap1_counts","hap2_counts"],
						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
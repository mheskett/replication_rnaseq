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

	# df["fdr_pval"]=mt.multipletests(pvals=df["binom_pval"], 
	# 							alpha=0.01,
	# 							method="fdr_bh")[1]
	# df["fdr_reject"] =  mt.multipletests(pvals=df["binom_pval"], 
	# 								alpha=0.01,
	# 								method="fdr_bh")[0]
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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	parser.add_argument("--peaks",
		type=str,
		metavar="[df of allele counts plus]",
		required=True,
		help="in bed format")
	parser.add_argument("--allele_counts",
		type=str,
		metavar="[allele specific counts]",
		required=True,
		help="in bed format")
	arguments = parser.parse_args()
	df = get_windows(peaks = arguments.peaks, 
					allele_counts = arguments.allele_counts)
	add_binom_pval(df)
	df_significant = df[df["binom_pval"] <= 0.05]

	df_significant.to_csv(os.path.basename(arguments.peaks.replace(".merged.peaks.bed",'.significant.monoallelic.peaks.bed')),
							sep="\t",index=None)
	print(df_significant)




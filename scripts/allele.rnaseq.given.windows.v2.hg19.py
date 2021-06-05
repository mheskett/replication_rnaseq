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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import statsmodels.stats.multitest as mt

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
		
def get_windows(window_file, read_counts_file, is_file=True): 
	
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles

	a = a.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=a,b=b,c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	### invalid int for literal...
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "name", "score", "strand", 
										"fraction_l1","hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"name":str, "score":float, "strand":str, "fraction_l1":float,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df

def get_windows_from_merged(window_file, read_counts_file): 
	a = window_file
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles

	a = a.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=a,b=b,c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	### invalid int for literal...
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", 
										"hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df


def get_tiling_windows(read_counts_file, size):
	a = pybedtools.BedTool()
	a = a.window_maker(w=size, s=size, g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
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
# hg19
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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	parser.add_argument("--dfplus",
		type=str,
		metavar="[df of allele counts plus]",
		required=True,
		help="in bed format")
	parser.add_argument("--dfminus",
		type=str,
		metavar="[df of allele counts minus]",
		required=True,
		help="in bed format")
	parser.add_argument("--window_file_plus",
		type=str,
		metavar="[bed file of windows]",
		required=False,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--window_file_minus",
		type=str,
		metavar="[bed file of windows]",
		required=False,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--out_directory",
		type=str,
		metavar="[out directory]",
		required=True,
		help="full path to output results")
	parser.add_argument("--repli_seq",
		type=str,
		metavar="[repliseq_files]",
		required=False,
		help="includ repliseq files if you have them")
	parser.add_argument("--tiling_windows",
		required=False,
		action="store_true")
	parser.add_argument("--single_chrom",
		type=str,
		metavar="[single_chrom]",
		required=False,
		help="plot a single chromosome in addition to all of them")

	###############
	# Main program
	###############

	arguments = parser.parse_args()
	############################
	###########################
	############################

	#####
	## tiling window specific program
	if arguments.tiling_windows:

	### Old method used this equation for skew: = df_plus.apply(lambda x: np.log2(x["hap1_reads"]  / x["total_reads"] / 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
	##														-np.log2(x["hap2_reads"]  / x["total_reads"] / 0.5)	, axis = 1)

	# new method just uses percentage

		df_plus = get_tiling_windows(read_counts_file=arguments.dfplus, size = 50000)
		df_minus = get_tiling_windows(read_counts_file=arguments.dfminus, size = 50000)
		df_plus.loc[:,"strand"] = "+"
		df_minus.loc[:,"strand"] = "-"
		print("tiling windows selected")
		add_binom_pval(df_plus)
		add_binom_pval(df_minus)
		## treat plus and minus separately still
		df_plus["total_reads"] = df_plus["hap1_reads"] + df_plus["hap2_reads"]
		df_plus = df_plus[df_plus["total_reads"]>=20] # or do at least min 1 read per kb
		df_plus.loc[:,"skew"] = df_plus.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
															(-x["hap2_reads"]  / x["total_reads"] + 0.5), axis = 1)
		df_plus_skewed = df_plus[(df_plus["binom_pval"]<=10**-6) & (abs(df_plus["skew"])>=0.1)]
		df_plus_skewed_merged = pybedtools.BedTool.from_dataframe(df_plus_skewed).merge(s=True) ## crucial merging step right here.....
		### do minus
		df_minus["total_reads"] = df_minus["hap1_reads"] + df_minus["hap2_reads"]
		df_minus = df_minus[df_minus["total_reads"]>=20] # or do at least min 1 read per kb
		df_minus.loc[:,"skew"] = df_minus.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
															(-x["hap2_reads"]  / x["total_reads"] + 0.5)	, axis = 1)
		df_minus_skewed = df_minus[(df_minus["binom_pval"]<=10**-6) & (abs(df_minus["skew"])>=0.1)]
		df_minus_skewed_merged = pybedtools.BedTool.from_dataframe(df_minus_skewed).merge(s=True)
		df_combined = pd.concat([df_plus,df_minus]) # non merged, but still min read filtered. use this for plotting. 
		# df_combined = df_combined[df_combined["chrom"]!="X"]

		## then put it back in to get windows, recount the number of reads, combine, add binom pvals.
		## basically this allows you to do tiling windows, then merge
		df_plus = get_windows_from_merged(window_file = df_plus_skewed_merged, read_counts_file=arguments.dfplus )
		df_minus = get_windows_from_merged(window_file = df_minus_skewed_merged, read_counts_file=arguments.dfminus)
		df_plus.loc[:,"strand"] = "+"
		df_minus.loc[:,"strand"] = "-"
		df_combined_merged = pd.concat([df_plus, df_minus])
		add_binom_pval(df_combined_merged)
		df_combined_merged["total_reads"] = df_combined_merged["hap1_reads"] + df_combined_merged["hap2_reads"]
		df_combined_merged.loc[:,"skew"] = df_combined_merged.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
															(-x["hap2_reads"]  / x["total_reads"] + 0.5)	, axis = 1)
		df_combined_merged = df_combined_merged[(df_combined_merged["binom_pval"]<=10**-6) & (abs(df_combined_merged["skew"])>=0.1)] # use this for actual stats on 
		# non overlapping mono allelic regions

		# df_combined_merged = df_combined_merged[df_combined_merged["chrom"]!="X"]

	else:
		df_plus = get_windows(window_file = arguments.window_file_plus,
							read_counts_file=arguments.dfplus)
		df_minus = get_windows(window_file = arguments.window_file_minus,
							read_counts_file=arguments.dfminus)
		df_combined = pd.concat([df_plus, df_minus])
		add_binom_pval(df_combined)
		df_combined["total_reads"] = df_combined["hap1_reads"] + df_combined["hap2_reads"]

		df_combined = df_combined[df_combined["total_reads"]>=20] # or do at least min 1 read per kb
		df_combined.loc[:,"skew"] = df_combined.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
															(-x["hap2_reads"]  / x["total_reads"])	+ 0.5, axis = 1)

	# df_combined = df_plus[df_plus["chrom"]!="X"]

	## hap1 is positive skew hap2 is negative skew
	print(df_combined)
	f, ax = plt.subplots(1, len(chromosomes), sharex=False,
								sharey=False,
								figsize=(15,1),
								gridspec_kw={'width_ratios': ratios})

	for i in range(len(chromosomes)):
		# if theres no data then just set the formatting and skip to next one
		if chromosomes[i] not in list(df_combined.chrom):
			ax[i].set_yticks([])
			ax[i].set_xticks([])
			ax[i].margins(x=0,y=0)

			# formatting
			ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
			ax[i].axvline(x=int(centromere[chromosomes[i]])/10**6, linestyle = "--", lw = 0.5,color="black")
			ax[i].set_xlim([0,lengths[i]/10**6])
			continue
		### still need two to achieve different cmaps for each point
		df_combined_plus = df_combined[(df_combined["strand"]=="+") & (df_combined["chrom"]==chromosomes[i]) ]
		ax[i].scatter(df_combined_plus["start"]/10**6,
		df_combined_plus["skew"],
		s=6,
		lw=0.1,
		edgecolor="black",
		c=-np.log10(df_combined_plus["binom_pval"]),
		vmin=4,
		vmax=30,
		cmap="Greens",
		alpha=0.6)
		#####
		df_combined_minus = df_combined[(df_combined["strand"]=="-") & (df_combined["chrom"]==chromosomes[i]) ]
		ax[i].scatter(df_combined_minus["start"]/10**6,
		df_combined_minus["skew"],
		s=6,
		lw=0.1,
		edgecolor="black",
		c=-np.log10(df_combined_minus["binom_pval"]),
		vmin=4,
		vmax=30,
		cmap="Reds",
		alpha=0.6)

		### formatting	
		ax[i].set_yticks([])
		ax[i].set_xticks([])
		ax[i].margins(x=0,y=0)
		ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
		ax[i].axvline(x=int(centromere[chromosomes[i]])/10**6, linestyle = "--", lw = 0.5,color="black")
		ax[i].set_xlim([0,lengths[i]/10**6])
		ax[i].set_ylim([-0.5, 0.5])
	## make some chromosomes gray
	for i in range(len(chromosomes)):
		if chromosomes[i] in gray_chromosomes:
		    ax[i].axvspan(xmin=0, xmax=lengths[i]/10**6, ymin=0, ymax=1,
		     alpha=0.2,facecolor="gray")
	plt.subplots_adjust(wspace=0, hspace=0)

	if arguments.tiling_windows:
		out_string = arguments.out_directory+os.path.basename(arguments.dfplus.replace(".plus.overlap.platinum.haplotypes.bed",""))+"."+"tiling"
		df_combined_merged.to_csv(out_string + ".skewed.bed", sep="\t", index=None, header=None)

	else: 
		out_string = arguments.out_directory+os.path.basename(arguments.dfplus.replace(".plus.overlap.platinum.haplotypes.bed",""))+os.path.basename(arguments.window_file_plus).replace(".bed","").replace(".plus","")
		df_combined[(df_combined["binom_pval"]<=10**-6) & (abs(df_combined["skew"])>=0.1)].to_csv(out_string + ".skewed.bed", sep="\t", index=None, header=None)

	plt.savefig(out_string+".png", dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)

	# combined_df = pd.concat([df_plus, df_minus])
	df_combined.to_csv(out_string + ".all.bed", sep="\t", index=None, header=None)
	if not arguments.tiling_windows:
		### dont need to make genome browser files for just random windows
		### output files with proper filenames
		browser_df = df_combined[df_combined["binom_pval"] <= 0.000001].loc[:,:"strand"]
		browser_df["chrom"] = "chr"+browser_df["chrom"].astype(str)
		browser_df.to_csv(out_string + ".browser.bed", sep="\t",index=None,header=None)
	
	####################
	####################
	####################
	#################### 
	# single chrom
	if arguments.single_chrom:
		f, ax = plt.subplots(1, 1, sharex=False,
									sharey=False,
									figsize=(10,1))
		chrom = str(arguments.single_chrom)
		### still need two to achieve different cmaps for each point
		df_combined_plus = df_combined[(df_combined["strand"]=="+") & (df_combined["chrom"]==chrom) ]
		ax.scatter(df_combined_plus["start"]/10**6,
		df_combined_plus["skew"],
		s=14,
		lw=0.1,
		edgecolor="black",
		c=-np.log10(df_combined_plus["binom_pval"]),
		vmin=4,
		vmax=30,
		cmap="Greens",
		alpha=0.6)
		#########
		df_combined_minus = df_combined[(df_combined["strand"]=="-") & (df_combined["chrom"]==chrom) ]
		ax.scatter(df_combined_minus["start"]/10**6,
		df_combined_minus["skew"],
		s=14,
		lw=0.1,
		edgecolor="black",
		c=-np.log10(df_combined_minus["binom_pval"]),
		vmin=4,
		vmax=30,
		cmap="Reds",
		alpha=0.6)
		### formatting	
		ax.set_yticks([])
		ax.margins(x=0,y=0)
		ax.set(xlabel=chrom) # x axis labels or no
		ax.axvline(x=int(centromere[chrom]) / 10**6, linestyle = "--", lw = 0.5,color="black")
		ax.set_xlim([0, lengths[i] / 10**6] )
		ax.set_ylim([-0.5, 0.5])
		if chrom=="X":
			chrom_index=22
		else:
			chrom_index=int(chrom) - 1
		plt.savefig(out_string+chrom+".only"+".png", dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)

	####################
	####################
	####################
	####################
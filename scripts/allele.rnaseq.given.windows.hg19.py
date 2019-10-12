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

# hg19
ratios=[249250621/249250621,243199373/249250621,198022430/249250621,191154276/249250621,180915260/249250621,
171115067/249250621,159138663/249250621,146364022/249250621,141213431/249250621,135534747/249250621,
135006516/249250621,133851895/249250621,115169878/249250621,107349540/249250621,102531392/249250621,
90354753/249250621,81195210/249250621,78077248/249250621,59128983/249250621,63025520/249250621,
48129895/249250621,51304566/249250621,155270560/249250621]

lengths = [249250621,243199373,198022430,191154276,180915260,
171115067,159138663,146364022,141213431,135534747,
135006516,133851895,115169878,107349540,102531392,
90354753,81195210,78077248,59128983,63025520,
48129895,51304566,155270560]
# hg19
centromere = {"1":"124535434",
"2":"95326171",
"3":"93504854",
"4":"52660117",
"5":"49405641",
"6":"61830166",
"7":"61054331",
"8":"46838887",
"9":"50367679",
"X":"61632012",
"Y":"13104553",
"10":"42254935",
"11":"54644205",
"12":"37856694",
"13":"19000000",
"14":"19000000",
"15":"20000000",
"16":"38335801",
"17":"25263006",
"18":"18460898",
"19":"27681782",
"20":"29369569",
"21":"14288129",
"22":"16000000"}
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]

def get_windows(window_file, read_counts_file): #
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=window_file,b=b,c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "name", "score", "strand_of_window",
										"hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"name":str, "score":float, "strand_of_window":str,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df

### hg38
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
		required=True,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--window_file_minus",
		type=str,
		metavar="[bed file of windows]",
		required=True,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--out_directory",
		type=str,
		metavar="[out directory]",
		required=True,
		help="full path to output results")
########## funcs
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

	def color_maker(x):
		result = [0,0,0,0]
		if x["strand_of_window"]=="-":
			result[0]=1 # minus strand red
		if x["strand_of_window"]=="+":
			result[1]=1 # plus strand green
		if x["strand_of_window"]==".":
			result[2]=1 ## no strand info is blue
		if x["fdr_reject"]==True:
			result[3]=1
		if x["fdr_reject"]==False:
			result[3]=0.1
		return result
	def marker_size(x):
		if x["total_reads"] <=100:
			return 10
		if 100 < x["total_reads"] <= 1000:
			return 40
		if 1000 < x["total_reads"]:
			return 75
###############
# Main program
###############

	arguments = parser.parse_args()
	df_plus = get_windows(window_file = arguments.window_file_plus,
		read_counts_file=arguments.dfplus)
	df_minus = get_windows(window_file = arguments.window_file_minus,
		read_counts_file=arguments.dfminus)
	add_binom_pval(df_plus)
	add_binom_pval(df_minus)

	f, ax = plt.subplots(1, len(chromosomes), sharex=False,
											sharey=False,
											figsize=(15,1),
											gridspec_kw={'width_ratios': ratios})

	for i in range(len(chromosomes)):
		if chromosomes[i] not in list(df_plus.chrom):
			ax[i].set_yticks([])
			ax[i].set_xticks([])
			ax[i].margins(x=0,y=0)

			# formatting
			ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
			ax[i].axvline(x=int(centromere[chromosomes[i]])/10**6, linestyle = "--", lw = 0.5,color="black")
			ax[i].set_xlim([0,lengths[i]/10**6])
			continue
		## clean up this bullcrap
		#### make a df for each
		### hap1 plus
		hap1_plus = df_plus[(df_plus["hap1_reads"] >= df_plus["hap2_reads"]) & (df_plus["chrom"]==chromosomes[i]) ]
		if len(hap1_plus)>0:
			hap1_plus["total_reads"] = hap1_plus["hap1_reads"] + hap1_plus["hap2_reads"]
			hap1_plus["color"] = hap1_plus.apply(color_maker,axis=1)
			hap1_plus["size"] = hap1_plus.apply(marker_size,axis= 1 )
			ax[i].scatter(hap1_plus["start"]/10**6,
					[-np.log10(x) if ( x >= 10**-12 ) else 12 for x in hap1_plus["binom_pval"]],
			 		s=hap1_plus["size"], # this scaling sucks still
					lw=0.1,
					edgecolor="black",
					c=hap1_plus["color"])
			### text annotation
			# for index,row in hap1_plus[hap1_plus["fdr_reject"]==True].iterrows():
			# 	ax[i].annotate(s=row["name"],
			# 	xy=(row["start"]/10**6, -np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12)).set_fontsize(6)
		### hap1 minus
		hap1_minus = df_minus[(df_minus["hap1_reads"] >= df_minus["hap2_reads"]) & (df_minus["chrom"]==chromosomes[i]) ]
		if len(hap1_minus)>0:
			hap1_minus["total_reads"] = hap1_minus["hap1_reads"] + hap1_minus["hap2_reads"]
			hap1_minus["color"] = hap1_minus.apply(color_maker,axis=1)
			hap1_minus["size"] = hap1_minus.apply(marker_size,axis= 1 )
			ax[i].scatter(hap1_minus["start"]/10**6,
					[-np.log10(x) if ( x >= 10**-12 ) else 12 for x in hap1_minus["binom_pval"]],
			 		s=hap1_minus["size"], # this scaling sucks still
					lw=0.1,
					edgecolor="black",
					c=hap1_minus["color"])
			### text annotation
			# for index,row in hap1_minus[hap1_minus["fdr_reject"]==True].iterrows():
			# 	ax[i].annotate(s=row["name"],
			# 	xy=(row["start"]/10**6, -np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12)).set_fontsize(6)
		### this breaks if the df is empty...thecrappy fix is to assign windows where hap1 reads = hap2 reads to both hap1 and hap2.
		### it may double plot them, but they will be non significant anyways so who cares?
		#### hap2 plus
		hap2_plus = df_plus[(df_plus["hap1_reads"] <= df_plus["hap2_reads"]) & (df_plus["chrom"]==chromosomes[i]) ]
		if len(hap2_plus)>0:
			hap2_plus["total_reads"] = hap2_plus["hap1_reads"] + hap2_plus["hap2_reads"]
			hap2_plus["color"] = hap2_plus.apply(color_maker,axis=1)
			hap2_plus["size"] = hap2_plus.apply(marker_size,axis=1)
			ax[i].scatter(hap2_plus["start"]/10**6,
					[np.log10(x) if ( x >= 10**-12 ) else -12 for x in hap2_plus["binom_pval"]],
			 		s=hap2_plus["size"],
					lw=0.1,
					edgecolor="black",
					c=hap2_plus["color"])
			### annotation by text. doesnt show up well.
			# for index,row in hap2_plus[hap2_plus["fdr_reject"]==True].iterrows():
			# 	ax[i].annotate(s=row["name"],
			# 	xy=(row["start"]/10**6, -1*(-np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12))).set_fontsize(6)
		### hap2 minus
		hap2_minus = df_minus[(df_minus["hap1_reads"] <= df_minus["hap2_reads"]) & (df_minus["chrom"]==chromosomes[i]) ]
		if len(hap2_minus)>0:
			hap2_minus["total_reads"] = hap2_minus["hap1_reads"] + hap2_minus["hap2_reads"]
			hap2_minus["color"] = hap2_minus.apply(color_maker,axis=1)
			hap2_minus["size"] = hap2_minus.apply(marker_size,axis=1)
			ax[i].scatter(hap2_minus["start"]/10**6,
					[np.log10(x) if ( x >= 10**-12 ) else -12 for x in hap2_minus["binom_pval"]],
			 		s=hap2_minus["size"],
					lw=0.1,
					edgecolor="black",
					c=hap2_minus["color"])
			### annotation by text. doesnt show up well.
			# for index,row in hap2_minus[hap2_minus["fdr_reject"]==True].iterrows():
			# 	ax[i].annotate(s=row["name"],
			# 	xy=(row["start"]/10**6, -1*(-np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12))).set_fontsize(6)

		ax[i].set_yticks([])
		ax[i].set_xticks([])
		ax[i].margins(x=0,y=0)

		# formatting
		ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
		ax[i].axvline(x=int(centromere[chromosomes[i]])/10**6, linestyle = "--", lw = 0.5,color="black")
		ax[i].set_xlim([0,lengths[i]/10**6])
		#ax[i].set_yticks([-12,-8,-4,0,4,8,12])
		#plt.grid(True)
		ax[i].set_ylim([-14,14]) # for all plots
	for i in range(len(chromosomes)):
		if chromosomes[i] in gray_chromosomes:
		    ax[i].axvspan(xmin=0, xmax=lengths[i]/10**6, ymin=0, ymax=1,
		     alpha=0.2,facecolor="gray") 
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.show()
	# plt.savefig(os.path.basename(arguments.dfp.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".png",dpi=400,transparent=True,bbox_inches='tight',pad_inches = 0)

	# ### output files with proper filenames
	# df.to_csv(os.path.basename(arguments.df.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".bed",sep="\t",index=None,header=None)
	# browser_df = df[df["fdr_reject"]==True].loc[:,:"strand_of_window"]
	# browser_df["chrom"] = "chr"+browser_df["chrom"].astype(str)
	# browser_df.to_csv(os.path.basename(arguments.df.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".browser.bed",sep="\t",index=None,header=None)
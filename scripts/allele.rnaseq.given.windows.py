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
## hg38
lengths = [248956422,#1
242193529,#2
198295559,#3
190214555,#4
181538259,# etc
170805979,
159345973,
145138636,
138394717,
133797422,
135086622,
133275309,
114364328,
107043718,
101991189,
90338345,
83257441,
80373285,
58617616,
64444167,
46709983,
50818468,
156040895] #X
#hg38
ratios = [x/lengths[0] for x in lengths]
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]


centromere={"1":123400000, # hg38
"2":93900000,
"3":90900000,
"4":50000000,
"5":48750000,
"6":60550000,
"7":60100000,
"8":45200000,
"9":43850000,
"10":39800000,
"11":53400000,
"12":35500000,
"13":17700000,
"14":17150000,
"15":19000000,
"16":36850000,
"17":25050000,
"18":18450000,
"19":26150000,
"20":28050000,
"21":11950000,
"22":15550000,
"X":60950000}

def get_windows(window_file, read_counts_file): #
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=window_file,b=b,c=[6,7],o="sum",g="hg38.10x.nochr.fa.fai") ## this can fail?? somehow dos2unix helped?
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "name", "score", "strand_of_window",
										"hap1_reads", "hap2_reads"],dtype={"chrom":str, "start":int, "stop":int,
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
	parser.add_argument("--df",
		type=str,
		metavar="[df of allele counts]",
		required=True,
		help="in bed format")
	parser.add_argument("--window_file",
		type=str,
		metavar="[bed file of windows]",
		required=True,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--out_directory",
		type=str,
		metavar="[out directory]",
		required=True,
		help="full path to output results")

	arguments = parser.parse_args()

	df = get_windows(window_file = arguments.window_file, read_counts_file=arguments.df)

	# for binomial pval analysis. dont need plot that shows significant p values, just do the same grid plot
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],row["hap1_reads"]+row["hap2_reads"],
		p=0.5,
		alternative="two-sided"), # v slow for some reason 
		axis=1)
	df["fdr_pval"]=mt.multipletests(pvals=df["binom_pval"], alpha=0.01,
													method="fdr_bh")[1]
	df["fdr_reject"] =  mt.multipletests(pvals=df["binom_pval"], alpha=0.01,
													method="fdr_bh")[0]
	f, ax = plt.subplots(1, len(chromosomes), sharex=False,
							sharey=False,
							figsize=(15,1),
							gridspec_kw={'width_ratios': ratios})
	def color_maker(x):
		result = [0,0,0,0]
		if x["strand_of_window"]=="-":
			result[0]=1 # minus strand red
		if x["strand_of_window"]=="+":
			result[1]=1 # plus strand green
		if x["strand_of_window"]==".":
			result[2]=1 ## no strand info is blue
		if x["fdr_reject"]<=0.05:
			result[3]=1
		if x["fdr_reject"]>0.05:
			result[3]=0.1
		return result
	def marker_size(x):
		if x["total_reads"] <=100:
			return 10
		if 100 < x["total_reads"] <= 1000:
			return 40
		if 1000 < x["total_reads"]:
			return 75
	print(df)
	for i in range(len(chromosomes)):
		if chromosomes[i] not in list(df.chrom):
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
		hap1 = df[(df["hap1_reads"] >= df["hap2_reads"]) & (df["chrom"]==chromosomes[i]) ]
		if len(hap1)>0:
			hap1["total_reads"] = hap1["hap1_reads"] + hap1["hap2_reads"]
			hap1["color"] = hap1.apply(color_maker,axis=1)
			hap1["size"] = hap1.apply(marker_size,axis= 1 )
			ax[i].scatter(hap1["start"]/10**6,
					[-np.log10(x) if ( x >= 10**-12 ) else 12 for x in hap1["binom_pval"]],
			 		s=hap1["size"], # this scaling sucks still
					lw=0.1,
					edgecolor="black",
					c=hap1["color"])
			for index,row in hap1[hap1["fdr_reject"]==True].iterrows():
				ax[i].annotate(s=row["name"],
				xy=(row["start"]/10**6, -np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12)).set_fontsize(6)
		### this breaks if the df is empty...thecrappy fix is to assign windows where hap1 reads = hap2 reads to both hap1 and hap2.
		### it may double plot them, but they will be non significant anyways so who cares?
		hap2 = df[(df["hap1_reads"] <= df["hap2_reads"]) & (df["chrom"]==chromosomes[i]) ]
		if len(hap2)>0:
			hap2["total_reads"] = hap2["hap1_reads"] + hap2["hap2_reads"]
			hap2["color"] = hap2.apply(color_maker,axis=1)
			hap2["size"] = hap2.apply(marker_size,axis=1)
			ax[i].scatter(hap2["start"]/10**6,
					[np.log10(x) if ( x >= 10**-12 ) else -12 for x in hap2["binom_pval"]],
			 		s=hap2["size"],
					lw=0.1,
					edgecolor="black",
					c=hap2["color"])
			for index,row in hap2[hap2["fdr_reject"]==True].iterrows():
				ax[i].annotate(s=row["name"],
				xy=(row["start"]/10**6, -1*(-np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12))).set_fontsize(6)

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
	plt.savefig(os.path.basename(arguments.df.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".png",dpi=400,transparent=True,bbox_inches='tight',pad_inches = 0)

	### output files with proper filenames
	df.to_csv(os.path.basename(arguments.df.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".bed",sep="\t",index=None,header=None)
	browser_df = df[df["fdr_reject"]==True].loc[:,:"strand_of_window"]
	browser_df["chrom"] = "chr"+browser_df["chrom"].astype(str)
	browser_df.to_csv(os.path.basename(arguments.df.rstrip(".bed"))+"."+os.path.basename(arguments.window_file.rstrip(".bed"))+".browser.bed",sep="\t",index=None,header=None)

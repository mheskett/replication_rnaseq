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

def get_windows(window_file, read_counts_file): #
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest

	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=window_file,b=b,c=[6,7],o="sum",g="hg38.10x.nochr.sort.fa.fai") ## this can fail?? somehow dos2unix helped?
	df =  window_read_counts.to_dataframe(names=["chrom","start","stop","name","score","strand_of_window",
										"hap1_reads","hap2_reads"])
	df = df[(df["hap1_reads"]!=".") & (df["hap2_reads"]!=".")]
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
	f, ax = plt.subplots(4, 6, sharex=False,
							sharey=False,
							figsize=(22,8))
	positions = ((0,0),(0,1),(0,2),(0,3),(0,4),(0,5),
		(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),
		(2,0),(2,1),(2,2),(2,3),(2,4),(2,5),
		(3,0),(3,1),(3,2),(3,3),(3,4))
	def color_maker(x):
		result = [0,0,0,0]
		if x["strand_of_window"]=="-":
			result[0]=1
		if x["strand_of_window"]=="+":
			result[1]=1
		if x["binom_pval"]<=0.0001:
			result[3]=1
		if x["binom_pval"]>0.0001:
			result[3]=0.2
		return result
	def marker_size(x):
		if x["total_reads"] <=100:
			return 10
		if 100 < x["total_reads"] <= 1000:
			return 75
		if 1000 < x["total_reads"]:
			return 140
	for i in range(len(chromosomes)):
		if chromosomes[i] not in list(df.chrom):
			continue
		## clean up this bullcrap
		#### make a df for each

		hap1 = df[(df["hap1_reads"] > df["hap2_reads"]) & (df["chrom"]==chromosomes[i]) ]
		hap1["total_reads"] = hap1["hap1_reads"] + hap1["hap2_reads"]
		hap1["color"] = hap1.apply(color_maker,axis=1)
		hap1["size"] = hap1.apply(marker_size,axis= 1 )
		#print(   [-np.log10(x) if x>10**-12 else 12for x in list(hap1["binom_pval"])]   )
		hap2 = df[(df["hap1_reads"] < df["hap2_reads"]) & (df["chrom"]==chromosomes[i]) ]
		hap2["total_reads"] = hap2["hap1_reads"] + hap2["hap2_reads"]
		hap2["color"] = hap2.apply(color_maker,axis=1)
		hap2["size"] = hap2.apply(marker_size,axis=1)

		####below already works
		# x1 = df[(df["hap1_reads"] > df["hap2_reads"]) & (df["chrom"]==chromosomes[i])]["start"]/10**6
		# x2 = df[(df["hap1_reads"] < df["hap2_reads"]) & (df["chrom"]==chromosomes[i])]["start"]/10**6
		# y1 = -np.log10(df[(df["hap1_reads"] > df["hap2_reads"]) & (df["chrom"]==chromosomes[i])]["binom_pval"])
		# y2 = -np.log10(df[(df["hap1_reads"] < df["hap2_reads"]) & (df["chrom"]==chromosomes[i])]["binom_pval"])*-1
		# y1 = [x if x < 12 else 12 for x in y1] # so we can see outliers
		# y2 = [x if x > -12 else -12 for x in y2]
		# colors1 = [[1,0,0,1] if x >=4 else [1,1,1,0.1] for x in y1] ## pval for significance here set to 10e-4
		# colors2 = [[1,0,0,1] if x <= -4 else [1,1,1,0.1] for x in y2]
		######## 
		### repliseq data
		# if arguments.repli_seq:
		# 	repliseq = df_repliseq[df_repliseq["chrom"]==chromosomes[i]]
		# 	for index,row in repliseq.iterrows():
		# 		ax[positions[i][0],positions[i][1]].axvspan(xmin = row["start"]/10**6, xmax = row["stop"]/10**6,
		# 													ymin=0, ymax=1, alpha=0.5, facecolor="blue")
		## annotate
		for index,row in hap1[hap1["binom_pval"]<=0.0001].iterrows():
			ax[positions[i][0],positions[i][1]].annotate(s=row["name"],
				xy=(row["start"]/10**6, -np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12)).set_fontsize(6)

		for index,row in hap2[hap2["binom_pval"]<=0.0001].iterrows():
			ax[positions[i][0],positions[i][1]].annotate(s=row["name"],
			xy=(row["start"]/10**6, -1*(-np.log10(row["binom_pval"]) if -np.log10(row["binom_pval"]) <= 12 else 12))).set_fontsize(6)

		## plot	## need to address values that go above axes
		ax[positions[i][0],positions[i][1]].scatter(hap1["start"]/10**6,
					-np.log10(hap1["binom_pval"]),
			 		s=hap1["size"], # this scaling sucks still
					lw=0.1,
					edgecolor="black",
					c=hap1["color"])
		ax[positions[i][0],positions[i][1]].scatter(hap2["start"]/10**6,
					np.log10(hap2["binom_pval"]) ,
			 		s=hap2["size"],
					lw=0.1,
					edgecolor="black",
					c=hap2["color"])
		# ax[positions[i][0],positions[i][1]].scatter(hap1[hap1["chrom"]==chromosomes[i]]["start"],y1,
		# 	 		s=hap1["total_reads"]/10,
		# 			lw=0.1,
		# 			edgecolor="black",
		# 			c=hap1["strand_color"])
		# ax[positions[i][0],positions[i][1]].scatter(x2,y2,
		# 	 		s=hap2["total_reads"]/10,
		# 			lw=0.1,
		# 			edgecolor="black",
		# 			c=colors2)
		ax[positions[i][0],positions[i][1]].set(xlabel=chromosomes[i]) # x axis labels or no
		#ax[positions[i][0],positions[i][1]].set_xticks([])
		ax[positions[i][0],positions[i][1]].set_xlim([min(df[df["chrom"]==chromosomes[i]]["start"].values)/10**6,
				max(df[df["chrom"]==chromosomes[i]]["stop"])/10**6])
		ax[positions[i][0],positions[i][1]].set_yticks([-12,-8,-4,0,4,8,12])
		plt.grid(True)
		ax[positions[i][0],positions[i][1]].xaxis.set_minor_locator(AutoMinorLocator(4))
		ax[positions[i][0],positions[i][1]].set_ylim([-14,14]) # for all plots

		### try to annotate all significant links
		#params = {'axes.labelsize': 18,'axes.titlesize':20, 'text.fontsize': 20, 
				#'legend.fontsize': 20, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
		#matplotlib.rcParams.update(params)
	plt.subplots_adjust(wspace=.2, hspace=.4)

	plt.savefig(arguments.df.rstrip(".bed")+"."+arguments.window_file.rstrip(".bed")+".png",dpi=400,transparent=True,bbox_inches='tight',pad_inches = 0)

	df.to_csv(arguments.df.rstrip(".bed")+"."+arguments.window_file.rstrip(".bed")+".txt",sep="\t",index=None,header=None)

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

### Separate script for window version since this is a totally different procedure and windows are made on command line?
### Or, use pybedtools to make windows within this script, but then must provide the reference genome location

def get_windows(length,bed_file): #ezpz
	a=pybedtools.BedTool()

	b = pybedtools.BedTool(bed_file)
	windows=a.window_maker(g="/Users/mike/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
						w=length,s=length/2)

	c = pybedtools.BedTool()
	window_read_counts = c.map(a=windows,b=b,c=[6,7],o="sum")
	df =  window_read_counts.to_dataframe(names=["chrom","start","stop",
										"hap1_reads","hap2_reads"])

	df = df[(df["hap1_reads"]!=".") & (df["hap2_reads"]!=".")]
	return df

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	parser.add_argument("--df",
		type=str,
		metavar="[df of allele counts]",
		required=True,
		help="in bed format")
	parser.add_argument("--out_directory",
		type=str,
		metavar="[out directory]",
		required=True,
		help="full path to output results")
	parser.add_argument("--window_size",
		type=int,
		metavar="[window_size]",
		required=False,
		default=100000,
		help="window size for window analysis")
	parser.add_argument("--binomial_test",
		action="store_true",
		metavar="[binom test?]",
		required=False,
		help="use binomial test instead of log ratio")
	parser.add_argument("--windows",
		metavar="[windows?]",
		action="store_true",
		required=False,
		default=100000,
		help="use windows instead of single snps")

	arguments = parser.parse_args()
	if not arguments.windows:
		df = pd.read_table(arguments.df, # for original data
							header=None,
							names=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"],
							dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"hap1_reads":int,"hap2_reads":int})
		df = df[df["hap1_reads"]+df["hap2_reads"] >= 10] # min 15 reads per site

	if arguments.windows:
		df = get_windows(length=arguments.window_size,bed_file=agruments.df)
		df = df[df["hap1_reads"]+df["hap2_reads"] >= 30] # for windows

	if not arguments.binomial_test:
	    ## happens for non-binomial test version
		df["hap1_logR"] = np.log2((df["hap1_reads"] / (df["hap1_reads"]+df["hap2_reads"])) / 0.5 )
		df["hap2_logR"] = np.log2((df["hap2_reads"] / (df["hap1_reads"]+df["hap2_reads"])) / 0.5 )
		f, ax = plt.subplots(1,len(chromosomes),sharex=False,
												sharey=False,
												figsize=(14,1))
		for i in range(len(chromosomes)):

			x1 = df[(df["hap1_logR"] > 0) & (df["chrom"]==chromosomes[i])]["start"] # this is for read ratios
			x2 = df[(df["hap2_logR"] > 0) & (df["chrom"]==chromosomes[i])]["start"]
			y1 = df[(df["hap1_logR"] > 0) & (df["chrom"]==chromosomes[i])]["hap1_logR"] 
			y2 = -df[(df["hap2_logR"] > 0) & (df["chrom"]==chromosomes[i])]["hap2_logR"]
			colors1 = [[1,0,0,x**3] for x in y1] # now you have linear alpha
			colors2 = [[0,1,0,abs(x)**3] for x in y2]
			ax[i].scatter(x1,y1,
				 		s=4,
						lw=0.1,
						edgecolor="black",
						c=colors1)
			ax[i].scatter(x2,y2,
				 		s=4,
						lw=0.1,
						edgecolor="black",
						c=colors2)
			ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
			ax[i].set_xticks([])
			ax[i].set_xlim([min(df[df["chrom"]==chromosomes[i]]["start"].values),
					max(df[df["chrom"]==chromosomes[i]]["stop"])])
			if chromosomes[i]=="1":
				ax[i].set_yticks([-2,-1,0,1,2])
			else:
				ax[i].set_yticks([])
			plt.grid(True)
			ax[i].axhline(y=2,linestyle="--",color="red")
			ax[i].set_ylim([-1.05,1.05]) # for all plots
		#ax[i].axvline(x=int(centromere[chromosomes[i]]),linestyle = "--", lw = 0.5,color="black")
		f.subplots_adjust(wspace=0, hspace=0)
		#plt.show()


	if arguments.binomial_test:
    # for binomial pval analysis
		df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],row["hap1_reads"]+row["hap2_reads"],
			p=0.5,
			alternative="two-sided"), # v slow for some reason 
			axis=1)
		f, ax = plt.subplots(1,len(chromosomes),sharex=False,
												sharey=False,
												figsize=(14,1))
		for i in range(len(chromosomes)):

			x1 = df[(df["chrom"]==chromosomes[i])]["start"] # this is for read ratios
			y1 = -np.log10(df[(df["chrom"]==chromosomes[i])]["binom_pval"])
			colors = ['red' if val >= 2 else 'blue' for val in y ] for log p val
			ax[i].scatter(x1,y1,
				 		s=4,
						lw=0.1,
						edgecolor="black",
						alpha=0.6,
						c=colors)
			ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
			ax[i].set_xticks([])
			ax[i].set_xlim([min(df[df["chrom"]==chromosomes[i]]["start"].values),
					max(df[df["chrom"]==chromosomes[i]]["stop"])])
			if chromosomes[i]=="1":
				ax[i].set_yticks([0,1,2,3,4,5,6,7,8])
			else:
				ax[i].set_yticks([])
			plt.grid(True)
			ax[i].axhline(y=2,linestyle="--",color="red")
			ax[i].set_ylim([0,8]) # for all plots
			ax[i].axvline(x=int(centromere[chromosomes[i]]),linestyle = "--", lw = 0.5,color="black")
		f.subplots_adjust(wspace=0, hspace=0)

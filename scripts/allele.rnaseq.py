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

def get_windows(length): #ezpz
	a=pybedtools.BedTool()

	b = pybedtools.BedTool(arguments.df)
	windows=a.window_maker(g="/Users/mike/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
						w=length,s=length/2)
	
	c = pybedtools.BedTool()
	window_read_counts = c.map(a=windows,b=b,c=[6,7],o="sum")
	df =  window_read_counts.to_dataframe(names=["chrom","start","stop",
										"hap1_reads","hap2_reads"])

	df = df[(df["hap1_reads"]!=".") & (df["hap2_reads"]!=".")]
	return df

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	
	parser.add_argument("--df",
		type=str,
		metavar="[df of allele counts]",
		required=True,
		help="")
	
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
		help="full path to output results")
	
	arguments = parser.parse_args()
	print(get_windows(length=arguments.window_size))
	chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
	"13","14","15","16","17","18","19","20","21","22","X"]

	df = pd.read_table(arguments.df, # for original data
						header=None,
						names=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"],
						dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"hap1_reads":int,"hap2_reads":int})

	# df = pd.read_table(arguments.df, # for window
	# 					header=None,
	# 					names=["chrom","start","stop","hap1_reads","hap2_reads"],
	# 					dtype={"chrom":str,"start":int,"stop":int,"hap1_reads":int,"hap2_reads":int})

	# df["hap1_reads_log"] = np.log2(df["hap1_reads"]+1)
	# df["hap2_reads_log"] = np.log2(df["hap2_reads"]+1)
    

    ## for binomial pval analysis
	# df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],row["hap1_reads"]+row["hap2_reads"],
	# 	p=0.5,
	# 	alternative="two-sided"), # v slow for some reason 
	# 	axis=1)
	df = df[df["hap1_reads"]+df["hap2_reads"] >= 15] # not required for window version, this is filtering
	df["hap1_logR"] = np.log2((df["hap1_reads"] / (df["hap1_reads"]+df["hap2_reads"])) / 0.5 )
	df["hap2_logR"] = np.log2((df["hap2_reads"] / (df["hap1_reads"]+df["hap2_reads"])) / 0.5 )
	#print(df)
	f, ax = plt.subplots(1,len(chromosomes),sharex=False,
											sharey=False,
											figsize=(14,1))
	
	for i in range(len(chromosomes)):
		#x = df[df["chrom"]==chromosomes[i]]["start"] for doing binom pval
		#y = -np.log10(df[df["chrom"]==chromosomes[i]]["binom_pval"])
		x1 = df[(df["hap1_logR"] > 0) & (df["chrom"]==chromosomes[i])]["start"] # this is for read ratios
		x2 = df[(df["hap2_logR"] > 0) & (df["chrom"]==chromosomes[i])]["start"]

		y1 = df[(df["hap1_logR"] > 0) & (df["chrom"]==chromosomes[i])]["hap1_logR"] 
		y2 = -df[(df["hap2_logR"] > 0) & (df["chrom"]==chromosomes[i])]["hap2_logR"]
		#colors = ['red' if val >= 2 else 'blue' for val in y ] for log p val
		colors1 = [[1,0,0,x**3] for x in y1] # now you have linear alpha
		colors2 = [[0,1,0,abs(x)**3] for x in y2]
		ax[i].scatter(x1,y1,# range on x.size here for ultra smooth TM
			 		s=4,
					lw=0.1,
					edgecolor="black",
					c=colors1)
		ax[i].scatter(x2,y2,# range on x.size here for ultra smooth TM
			 		s=4,
					lw=0.1,
					edgecolor="black",
					c=colors2)
		ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
		ax[i].set_xticks([])
		#ax[i].set_ylim([0,8])

		ax[i].set_xlim([min(df[df["chrom"]==chromosomes[i]]["start"].values),
				max(df[df["chrom"]==chromosomes[i]]["stop"])])
		if chromosomes[i]=="1":
			#ax[i].set_yticks([0,1,2,3,4,5,6,7,8])
			ax[i].set_yticks([-2,-1,0,1,2])

		else:
			ax[i].set_yticks([])

		#ax[i].grid(True)
		plt.grid(True)
		ax[i].axhline(y=2,linestyle="--",color="red")
		ax[i].set_ylim([-1.05,1.05]) # for all plots

	#ax[i].axvline(x=int(centromere[chromosomes[i]]),linestyle = "--", lw = 0.5,color="black")

	f.subplots_adjust(wspace=0, hspace=0)
	plt.show()

    #### create pybedtools object with windows
    #### use bedtools map to summarize read counts
    #### or do this on command line first

    #### bedtools map -a hg38.nochr.25kb.windows.bed -b gm12878.rep1Aligned.hcgatk3.overlap.platinum.haplotype.resolved.bed -c 6,7 -o sum | grep -Fv "." |

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
	
	arguments = parser.parse_args()

	chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
	"13","14","15","16","17","18","19","20","21","22","X"]

	df = pd.read_table(arguments.df,
						header=None,
						names=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"],
						dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"hap1_reads":int,"hap2_reads":int})

	df["hap1_reads_log"] = np.log2(df["hap1_reads"]+1)
	df["hap2_reads_log"] = np.log2(df["hap2_reads"]+1)

	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],row["hap1_reads"]+row["hap2_reads"],
		p=0.5,alternative="two-sided"),
		axis=1)
	print(df)
	# plt.scatter(df[df["chrom"]=="1"]["start"],df[df["chrom"]=="1"]["hap1_reads"],
	# 	s=3,lw="0.2",edgecolor="black",c="blue")
	# plt.scatter(df[df["chrom"]=="1"]["start"],-(df[df["chrom"]=="1"]["hap2_reads"]),
	# 	s=3,lw="0.2",edgecolor="black",c="red")

	f, ax = plt.subplots(1,len(chromosomes),sharex=False,sharey=False,figsize=(14,1))
	
	for i in range(len(chromosomes)):
		x = df[df["chrom"]==chromosomes[i]]["start"]
		y = -np.log10(df[df["chrom"]==chromosomes[i]]["binom_pval"])
		ax[i].scatter(x,y,# range on x.size here for ultra smooth TM
			 		s=2,
					lw=0.2,
					edgecolor="black",
					alpha=0.6)
		ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
		ax[i].set_xticks([])
		ax[i].set_ylim([0,6])

		if chromosomes[i]=="1":
			ax[i].set_yticks([0,1,2,3,4,5,6])

		else:
			ax[i].set_yticks([])
		ax[i].set_xlim([min(df[df["chrom"]==chromosomes[i]]["start"].values),
						max(df[df["chrom"]==chromosomes[i]]["stop"])])
		#ax[i].grid(True)
		plt.grid(True)
		ax[i].axhline(y=2,linestyle="--",color="red")
	#ax[i].axvline(x=int(centromere[chromosomes[i]]),linestyle = "--", lw = 0.5,color="black")

	f.subplots_adjust(wspace=0, hspace=0)
	plt.show()
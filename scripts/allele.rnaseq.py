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

	df = pd.read_table(arguments.df,
						header=None,
						names=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"],
						dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"hap1_reads":int,"hap2_reads":int})
	df["hap1_reads_log"] = np.log2(df["hap1_reads"]+1)
	df["hap2_reads_log"] = np.log2(df["hap2_reads"]+1)

	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],row["hap1_reads"]+row["hap2_reads"],p=0.5,alternative="two-sided"),
		axis=1)
	print(df)
	# plt.scatter(df[df["chrom"]=="1"]["start"],df[df["chrom"]=="1"]["hap1_reads"],
	# 	s=3,lw="0.2",edgecolor="black",c="blue")
	# plt.scatter(df[df["chrom"]=="1"]["start"],-(df[df["chrom"]=="1"]["hap2_reads"]),
	# 	s=3,lw="0.2",edgecolor="black",c="red")

	plt.scatter(df[df["chrom"]=="1"]["start"],-np.log10(df[df["chrom"]=="1"]["binom_pval"]))
	plt.xlim([0,max(df[df["chrom"]=="1"]["stop"].values)])
	plt.ylim([-9,9])
	plt.ylim([-100,100])
	plt.ylim([0,6])
	plt.grid(True)


	plt.show()
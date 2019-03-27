import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
import seaborn as sns
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
						names=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"])
	df["hap1_reads_log"] = np.log2(df["hap1_reads"]+1)
	df["hap2_reads_log"] = np.log2(df["hap2_reads"]+1)
	print(df)

	plt.scatter(df[df["chrom"]==1]["start"],df[df["chrom"]==1]["hap1_reads"],
		s=3,lw="0.2",edgecolor="black",c="blue")
	plt.scatter(df[df["chrom"]==1]["start"],-(df[df["chrom"]==1]["hap2_reads"]),
		s=3,lw="0.2",edgecolor="black",c="red")
	plt.xlim([0,max(df[df["chrom"]==1]["stop"])])
	plt.ylim([-9,9])
	plt.ylim([-100,100])
	plt.grid(True)


	plt.show()
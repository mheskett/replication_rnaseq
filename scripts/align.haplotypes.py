import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="takes mpileup VCF file turns into bed")

	parser.add_argument("--bed",
	type=str,
	metavar="[table file]",
	required=True,
	help="")
	parser.add_argument("--out_directory", #differentiate between out directory and snapshot directory
	type=str,
	metavar="[out directory]",
	required=True,
	help="full path to output results")
	### warning: this program expects a very specific file format
	### must have done bedtools intersect -wa -wb where a is a bed file
	### remove triallelic sites and separate the AD string into ref and alt counts
	arguments = parser.parse_args()
	result = []
	print("file",arguments.bed)
	with open(arguments.bed,'r') as f:
		for line in f:
			try:
				tmp = line.rstrip().split()
				ref = tmp[3]
				alt = tmp[4].split(",")

				if len(alt) == 1 and alt[0] == "<*>":
					alt = "<*>"
				elif len(alt) == 2:
					alt = alt[0]
				elif len(alt) == 3:
					continue ## continue to next line to filter out triallelic
				counts = list(map(int,tmp[5].split(",")))[0:2]

				hap1 = tmp[9]
				hap2 = tmp[10]

				ref_count = counts[0]
				alt_count = counts[1]

				if (ref == hap1) and (alt == hap2): ## no switching necessary
					result += [[tmp[6], tmp[7], tmp[8], hap1, hap2, ref_count, alt_count]]
				elif (ref == hap2) and (alt == hap1): ## switch
					result += [[tmp[6], tmp[7], tmp[8], hap1, hap2, alt_count, ref_count]]
				elif ref == hap1 and alt == "<*>":
					result += [[tmp[6], tmp[7], tmp[8], hap1, hap2, ref_count, alt_count]]
				elif ref == hap2 and alt == "<*>":
					result += [[tmp[6], tmp[7], tmp[8], hap1, hap2, alt_count, ref_count]]
				else:
					continue ### these are likely indels that weren't called correctly by mpileup
			except:
				print("error :",arguments.bed)
				print(line)
				continue

	df = pd.DataFrame(result,columns=["chrom","start","stop","hap1","hap2","hap1_reads","hap2_reads"])
	df.to_csv(arguments.out_directory + os.path.basename(arguments.bed.replace(".bed",'')+ ".haplotype.resolved.counts.bed"), 
		sep="\t", header=False, index=False)



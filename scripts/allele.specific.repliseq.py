import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="takes mpileup VCF file turns into bed")

	parser.add_argument("--table",
	type=str,
	metavar="[table file]",
	required=True,
	help="")
	parser.add_argument("--out_directory", #differentiate between out directory and snapshot directory
	type=str,
	metavar="[out directory]",
	required=True,
	help="full path to output results")


	arguments = parser.parse_args()

	data = []
	## just simply reformatting from table to bed. remove trialleleic
	with open(arguments.table,'r') as f:
		next(f) # skip header
		for line in f:
			tmp = line.rstrip().split()
			## gets rid of multi allelic sites
			if len(tmp[3].split(",")) > 2:
				continue
			counts = tmp[5].split(",")

			data += [[tmp[0], int(tmp[1])-1, int(tmp[1]), tmp[2], tmp[3].split(",")[0], counts[0], counts[1]]]

		### 
	### now arrange everything by haplotype 
	df = pd.DataFrame(data,columns=["chrom","start","stop","ref","alt","ref_reads","alt_reads"])
	df.to_csv(arguments.out_directory + os.path.basename(arguments.table.replace(".table",'')+".bed"),header=False, index=False, sep="\t")

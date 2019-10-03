import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
import multiprocessing as mp

def sample_introns(length):
	# creates a df of introns where the sum of lengths equals length
	tmp = df_introns.sample(n=1)
	while tmp["intron_length"].sum() < length:
		tmp = pd.concat([tmp,df_introns.sample(n=1)], ignore_index=True)
	remainder = tmp["intron_length"].sum() - length

	if remainder > 0:
		# print(tmp)
		###
		last_row = tmp.tail(1).reset_index(drop=True)
		remainder_start = last_row.at[0, "start"]
		remainder_stop = last_row.at[0, "stop"]
		remainder_chrom = last_row.at[0, "chrom"]
		remainder_intron_name = last_row.at[0, "intron_name"]
		remainder_score = 0
		remainder_strand = last_row.at[0, "strand"]
		new_remainder_stop = remainder_stop - remainder

		bedtool_string = str(remainder_chrom)+' '+str(remainder_start)+' '+str(new_remainder_stop) + ' ' + str(remainder_intron_name) + ' ' + str(remainder_score)+ ' ' + str(remainder_strand) 
		# convert new intron to bed tool, get line counts
		try:
			new_intron_fragment = pybedtools.BedTool(bedtool_string, from_string=True)
			# print(bedtool_string)
		except:
			print(bedtool_string)
			print(tmp)
			print("ERROR...handling error")
			return sample_introns(length) # this could go infite loop if errors keep happening, but since they are random i dont think it can
		new_intron_fragment_lines = new_intron_fragment.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","intron_name","score",
																									"strand","line_count","num_bases_covered",
																									"intron_length","fraction_l1"],
																							dtype={"chrom":str,"start":int,"stop":int,"intron_name":str,
																							"score":int,"strand":str,"line_count":int,
																							"num_bases_covered":int,"intron_length":int,"fraction_l1":float})
		tmp = tmp.drop(tmp.tail(1).index)
		tmp = pd.concat([tmp,new_intron_fragment_lines], ignore_index=True)

	sum_intron_length = tmp["intron_length"].sum()
	sum_covered_bases = tmp["num_bases_covered"].sum()
	intron_fraction_l1 = sum_covered_bases / sum_intron_length

	if sum_intron_length != length:
		print("DEBUG ALERT")
		print(tmp)
		print("sum intron length ", tmp["intron_length"].sum())
		print("length: ", length)
		print("remainder: ", remainder)
		print(bedtool_string)
	if intron_fraction_l1 == 0:
		return (sum_covered_bases + 1) / sum_intron_length

	return intron_fraction_l1

def get_line_fraction(chrom, start, stop):
	my_bedtool = pybedtools.BedTool(str(chrom)+' '+str(start)+' '+str(stop), from_string=True)
	tmp = my_bedtool.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","line_count","num_bases","input_length","fraction_l1"],
														dtype={"chrom":str,"start":int,"stop":int,"line_count":int,"num_bases":int,"input_length":int,"fraction_l1":float})
	return tmp["fraction_l1"].values[0]

def get_pvalue_line_content(vlinc_name, vlinc_length, line_fraction, num_simulations):
	results = []
	for i in range(0, num_simulations):
		results += [sample_introns(length = vlinc_length)]
	### fit beta to results
	a1,b1,loc1,scale1 = scipy.stats.beta.fit(results, floc=0, fscale=1)
	x = np.linspace(0,1,1000)
	y = scipy.stats.beta.pdf(x, a1, b1)
	beta_pvalue = 1 - scipy.stats.beta.cdf(line_fraction, a1, b1)
	plt.plot(x, y)
	plt.show()
	plt.hist(results,bins=20)
	plt.show()
	pvalue = (100 - scipy.stats.percentileofscore(results, line_fraction, kind='strict')) / 100
	if pvalue == 0:
		pvalue = 1/num_simulations

	### or center and scale data then fit a normal, or fit a poisson or whatever
	return vlinc_name, line_fraction, pvalue, beta_pvalue


introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs_file = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

introns = pybedtools.BedTool(introns_file)
vlincs = pybedtools.BedTool(vlincs_file)
lines = pybedtools.BedTool(lines_file)

df_introns = introns.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","intron_name","score","strand","line_count","num_bases_covered","intron_length","fraction_l1"],
														dtype={"chrom":str,"start":int,"stop":int,"intron_name":str,"strand":str,"line_count":int,"num_bases":int,"intron_length":int,"fraction_l1":float}) ## the line should be at least 90 % covered by the intron 
# df_introns.to_csv("test_introns.txt",sep="\t",index=None)
# test = pd.read_table("test_introns.txt",sep="\t",index_col=None)
# print(test)
df_vlincs = vlincs.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","vlinc_name","score","strand","line_count","num_bases_covered","vlinc_length","fraction_l1"],
														dtype={"chrom":str,"start":int,"stop":int,"vlinc_name":str,"strand":str,"line_count":int,"num_bases":int,"vlinc_length":int,"fraction_l1":float}) ## the line should be at least 90 % covered by the intron )

vlinc273_row = df_vlincs[df_vlincs["vlinc_name"]=="273,313,6.141125478.141219540"]

vlinc273_line_content = get_line_fraction(chrom=vlinc273_row["chrom"].values[0], start=vlinc273_row["start"].values[0], stop=vlinc273_row["stop"].values[0])
print(get_pvalue_line_content(vlinc_name = "asar6-141", vlinc_length=vlinc273_row["vlinc_length"].values[0], line_fraction=0.28, num_simulations=100))

asar6_line_fraction = get_line_fraction(6,96075000,96279595)
print(get_pvalue_line_content(vlinc_name = "asar6-96",vlinc_length=96279595 - 96075000 , line_fraction=0.38, num_simulations=100))


### 
# results_regression = []
# for i in range(50000,500000,10000):
# 	results_regression += [sample_introns(i)]


### main program
# pool = mp.Pool()
# results = pool.starmap(get_pvalue_line_content,[(row.vlinc_name, row.vlinc_length, row.fraction_l1, 1000) for row in df_vlincs.itertuples()]) ##  each link needs to be list of list in parallel call

# with open("vlinc_l1_significance.txt", "w") as f:
# 	writer = csv.writer(f,delimiter="\t") 
# 	for i in range(len(results)):
# 		writer.writerow(results[i])
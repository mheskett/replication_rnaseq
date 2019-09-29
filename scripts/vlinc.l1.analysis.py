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

def sample_introns(length):
	# creates a df of introns where the sum of lengths equals length
	tmp = df_introns.sample(n=1)
	while tmp["intron_length"].sum() < length:
		tmp = tmp.append(df_introns.sample(n=1))
	remainder = tmp["intron_length"].sum() - length
	if remainder > 0:
		# print(tmp)
		###
		remainder_start = tmp.at[tmp.index[-1], "start"]
		remainder_stop = tmp.at[tmp.index[-1], "stop"]
		remainder_chrom = tmp.at[tmp.index[-1], "chrom"]
		new_remainder_stop = remainder_stop - remainder
		# convert new intron to bed tool, get line counts
		new_intron_fragment = pybedtools.BedTool(remainder_chrom+' '+str(remainder_start)+' '+str(new_remainder_stop), from_string=True)
		new_intron_fragment_lines = new_intron_fragment.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","line_count","num_bases_covered","intron_length","fraction_l1"],
																				dtype={"chrom":str,"start":int,"stop":int,"line_count":int,"num_bases_covered":int,"intron_length":int,"fraction_l1":float})
		# can optimize this thing by adding name score strand to the bedtools string above, then just doing drop followed by append below...
		# new_intron_fragment_line_count = new_intron_fragment_lines.at[0, "line_count"]
		# # replace the values in the tmp data frame to reflect the shortened intron while keeping the name the same
		tmp.loc[tmp.index[-1],"stop"] = new_remainder_stop # dis works
		tmp.loc[tmp.index[-1],"line_count"] = new_intron_fragment_lines.at[0,"line_count"]
		tmp.loc[tmp.index[-1],"num_bases_covered"] = new_intron_fragment_lines.at[0,"num_bases_covered"]
		tmp.loc[tmp.index[-1],"intron_length"] = new_remainder_stop - remainder_start
		tmp.loc[tmp.index[-1],"fraction_l1"] = new_intron_fragment_lines.at[0,"fraction_l1"]
		if tmp["intron_length"].sum() != length:
			print("DEBUG ALERT")
	return tmp["num_bases_covered"].sum() / tmp["intron_length"].sum()

def get_line_fraction(chrom, start, stop):
	my_bedtool = pybedtools.BedTool(str(chrom)+' '+str(start)+' '+str(stop), from_string=True)
	tmp = my_bedtool.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","line_count","num_bases","input_length","fraction_l1"],
														dtype={"chrom":str,"start":int,"stop":int,"line_count":int,"num_bases":int,"input_length":int,"fraction_l1":float})
	return tmp["fraction_l1"].values[0]

def get_pvalue_line_content(vlinc_length, line_fraction, num_simulations):
	results = []
	for i in range(0,num_simulations):
		results += [sample_introns(length = vlinc_length)]
	# plt.hist(results,bins=20)
	# plt.show()
	pvalue = (100 - scipy.stats.percentileofscore(results, line_fraction, kind='strict')) / 100
	if pvalue == 0:
		pvalue = 1/num_simulations

	### or center and scale data then fit a normal, or fit a poisson or whatever
	return pvalue


introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs_file = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

introns = pybedtools.BedTool(introns_file)
vlincs = pybedtools.BedTool(vlincs_file)
lines = pybedtools.BedTool(lines_file)

df_introns = introns.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","intron_name","score","strand","line_count","num_bases_covered","intron_length","fraction_l1"],
																	dtype={"chrom":str,"start":int,"stop":int,"intron_name":str,"strand":str,"line_count":int,"num_bases":int,"intron_length":int,"fraction_l1":float}) ## the line should be at least 90 % covered by the intron 

df_vlincs = vlincs.coverage(lines, F=0.9).to_dataframe(names=["chrom","start","stop","vlinc_name","score","strand","line_count","num_bases_covered","vlinc_length","fraction_l1"],
																	dtype={"chrom":str,"start":int,"stop":int,"vlinc_name":str,"strand":str,"line_count":int,"num_bases":int,"vlinc_length":int,"fraction_l1":float}) ## the line should be at least 90 % covered by the intron )

vlinc273_row = df_vlincs[df_vlincs["vlinc_name"]=="273,313,6.141125478.141219540"]
# print(vlinc273_row["start"].values[0])

vlinc273_line_content = get_line_fraction(chrom=vlinc273_row["chrom"].values[0], start=vlinc273_row["start"].values[0], stop=vlinc273_row["stop"].values[0])
print(get_pvalue_line_content(vlinc_length=vlinc273_row["vlinc_length"].values[0], line_fraction=vlinc273_row["fraction_l1"].values[0], num_simulations=1000))

asar6_line_fraction = get_line_fraction(6,96075000,96279595)
print(get_pvalue_line_content(vlinc_length=96279595 - 96075000 , line_fraction=asar6_line_fraction, num_simulations=1000))



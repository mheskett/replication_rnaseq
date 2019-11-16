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
		new_intron_fragment_eclip = new_intron_fragment.coverage(eclip, F=0.9).to_dataframe(names=["chrom", "start", "stop", "intron_name","score",
																									"strand","eclip_count", "num_bases_covered",
																									"intron_length", "fraction_bound"],
																							dtype={"chrom":str,"start":int,"stop":int,"intron_name":str,
																							"score":int,"strand":str,"eclip_count":int,
																							"num_bases_covered":int,"intron_length":int,"fraction_bound":float})
		tmp = tmp.drop(tmp.tail(1).index)
		tmp = pd.concat([tmp,new_intron_fragment_eclip], ignore_index=True)

	sum_intron_length = tmp["intron_length"].sum()
	sum_eclip_counts = tmp["eclip_count"].sum()
	peaks_per_length = sum_eclip_counts / sum_intron_length

	if sum_intron_length != length:
		print("DEBUG ALERT")
		print(tmp)
		print("sum intron length ", tmp["intron_length"].sum())
		print("length: ", length)
		print("remainder: ", remainder)
		print(bedtool_string)
	if peaks_per_length == 0:
		return 1 / sum_intron_length

	return peaks_per_length

def get_pvalue_eclip_peaks(vlinc_chrom, vlinc_start, vlinc_stop, vlinc_name, vlinc_length, peaks_per_length, num_simulations):
	results = []
	for i in range(0, num_simulations):
		results += [sample_introns(length = vlinc_length)]
	### fit beta to results
	a1,b1,loc1,scale1 = scipy.stats.beta.fit(results, floc=0, fscale=1)
	x = np.linspace(0,1,1000)
	y = scipy.stats.beta.pdf(x, a1, b1)
	beta_pvalue = 1 - scipy.stats.beta.cdf(peaks_per_length, a1, b1)
	# plt.plot(x, y)
	# plt.show()
	# plt.hist(results,bins=20)
	# plt.show()
	# pvalue = (100 - scipy.stats.percentileofscore(results, line_fraction, kind='strict')) / 100
	# if pvalue == 0: # beta dist is better than this i think
	# 	pvalue = 1/num_simulations
	print(vlinc_name, "done")

	### or center and scale data then fit a normal, or fit a poisson or whatever
	return vlinc_chrom, vlinc_start, vlinc_stop, vlinc_name, line_fraction, beta_pvalue


## the main question is:
## Are any particular RNBPs enriched at ASAR L1s compared to either intron L1s or whatever
## count total number of binding events in all ASAR L1s, then compare that to the number of
## binding events for an equivolent amount of sequence in intron L1s, using the beta method from
## other script

## for local use
eclip_file = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.hepg2.nochr.sorted.bed"
introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"
skewed_vlincs = "/Users/heskett/replication_rnaseq/scripts/gm12878.rep1.hg19Aligned.outgm12878.rep1.expressed.vlincs.skewed.bed"
introns = pybedtools.BedTool(introns_file)
lines = pybedtools.BedTool(lines_file)
vlincs = pybedtools.BedTool(skewed_vlincs)
################
### get rid of the non-significant peaks
df_eclip = pd.read_csv(eclip_file,sep="\t",
					names=["chrom", "start", "stop", "name","score", 
					"strand", "log_fold_enrichment","log_pvalue","qvalue", "peak"],
					 dtype={"chrom":str, "start":int, "stop":int,
					 "name":str, "score":float, "strand":str,"log_fold_enrichment":float,"log_pvalue":float,"qvalue":int,"peak":int})
df_eclip = df_eclip[ (~df_eclip.name.str.contains("IDR",case=False)) & (df_eclip["log_fold_enrichment"] >= 1) & (df_eclip["log_pvalue"] >= 2) ] # filter out replication experiments
df_eclip = df_eclip.loc[:,["chrom", "start", "stop", "name","score", "strand", "log_fold_enrichment","log_pvalue"]]
eclip_experiments = list( df_eclip.name.unique() ) ## each protein is separate
#################
df_vlincs = pd.read_csv(skewed_vlincs,sep="\t",names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
																"binom_pval","fdr_pval","fdr_reject","total_reads","skew"],
															dtype={"chrom":str,"start":int, "stop":int,"name_lncrna":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
																"binom_pval":float,"fdr_pval":float,"fdr_reject":str,"total_reads":int,"skew":float}) 
vlincs = pybedtools.BedTool.from_dataframe(df_vlincs.loc[:,["chrom","start","stop","name_lncrna","rpkm","strand","hap1_reads","hap2_reads","binom_pval","skew"]])
################
lines_in_vlincs = lines.intersect(vlincs,f=0.5,wa=True)
###############
eclip = pybedtools.BedTool.from_dataframe(df_eclip)
binding_in_vlincs = eclip.intersect(vlincs,f=0.5,wa=True,s=True,wb=True) # this binding does need to be stranded
## this is all eclip within vlincs
df_binding_in_vlincs = binding_in_vlincs.to_dataframe(names=["chrom_eclip", "start_eclip", "stop_eclip", "name_eclip","score_eclip", "strand_eclip", "log_fold_enrichment","log_pvalue",
															"chrom_vlinc","startvlinc","stop_vlinc","name_vlinc","rpkm","strand_vlinc","hap1_reads","hap2_reads",
																"binom_pval","skew"],
													dtype={"chrom_eclip":str, "start_eclip":int, "stop_eclip":int,"name_eclip":str, "score_eclip":float, "strand_eclip":str,"log_fold_enrichment":float,"log_pvalue":float,
													"chrom_vlinc":str,"start_vlinc":int, "stop_vlinc":int,"name_vlinc":str,"rpkm":float,"strand_vlinc":str,"hap1_reads":int,"hap2_reads":int,"binom_pval":float,"skew":float})
## this is all eclip within lines within vlincs
binding_in_lines_in_vlincs = lines_in_vlincs.intersect(binding_in_vlincs,wa=True,wb=True)

print(len(binding_in_lines_in_vlincs))

##############



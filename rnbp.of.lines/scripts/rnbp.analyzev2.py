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


def line_sampler(df_lines_in_introns_coverage, eclip, length):
	tmp = [df_lines_in_introns_coverage.sample(n=1)]
	tmp_length = tmp[0]["line_length"].values[0]
	while tmp_length < length: 
		sample = df_lines_in_introns_coverage.sample(n=1)
		tmp.append( sample )
		tmp_length = tmp_length + sample["line_length"].values[0]
	remainder = tmp_length - length
	tmp = pd.concat(tmp, ignore_index=True) # this fixes the quadratic copy of doing a concat in a loop
	## above can still easily be optimized because its slow when you have small size of samples but big total sequence
	if remainder > 0:
		last_row = tmp.tail(1).reset_index(drop=True)
		remainder_chrom = last_row.at[0, "chrom"]
		remainder_start = last_row.at[0, "start"]
		remainder_stop = last_row.at[0, "stop"]
		remainder_line_name = last_row.at[0, "line_name"]
		remainder_score = 0
		remainder_strand = last_row.at[0, "strand"]
		new_remainder_stop = remainder_stop - remainder
		bedtool_string = str(remainder_chrom) + ' ' + str(remainder_start) + ' ' + str(new_remainder_stop) + ' ' + str(remainder_line_name) + ' ' + str(remainder_score)+ ' ' + str(remainder_strand) 
		# convert new intron to bed tool, get line counts
		try:
			new_line_fragment = pybedtools.BedTool(bedtool_string, from_string=True)
		except:
			print(bedtool_string)
			print(tmp)
			print("ERROR...handling error")
			return sample_introns(length) # this could go infite loop if errors keep happening, but since they are random i dont think it can
		new_line_fragment_coverage = new_line_fragment.coverage(eclip, F=0.5).to_dataframe(names=["chrom", "start", "stop", "line_name","score",
																									"strand","eclip_count", "num_bases_covered",
																									"line_length", "fraction_bound"],
																							dtype={"chrom":str,"start":int,"stop":int,"line_name":str,
																							"score":int,"strand":str,"eclip_count":int,
																							"num_bases_covered":int,"line_length":int,"fraction_bound":float})
		tmp = tmp.drop(tmp.tail(1).index)
		tmp = pd.concat([tmp,new_line_fragment_coverage], ignore_index=True)

	sum_line_length = tmp["line_length"].sum()
	sum_eclip_counts = tmp["eclip_count"].sum()
	peaks_per_length = sum_eclip_counts / sum_line_length
	if sum_line_length != length:
		print("DEBUG ALERT")
		print(tmp)
		print("sum intron length ", tmp["intron_length"].sum())
		print("length: ", length)
		print("remainder: ", remainder)
		print(bedtool_string)
	if peaks_per_length == 0:
		return 1 / sum_line_length

	return peaks_per_length


def get_pvalue_eclip_peaks(length, rnbp_name, peaks_per_length_observed, num_simulations, df_lines_in_introns_coverage, eclip):
	## going to do this once for each protein
	results = []
	for i in range(0, num_simulations):
		results += [line_sampler(df_lines_in_introns_coverage = df_lines_in_introns_coverage, eclip = eclip, length = length)]
	### fit beta to results
	a1,b1,loc1,scale1 = scipy.stats.beta.fit(results, floc=0, fscale=1)
	x = np.linspace(0,1,1000)
	y = scipy.stats.beta.pdf(x, a1, b1)
	beta_pvalue = 1 - scipy.stats.beta.cdf(peaks_per_length_observed, a1, b1)
	# plt.plot(x, y)
	# plt.show()
	# plt.hist(results,bins=20)
	# plt.show()
	print(rnbp_name, "pval: ", beta_pvalue, "done")
	return rnbp_name, beta_pvalue


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
eclip = pybedtools.BedTool.from_dataframe(df_eclip)
#################
### create the vlincs data frame and the vlincs bed object
# df_vlincs = pd.read_csv(skewed_vlincs,sep="\t",names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
# 																"binom_pval","fdr_pval","fdr_reject","total_reads","skew"],
# 															dtype={"chrom":str,"start":int, "stop":int,"name_lncrna":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
# 																"binom_pval":float,"fdr_pval":float,"fdr_reject":str,"total_reads":int,"skew":float}) 
# vlincs = pybedtools.BedTool.from_dataframe(df_vlincs.loc[:,["chrom","start","stop","name_lncrna","rpkm","strand","hap1_reads","hap2_reads","binom_pval","skew"]])
# ################
# ### Get all the L1 elements within vlincs, non strand specific, and create dataframe
# lines_in_vlincs = lines.intersect(vlincs,f=0.5,wa=True)
# df_lines_in_vlincs = lines_in_vlincs.to_dataframe(names=["chrom","start","stop","name","score","strand"])
# ###############
# ## get all the eclip binding within vlincs, strand specific
# binding_in_vlincs = eclip.intersect(vlincs,f=0.5,wa=True,s=True,wb=True) # this binding does need to be stranded
# ## make a dataframe that has eclip within vlincs
# df_binding_in_vlincs = binding_in_vlincs.to_dataframe(names=["chrom_eclip", "start_eclip", "stop_eclip", "name_eclip","score_eclip", "strand_eclip", "log_fold_enrichment","log_pvalue",
# 															"chrom_vlinc","startvlinc","stop_vlinc","name_vlinc","rpkm","strand_vlinc","hap1_reads","hap2_reads",
# 																"binom_pval","skew"],
# 													dtype={"chrom_eclip":str, "start_eclip":int, "stop_eclip":int,"name_eclip":str, "score_eclip":float, "strand_eclip":str,"log_fold_enrichment":float,"log_pvalue":float,
# 													"chrom_vlinc":str,"start_vlinc":int, "stop_vlinc":int,"name_vlinc":str,"rpkm":float,"strand_vlinc":str,"hap1_reads":int,"hap2_reads":int,"binom_pval":float,"skew":float})
# ## now intersect eclip
# binding_in_lines_in_vlincs = lines_in_vlincs.intersect(binding_in_vlincs,wa=True,wb=True)
# df_binding_in_lines_in_vlincs = binding_in_lines_in_vlincs.to_dataframe(names = ["chrom_line","start_line","stop_line","name_line","score_line","strand_line",
# 																				"chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip","strand_eclip","log_fold_enrichment","log_pvalue",
# 																				"chrom_vlinc","start_vlinc","stop_vlinc","name_vlinc", "rpkm_vlinc","strand_vlinc","hap1_reads","hap2_reads","binom_pval","skew"])
##############
## get all lines within introns
lines_in_introns = lines.intersect(introns,f=0.5,wa=True)
df_lines_in_introns = lines_in_introns.to_dataframe(names = ["chrom_line","start_line","stop_line","name_line","score_line","strand_line"])
## get all eclip peaks within introns, strand specific, and make df
binding_in_introns = eclip.intersect(introns,f=0.5,wa=True,s=True,wb=True) # must be stranded
df_binding_in_introns = binding_in_introns.to_dataframe(names=["chrom_eclip", "start_eclip", "stop_eclip", "name_eclip","score_eclip", "strand_eclip", "log_fold_enrichment","log_pvalue",
															"chrom_intron","start_intron","stop_intron","name_intron","score_intron","strand_intron"],
													dtype={"chrom_eclip":str, "start_eclip":int, "stop_eclip":int,"name_eclip":str, "score_eclip":float, "strand_eclip":str,"log_fold_enrichment":float,"log_pvalue":float,
													"chrom_intron":str,"start_intron":int, "stop_intron":int,"name_intron":str,"score_intron":int,"strand_intron":str})

df_binding_in_introns_subtract_lines = eclip.intersect(introns.subtract(lines_in_introns),wa=True,wb=True,s=True,f=0.5).to_dataframe(names=["chrom_eclip", "start_eclip", "stop_eclip", "name_eclip","score_eclip", "strand_eclip", "log_fold_enrichment","log_pvalue",
															"chrom_intron","start_intron","stop_intron","name_intron","score_intron","strand_intron"],
													dtype={"chrom_eclip":str, "start_eclip":int, "stop_eclip":int,"name_eclip":str, "score_eclip":float, "strand_eclip":str,"log_fold_enrichment":float,"log_pvalue":float,
													"chrom_intron":str,"start_intron":int, "stop_intron":int,"name_intron":str,"score_intron":int,"strand_intron":str})

# df_lines_in_introns_coverage = lines_in_introns.coverage(binding_in_introns,F=0.5).to_dataframe(names=["chrom", "start", "stop", "line_name","score",
# 																									"strand","eclip_count", "num_bases_covered",
# 																									"line_length", "fraction_bound"],
# 																							dtype={"chrom":str,"start":int,"stop":int,"line_name":str,
# 																							"score":int,"strand":str,"eclip_count":int,
# 																							"num_bases_covered":int,"line_length":int,"fraction_bound":float})
binding_in_lines_in_introns = lines_in_introns.intersect(binding_in_introns,wa=True,wb=True)
df_binding_in_lines_in_introns = binding_in_lines_in_introns.to_dataframe(names = ["chrom_line","start_line","stop_line","name_line","score_line","strand_line",
																				"chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip","strand_eclip","log_fold_enrichment","log_pvalue",
																				"chrom_intron","start_intron","stop_intron","name_intron", "score_intron","strand_intron"])
#####
## KEEP BELOW.
# Just the total summed length of all the lines in all the vlincs
# total_line_length_in_vlincs = (df_lines_in_vlincs["stop"] - df_lines_in_vlincs["start"]).sum()

# pool = mp.Pool()
# print(mp.cpu_count(),"cpu's available")
# results = pool.starmap(get_pvalue_eclip_peaks,[(10**5,
# 													x,
# 													len(df_binding_in_lines_in_vlincs[df_binding_in_lines_in_vlincs["name_eclip"] == x].index) / total_line_length_in_vlincs,
# 													30,
# 													df_lines_in_introns_coverage,
# 													pybedtools.BedTool.from_dataframe(df_eclip[df_eclip["name"] == x]),
# 													 ) for x in eclip_experiments])

# with open("rnbp_l1_binding_in_asars_results.txt", "w") as f:
# 	writer = csv.writer(f,delimiter="\t") 
# 	for i in range(len(results)):
# 		writer.writerow(results[i])

###################################
## check to see peaks / bp of L1 within asars vs peaks / bp of non-L1
# vlincs_subtract_lines_eclip_coverage = vlincs.subtract(lines_in_vlincs).coverage(binding_in_vlincs).to_dataframe(names=["chrom_vlinc","start_vlinc","stop_vlinc","name_vlinc", "rpkm_vlinc","strand_vlinc",
# 																														"eclip_count","num_bases_covered","vlinc_length","fraction_bound"])
# lines_coverage = lines_in_vlincs.coverage(binding_in_vlincs).to_dataframe(names=["chrom_line","start_line","stop_line","name_line","score_line","strand_line",
# 																				"eclip_count","num_bases_covered","line_length","fraction_bound"])

# vlincs_peaks_per_mb = vlincs_subtract_lines_eclip_coverage["eclip_count"].sum() / vlincs_subtract_lines_eclip_coverage["vlinc_length"].sum() * 10**6
# lines_peaks_per_mb = lines_coverage["eclip_count"].sum() / lines_coverage["line_length"].sum() * 10**6

# print("vlincs without lines peaks per base: ",vlincs_peaks_per_mb, "lines in vlincs peaks per base: ", lines_peaks_per_mb)
###########
## check to see FOR EACH PROTEIN if it binds more in vlinc L1s vs vlinc non-L1
# for rnbp in eclip_experiments:
# 	eclip_tmp = pybedtools.BedTool.from_dataframe(df_binding_in_vlincs[df_binding_in_vlincs["name_eclip"]==rnbp])
# 	vlincs_subtract_lines_eclip_coverage = vlincs.subtract(lines_in_vlincs).coverage(eclip_tmp).to_dataframe(names=["chrom_vlinc","start_vlinc","stop_vlinc","name_vlinc", "rpkm_vlinc","strand_vlinc",
# 																														"eclip_count","num_bases_covered","vlinc_length","fraction_bound"])
# 	lines_coverage = lines_in_vlincs.coverage(eclip_tmp).to_dataframe(names=["chrom_line","start_line","stop_line","name_line","score_line","strand_line",
# 																				"eclip_count","num_bases_covered","line_length","fraction_bound"])
# 	vlincs_peaks_per_100kb = vlincs_subtract_lines_eclip_coverage["eclip_count"].sum() / vlincs_subtract_lines_eclip_coverage["vlinc_length"].sum() * 10**5
# 	lines_peaks_per_100kb = lines_coverage["eclip_count"].sum() / lines_coverage["line_length"].sum() * 10**5
# 	print(rnbp,lines_peaks_per_100kb,vlincs_peaks_per_100kb)

# #############
# ## check to see antisense vs sense binding
# for rnbp in np.unique(df_binding_in_lines_in_vlincs["name_eclip"].values):
# 	antisense_binding = len(df_binding_in_lines_in_vlincs[(df_binding_in_lines_in_vlincs["strand_line"] != df_binding_in_lines_in_vlincs["strand_vlinc"]) & (df_binding_in_lines_in_vlincs["name_eclip"]==rnbp)].index)
# 	sense_binding = len(df_binding_in_lines_in_vlincs[(df_binding_in_lines_in_vlincs["strand_line"] == df_binding_in_lines_in_vlincs["strand_vlinc"]) & (df_binding_in_lines_in_vlincs["name_eclip"]==rnbp)].index)
# 	print(rnbp, antisense_binding + sense_binding, antisense_binding / (antisense_binding + sense_binding), sense_binding / (antisense_binding + sense_binding) )

# ############

## Look at intronic L1 non-L1 ratio same way as we did with ASARs but in introns.																																])
df_introns_subtract_lines = introns.subtract(lines).to_dataframe(names=["chrom_intron","start_intron","stop_intron","name_intron", "score_intron","strand_intron"])
for rnbp in eclip_experiments:
	### pybedtools stupidly writes tmp files that don't get deleted until the end of the entire script 
	### even if you overwrite the variable in python.
	### however, you are stupid for doing this loop when you could just do it once before the loop
	# eclip_tmp = pybedtools.BedTool.from_dataframe(df_binding_in_introns[df_binding_in_introns["name_eclip"]==rnbp])
	# def be better off working with the big intersect file, then indexing by name and counting rows, instead of doing .coverage()
	tmp = df_binding_in_introns_subtract_lines[df_binding_in_introns_subtract_lines["name_eclip"]==rnbp]
	tmp2 = df_binding_in_lines_in_introns[df_binding_in_lines_in_introns["name_eclip"]==rnbp]
	introns_peaks_per_100kb = len(tmp.index) / (df_introns_subtract_lines["stop_intron"] - df_introns_subtract_lines["start_intron"]).sum() * 10**5

	lines_peaks_per_100kb = len(tmp2.index) / (df_lines_in_introns["stop_line"] - df_lines_in_introns["start_line"]).sum() * 10**5
	print(rnbp,lines_peaks_per_100kb,introns_peaks_per_100kb)



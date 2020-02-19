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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import statsmodels.stats.multitest as mt


## example of how to call this script

"""
python allele.rnaseq.protein.coding.py --dfplus ../data/gm12878.rep1.hg19Aligned.out.plus.overlap.platinum.haplotypes.bed 
	--dfminus ../data/gm12878.rep1.hg19Aligned.out.minus.overlap.platinum.haplotypes.bed 
	--window_file_plus ../data/gm12878.rep1.hg19Aligned.out.samtool.rmdup.plus.1000.10000.50000.vlinc.discovery.bed 
	--window_file_minus ../data/gm12878.rep1.hg19Aligned.out.samtool.rmdup.minus.1000.10000.50000.vlinc.discovery.bed 
	--ref_seq_plus ../annotation.files/refseq.unique.genes.merged.plus.nochr.hg19.bed 
	--ref_seq_minus ../annotation.files/refseq.unique.genes.merged.minus.nochr.hg19.bed 
	--out_directory ./
"""

ratios=[249250621/249250621,
	243199373/249250621,
	198022430/249250621,
	191154276/249250621,
	180915260/249250621,
	171115067/249250621,
	159138663/249250621,
	146364022/249250621,
	141213431/249250621,
	135534747/249250621,
	135006516/249250621,
	133851895/249250621,
	115169878/249250621,
	107349540/249250621,
	102531392/249250621,
	90354753/249250621,
	81195210/249250621,
	78077248/249250621,
	59128983/249250621,
	63025520/249250621,
	48129895/249250621,
	51304566/249250621,
	155270560/249250621]

lengths = [249250621,
	243199373,
	198022430,
	191154276,
	180915260,
	171115067,
	159138663,
	146364022,
	141213431,
	135534747,
	135006516,
	133851895,
	115169878,
	107349540,
	102531392,
	90354753,
	81195210,
	78077248,
	59128983,
	63025520,
	48129895,
	51304566,
	155270560]

centromere = {"1":124535434,
				"2":95326171,
				"3":93504854,
				"4":52660117,
				"5":49405641,
				"6":61830166,
				"7":61054331,
				"8":46838887,
				"9":50367679,
				"X":61632012,
				"Y":13104553,
				"10":42254935,
				"11":54644205,
				"12":37856694,
				"13":19000000,
				"14":19000000,
				"15":20000000,
				"16":38335801,
				"17":25263006,
				"18":18460898,
				"19":27681782,
				"20":29369569,
				"21":14288129,
				"22":16000000}
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/data/cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

def add_binom_pval(df):
	df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_reads"],
							row["hap1_reads"]+row["hap2_reads"],
							p=0.5,
							alternative="two-sided"), # v slow for some reason 
							axis=1)

	df["fdr_pval"]=mt.multipletests(pvals=df["binom_pval"], 
								alpha=0.01,
								method="fdr_bh")[1]
	df["fdr_reject"] =  mt.multipletests(pvals=df["binom_pval"], 
									alpha=0.01,
									method="fdr_bh")[0]
	return
def get_windows_refseq(window_file, read_counts_file, is_file=True): 
	
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles

	a = a.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=a,b=b,c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	### invalid int for literal...
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "name", "score", "strand", 
										"hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"name":str, "score":float, "strand":str,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df

def get_windows_vlincs(window_file, read_counts_file, is_file=True): 
	
	a = pybedtools.BedTool(window_file) # change this to bed file of previously determined windows of interest
	b = pybedtools.BedTool(read_counts_file) ## read counts at alleles

	a = a.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	b = b.sort(faidx="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai")
	#/home/groups/Spellmandata/heskett/refs/
	#windows=a.window_maker(g="/Users/heskett/replication_rnaseq/scripts/hg38.10x.nochr.fa.fai",
	#					w=length,s=length/2)

	c = pybedtools.BedTool()
	## genome file to specify sort order

	window_read_counts = c.map(a=a,b=b,c=[6,7],o="sum",g="/Users/heskett/replication_rnaseq/annotation.files/human_g1k_v37.fasta.fai") ## this can fail?? somehow dos2unix helped?
	### invalid int for literal...
	df =  window_read_counts.to_dataframe(names=["chrom", "start", "stop", "name", "score", "strand", 
										"fraction_l1","hap1_reads", "hap2_reads"],
										dtype={"chrom":str, "start":int, "stop":int,
										"name":str, "score":float, "strand":str, "fraction_l1":float,
										"hap1_reads":str, "hap2_reads":str})
	df = df[ (df["hap1_reads"]!=".") & (df["hap2_reads"]!=".") ]
	df["hap1_reads"] = df["hap1_reads"].astype(int)
	df["hap2_reads"] = df["hap2_reads"].astype(int)
	df["chrom"] = df["chrom"].astype(str)
	return df




if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	parser.add_argument("--dfplus",
		type=str,
		metavar="[df of allele counts plus]",
		required=True,
		help="in bed format")
	parser.add_argument("--dfminus",
		type=str,
		metavar="[df of allele counts minus]",
		required=True,
		help="in bed format")

	parser.add_argument("--out_directory",
		type=str,
		metavar="[out directory]",
		required=True,
		help="full path to output results")
	parser.add_argument("--ref_seq_plus",
		type=str,
		metavar="[bed file of genes with gene names]",
		required=False,
		help="use protein coding gene windows")
	parser.add_argument("--ref_seq_minus",
		type=str,
		metavar="[bed file of genes with gene names]",
		required=False,
		help="use protein coding gene windows")
	parser.add_argument("--window_file_plus",
		type=str,
		metavar="[bed file of windows]",
		required=False,
		help="use previously determined windows instead of tiling")
	parser.add_argument("--window_file_minus",
		type=str,
		metavar="[bed file of windows]",
		required=False,
		help="use previously determined windows instead of tiling")

	arguments = parser.parse_args()

	df_plus_coding = get_windows_refseq(window_file = arguments.ref_seq_plus,
						read_counts_file=arguments.dfplus)
	df_minus_coding = get_windows_refseq(window_file = arguments.ref_seq_minus,
						read_counts_file=arguments.dfminus)
	df_combined_coding = pd.concat([df_plus_coding, df_minus_coding])
	add_binom_pval(df_combined_coding)
	df_combined_coding["total_reads"] = df_combined_coding["hap1_reads"] + df_combined_coding["hap2_reads"]

## definition of skewing here...can try a few different things.
	df_combined_coding = df_combined_coding[df_combined_coding["total_reads"]>=15] # or do at least min 1 read per kb
	df_combined_coding.loc[:,"skew"] = df_combined_coding.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
														(-x["hap2_reads"]  / x["total_reads"])	+ 0.5, axis = 1)
	skew_thresholds = np.percentile(df_combined_coding[df_combined_coding["chrom"]!="X"]["skew"],[5,95])

#### now do this for the vlincs
	df_plus_vlinc = get_windows_vlincs(window_file = arguments.window_file_plus,
						read_counts_file=arguments.dfplus)
	df_minus_vlinc = get_windows_vlincs(window_file = arguments.window_file_minus,
						read_counts_file=arguments.dfminus)
	df_combined_vlinc = pd.concat([df_plus_vlinc, df_minus_vlinc])
	add_binom_pval(df_combined_vlinc)
	df_combined_vlinc["total_reads"] = df_combined_vlinc["hap1_reads"] + df_combined_vlinc["hap2_reads"]

	df_combined_vlinc = df_combined_vlinc[df_combined_vlinc["total_reads"]>=20] # or do at least min 1 read per kb
	df_combined_vlinc.loc[:,"skew"] = df_combined_vlinc.apply(lambda x: (x["hap1_reads"]  / x["total_reads"] - 0.5) if (x["hap1_reads"] >= x["hap2_reads"]) else 
														(-x["hap2_reads"]  / x["total_reads"])	+ 0.5, axis = 1)

	df_combined_vlinc = df_combined_vlinc[(df_combined_vlinc["skew"] <= skew_thresholds[0]) | (df_combined_vlinc["skew"] >= skew_thresholds[1])]
	df_combined_vlinc = df_combined_vlinc[(df_combined_vlinc["binom_pval"] <= 0.001) & (df_combined_vlinc["chrom"] != "X")]
	print(df_combined_vlinc)

	plt.hist(df_combined_coding[df_combined_coding["chrom"]!="X"]["skew"],bins=30)
	plt.axvline(x=skew_thresholds[1],linestyle="--",c="red")
	plt.axvline(x=skew_thresholds[0],linestyle="--",c="red")
	# plt.show()
	# plt.close()

out_string = arguments.out_directory + os.path.basename(arguments.dfplus.replace(".plus.overlap.platinum.haplotypes.bed","")) + "significant.skewed.vlincs.txt"
df_combined_vlinc =  df_combined_vlinc.sort_values(by = ["chrom","start"])
df_combined_vlinc.to_csv(out_string,index=None,sep="\t")



out_string_bed =  arguments.out_directory + os.path.basename(arguments.dfplus.replace(".plus.overlap.platinum.haplotypes.bed","")) + "significant.skewed.vlincs.bed"
df_combined_vlinc.to_csv(out_string_bed,index=None,header=None, sep="\t")




# ## do the same thing with RT using protein coding autosomal genes?
# ## way too noisy to do this at gene level without smoothing first...
# 	df_early =	pd.read_csv("4e.combined.overlap.na12878.hg19.haplotype.resolved.counts.sorted.refseq.protein.windows.bed",sep="\t",
# 		names=["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
# 		dtype={"chrom":str, "start":int, "stop":int, "name":str, "score":float, "strand":str,"hap1_counts":str, "hap2_counts":str})

# 	df_late =	pd.read_csv("4l.combined.overlap.na12878.hg19.haplotype.resolved.counts.sorted.refseq.protein.windows.bed",sep="\t",
# 		names=["chrom","start","stop","name","score","strand","hap1_counts","hap2_counts"],
# 		dtype={"chrom":str, "start":int, "stop":int, "name":str, "score":float, "strand":str,"hap1_counts":str, "hap2_counts":str})

	
# 	df_early = df_early[ (df_early["hap1_counts"]!=".") & (df_early["hap2_counts"]!=".") ]
# 	df_early["hap1_counts"] = df_early["hap1_counts"].astype(int)
# 	df_early["hap2_counts"] = df_early["hap2_counts"].astype(int)
# 	df_early = df_early[df_early["hap1_counts"] + df_early["hap2_counts"] >= 15]


# 	df_late = df_late[ (df_late["hap1_counts"]!=".") & (df_late["hap2_counts"]!=".") ]
# 	df_late["hap1_counts"] = df_late["hap1_counts"].astype(int)
# 	df_late["hap2_counts"] = df_late["hap2_counts"].astype(int)
# 	df_late = df_late[df_late["hap1_counts"] + df_late["hap2_counts"] >= 15]

# 	df_early.loc[:,"logR"] = np.log2( (df_early["hap1_counts"]+1) / (df_early["hap2_counts"]+1) )
# 	df_late.loc[:,"logR"] = np.log2( (df_late["hap1_counts"]+1) / (df_late["hap2_counts"]+1) )
# 	plt.hist(df_early[df_early["chrom"]!="X"]["logR"],bins=100)
# 	plt.show()
# 	plt.close()
# 	print(df_early[df_early["chrom"]!="X"]["logR"].describe(percentiles=[0.001,0.05,0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99,.999]))

# 	plt.hist(df_late[df_late["chrom"]!="X"]["logR"],bins=100)
# 	plt.show()
# 	plt.close()
# 	print(df_late[df_late["chrom"]!="X"]["logR"].describe(percentiles=[0.001,0.05,0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99,.999]))
# 	print(df_late[(df_late["chrom"]!="X") & (abs(df_late["logR"]) >= 4)])

# #####
# arm_dict = get_arms(cytoband)

# df_windows = pd.read_csv("4e.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
# 						names=["chrom","start","stop","hap1_counts","hap2_counts"],
# 						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
# df2_windows = pd.read_csv("4l.100kb.slide.25kb.haplotype.counts.bed",sep="\t",header=None,index_col=None,
# 						names=["chrom","start","stop","hap1_counts","hap2_counts"],
# 						dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
# df_windows = df_windows[df_windows["hap1_counts"] + df_windows["hap2_counts"] >= 15]
# df2_windows = df2_windows[df2_windows["hap1_counts"] + df2_windows["hap2_counts"] >= 15]

# df_windows.loc[:,"logR"] = np.log2( (df_windows["hap1_counts"]+1) / (df_windows["hap2_counts"]+1) )

# df2_windows.loc[:,"logR"] = np.log2( (df2_windows["hap1_counts"]+1) / (df2_windows["hap2_counts"]+1) )

# df_windows["arm"] = df_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
# df2_windows["arm"] = df2_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)


# print(np.percentile(df_windows[df_windows["chrom"]!="X"]["logR"],[5,95]))
# print(np.percentile(df2_windows[df2_windows["chrom"]!="X"]["logR"],[5,95]))

# df_windows_asynch = df_windows[(df_windows["chrom"]!="X")]["logR"]




# plt.hist(df_windows["logR"],bins=30)
# plt.show()
# plt.close()
# plt.hist(df2_windows["logR"],bins=30)
# plt.show()
# plt.close()
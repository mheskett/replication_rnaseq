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


### repliseq
df_asynch_4 = pd.read_csv("4el.asynchronous.regions.bed",sep="\t",
	names=["chrom",  "start","stop", "hap1_reads", "hap2_reads","pval","type"])
df_asynch_5 = pd.read_csv("5el.asynchronous.regions.bed",sep="\t",
	names=["chrom",  "start","stop", "hap1_reads", "hap2_reads","pval","type"])
#### all windows
df_pool_tiling = pd.read_csv("gm12878.rep1.hg19Aligned.out.tiling.all.bed",sep="\t",
	names = ["chrom",  "start","stop", "hap1_reads", "hap2_reads", "strand" ,"binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "strand":str, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})
df_4_tiling = pd.read_csv("gm12878.4x.hg19Aligned.out.tiling.all.bed",sep="\t",
	names = ["chrom",  "start","stop", "hap1_reads", "hap2_reads", "strand" ,"binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "strand":str, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})
df_5_tiling = pd.read_csv("gm12878.5x.hg19Aligned.out.tiling.all.bed",sep="\t",
	names = ["chrom",  "start","stop", "hap1_reads", "hap2_reads", "strand" ,"binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "strand":str, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})
######### vlincs only
df_pool = pd.read_csv("gm12878.rep1.hg19Aligned.outall.vlincs.bed",sep="\t",
	names = ["chrom",  "start","stop", "name", "score", "strand","fraction_l1", "hap1_reads", "hap2_reads", "binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "name":str, "score":str, "strand":str,"fraction_l1":float, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})

df_clone4 = pd.read_csv("gm12878.4x.hg19Aligned.outall.vlincs.bed",sep="\t",
	names = ["chrom",  "start","stop", "name", "score", "strand","fraction_l1", "hap1_reads", "hap2_reads", "binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "name":str, "score":str, "strand":str,"fraction_l1":float, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})

df_clone5 = pd.read_csv("gm12878.5x.hg19Aligned.outall.vlincs.bed",sep="\t",
	names = ["chrom",  "start","stop", "name", "score", "strand","fraction_l1", "hap1_reads", "hap2_reads", "binom_pval", "fdr_pval", "fdr_reject", "total_reads","skew"],
	dtype= {"chrom":str,  "start":int,"stop":int, "name":str, "score":str, "strand":str,"fraction_l1":float, "hap1_reads":int, "hap2_reads":int, "binom_pval":float, "fdr_pval":float, "fdr_reject":str, "total_reads":int,"skew":float})

###
## used this command
### cat gm12878.4x.hg19Aligned.outall.vlincs.bed gm12878.5x.hg19Aligned.outall.vlincs.bed | bedtools sort -i stdin | bedtools merge -s -i stdin -c 6,14,14 -o distinct,count,distinct > all.vlincs.clone4.vs.clone5.skew.variance.bed
df_switchers = pd.read_csv("all.vlincs.clone4.vs.clone5.skew.variance.bed",sep="\t",header=None,names=["chrom","start","stop","strand","count","skews"],	dtype={"chrom":str,"start":int,"stop":int,"strand":str,"count":int})
# skews=[]
# for index,row in df_switchers.iterrows():
# 	skews += [[float(i) for i in row["skews"].split(',')]]
# df_switchers["skews"] = skews
# df_switchers = df_switchers[df_switchers["count"]>1]
# df_switchers.loc[:,"max_skew_distance"] =df_switchers.apply(lambda x: abs(np.min(x["skews"]) - np.max(x["skews"])),axis=1 )
# df_switchers.sort_values(by=["max_skew_distance"], ascending=False).to_csv("all.vlincs.clone4.vs.clone5.skew.switching.bed",sep="\t",header=None,index=None)
# plt.hist(df_switchers["max_skew_distance"],bins=30)
# plt.show()
# plt.close()
# ####


for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))
	## also plot tiling in gm12878 pool
	plt.scatter(df_pool_tiling[df_pool_tiling["chrom"]==chromosomes[i]]["start"],df_pool_tiling[df_pool_tiling["chrom"]==chromosomes[i]]["skew"],
				label="gm12878 pool windows",s=20,c="blue",lw=0.2,edgecolor="black",alpha=0.6 )
	# plt.scatter(df_4_tiling[df_4_tiling["chrom"]==chromosomes[i]]["start"],df_4_tiling[df_4_tiling["chrom"]==chromosomes[i]]["skew"],
	# 			label="gm12878 4 windows",s=20,c="orange",lw=0.2,edgecolor="black" ,alpha=0.6)
	# plt.scatter(df_5_tiling[df_5_tiling["chrom"]==chromosomes[i]]["start"],df_5_tiling[df_5_tiling["chrom"]==chromosomes[i]]["skew"],
	# 			label="gm12878 5 windows",s=20,c="green",lw=0.2,edgecolor="black",alpha=0.6 )
	# ### repliseq here
	# for index, row in df_asynch_4[df_asynch_4["chrom"]==chromosomes[i]].iterrows():
	# 	if (row["type"] == "early_hap1") or (row["type"] == "late_hap2"):
	# 		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="blue", alpha=0.15)
	# 	if (row["type"] == "late_hap1") or (row["type"] == "early_hap2"):
	# 		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="red", alpha=0.15)
	# # for index, row in df_asynch_5[df_asynch_5["chrom"]==chromosomes[i]].iterrows():
	# # 	if (row["type"] == "early_hap1") or (row["type"] == "late_hap2"):
	# # 		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="blue", alpha=0.15)
	# # 	if (row["type"] == "late_hap1") or (row["type"] == "early_hap2"):
	# # 		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="red", alpha=0.15)



	### significant vlincs here
	# plt.scatter(df_pool[df_pool["chrom"]==chromosomes[i]]["start"],df_pool[df_pool["chrom"]==chromosomes[i]]["skew"],label="gm12878 pool",s=20,c="green",lw=0.2,edgecolor="black")
	# plt.scatter(df_clone4[df_clone4["chrom"]==chromosomes[i]]["start"],df_clone4[df_clone4["chrom"]==chromosomes[i]]["skew"],label="gm12878 clone4",s=20,c="orange",lw=0.2,edgecolor="black")
	# plt.scatter(df_clone5[df_clone5["chrom"]==chromosomes[i]]["start"],df_clone5[df_clone5["chrom"]==chromosomes[i]]["skew"],label="gm12878 clone5",s=20,c="pink",lw=0.2,edgecolor="black")

	plt.axhline(y=0,linestyle="--",color="black")
	plt.legend()
	# plt.grid()
	# for index,row in df_switchers[df_switchers["chrom"]==chromosomes[i]].iterrows():
	# 	if row["max_skew_distance"] > 0.4:
	# 		plt.axvspan(xmin=row["start"], xmax=row["stop"],color="red",alpha=0.2)
	ax.margins(x=0,y=0)
	plt.show()
	# plt.savefig("vlinc_only_"+chromosomes[i],dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()




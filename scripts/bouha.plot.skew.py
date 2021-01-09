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
import scipy.stats
import seaborn as sns
import statsmodels.api as sm
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="scripts to combine samples")
	parser.add_argument("--df_expression",
		type=str,
		metavar="[df of allele counts plus]",
		required=True,
		help="in bed format")
	parser.add_argument("--dfminus",
		type=str,
		metavar="[df of allele counts minus]",
		required=True,
		help="in bed format")





	## list of imprinted genes
	df_imprinted = pd.read_table("/Users/heskett/replication_rnaseq/scripts/imprinted.genes.fixed.bed",
		sep="\t", names=["chrom","start","stop","gene_name"],dtype={"chrom":str},
		header=None,index_col=None)

	## list of chromosome names
	chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
					"13","14","15","16","17","18","19","20","21","22","X"]
	arms = ["p","q"]
	#### for arm level data to skip over centromeres				
	cytoband = pd.read_table("/Users/heskett/replication_rnaseq/data/cytoband.nochr.hg19.bed",sep="\t",
								names =["chrom","start","stop","arm","band"])
	chromosome_length = {"1":249250621,
	"2":243199373,
	"3":198022430,
	"4":191154276,
	"5":180915260,
	"6":171115067,
	"7":159138663,
	"8":146364022,
	"9":141213431,
	"10":135534747,
	"11":135006516,
	"12":133851895,
	"13":115169878,
	"14":107349540,
	"15":102531392,
	"16":90354753,
	"17":81195210,
	"18":78077248,
	"19":59128983,
	"20":63025520,
	"21":48129895,
	"22":51304566,
	"X":155270560}

	def get_arms(cytoband):
		## given a data frame with genome elements, add the arm information to a new column
		arm_dict = {}
		for i in range(len(chromosomes)):
			# should be (p end, q end)
			arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
			cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
		return arm_dict

	arm_dict = get_arms(cytoband)


	df_tiling_expression = pd.read_csv("gm12878.5x.hg19Aligned.out.tiling.all.bed",sep="\t",
							names=["chrom","start","stop","hap1_counts","hap2_counts","strand","pval","qval","reject","total_reads","skew"])

	df = pd.read_csv("gm12878.5x.hg19Aligned.outsignificant.skewed.vlincs.bed",sep="\t",
						names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
						dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})

	df_windows = df_windows[df_windows["hap1_counts"] + df_windows["hap2_counts"] >= 15]
	df2_windows = df2_windows[df2_windows["hap1_counts"] + df2_windows["hap2_counts"] >= 15]

	df_windows.loc[:,"logR"] = np.log2( (df_windows["hap1_counts"]+1) / (df_windows["hap2_counts"]+1) )

	df2_windows.loc[:,"logR"] = np.log2( (df2_windows["hap1_counts"]+1) / (df2_windows["hap2_counts"]+1) )

	df_windows["arm"] = df_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
	df2_windows["arm"] = df2_windows.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)

	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=20)
	plt.rc('figure', titlesize=15)

	# ### get asynchronous distributions
	# early_thresholds = np.percentile(df_windows[df_windows["chrom"]!="X"]["logR"],[5,95])
	# late_thresholds = np.percentile(df2_windows[df2_windows["chrom"]!="X"]["logR"],[5,95])


	# asynch_early =  df_windows[(df_windows["chrom"] == chromosomes[i]) & ((df_windows["logR"] >= early_thresholds[1]) | (df_windows["logR"] <= early_thresholds[0]))]
	# asynch_late =  df2_windows[(df2_windows["chrom"] == chromosomes[i]) & ((df2_windows["logR"] >= late_thresholds[1]) | (df2_windows["logR"] <= late_thresholds[0]))]

	# print(asynch_early)

	for i in range(len(chromosomes)):
		f,ax = plt.subplots(figsize=(12,2))
		## raw data repliseq
		# ax.scatter(df_windows[df_windows["chrom"]==chromosomes[i]]["start"], df_windows[df_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="pink",label="early hap1/hap2",alpha=0.6)
		# ax.scatter(df2_windows[df2_windows["chrom"]==chromosomes[i]]["start"], df2_windows[df2_windows["chrom"]==chromosomes[i]]["logR"],s=5,c="green",label="late hap1/hap2",alpha=0.6)
		## smoothed repliseq
		## tiling windows rna-seq
		df_tiling_expression.loc[:,"strand_color"] = df_tiling_expression.apply(lambda x: "pink" if x["strand"]=="+" else "blue", axis = 1)
		ax.scatter(df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["start"],
				df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["skew"],
				label="all RNA Expression Skew",lw=0.2,edgecolor="black",
				c=df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["strand_color"],zorder=3,s=8,alpha=0.4)
		ax.set_ylim([-.52,.52])
		ax.set_yticks([])
		ax.set_xlim([0, chromosome_length[chromosomes[i] ] ] )
		# ax.set_xticks([])
		ax.axhline(y=0,linestyle="--",c="black")

		if chromosomes[i]=="X":
			ax3.set_ylim([-2.5,2.5])
		else:
			ax3.set_ylim([-2,2])
		ax3.set_yticks([])
		ax3.set_xlim([0, chromosome_length[chromosomes[i] ] ] )

		ax3.axhline(y=0, linestyle="--", c="black")
		ax3.margins(x=0,y=0)

		### lncrna goes on top
		ax2 = ax.twinx()
		ax2.scatter(df[df["chrom"]==chromosomes[i]]["start"], df[df["chrom"]==chromosomes[i]]["skew"],c="orange",
			label="lncRNA Expression Skew",lw=0.2,edgecolor="black",zorder=2,s=15)
		ax2.set_ylim([-.5,.5])
		ax2.set_xlim([0, chromosome_length[chromosomes[i] ] ] )

		ax2.set_yticks([])
		# ax2.set_xticks([])

		ax.margins(x=0,y=0)
		ax2.margins(x=0,y=0)
		# ax.ticklabel_format(style='plain')
		# ax2.ticklabel_format(style='plain')
		# ax3.ticklabel_format(style='plain')
		# ax2.set_xticks([])
		### add known imprinted genes as green dots
		ax4 = ax.twinx()
		print(df_imprinted[df_imprinted["chrom"]==chromosomes[i]])
		for index, row in df_imprinted[df_imprinted["chrom"]==chromosomes[i]].iterrows():
			ax4.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="green", alpha=0.8)

		# plt.suptitle("chromosome " + chromosomes[i])	# plt.show()
		# plt.show()


		plt.savefig("5x_skew"+chromosomes[i]+".png", dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
		plt.close()
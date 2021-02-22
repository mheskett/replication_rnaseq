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
from sys import argv
 

def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict

## list of chromosome names
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
#### for arm level data to skip over centromeres				
cytoband = pd.read_table("cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)

## expression combined file "tiling.all.bed"
df_tiling_expression = pd.read_csv(argv[1],sep="\t",
						names=["chrom","start","stop","hap1_counts","hap2_counts",
						"strand","pval","qval","reject","total_reads","skew"])
## vlincs file "vlinc discovery skewed.bed"
df = pd.read_csv(argv[2],sep="\t",
					names= ["chrom","start","stop","name","rpkm","strand", "l1_fraction","hap1_counts","hap2_counts","pval","qval","reject","total_reads","skew"],
					dtype = {"chrom":str,"start":int,"stop":int,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_counts":int,"hap2_counts":int})
df_repli = pd.read_csv(argv[3],sep="\t",
					names = ["chrom","start","stop","pval"],
					dtype = {"chrom":str,"start":int,"stop":int,"pval":float})

for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))
	df_tiling_expression.loc[:,"strand_color"] = df_tiling_expression.apply(lambda x: "pink" if x["strand"]=="+" else "blue", axis = 1)
	ax.scatter(df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["start"],
			df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["skew"],
			label="all RNA Expression Skew",lw=0.2,edgecolor="black",c=df_tiling_expression[df_tiling_expression["chrom"]==chromosomes[i]]["strand_color"],zorder=3,s=8,alpha=0.4)
	ax.set_ylim([-.52,.52])
	# ax.set_xticks([])
	ax.axhline(y=0,linestyle="--",c="black")

	ax2 = ax.twinx()
	ax2.scatter(df[df["chrom"]==chromosomes[i]]["start"], df[df["chrom"]==chromosomes[i]]["skew"],c="orange",
		label="lncRNA Expression Skew",lw=0.2,edgecolor="black",zorder=2,s=15)
	ax2.set_ylim([-.5,.5])
	ax2.set_yticks([])
	ax.margins(x=0,y=0)
	ax2.margins(x=0,y=0)
	ax.ticklabel_format(style='plain')
	ax2.ticklabel_format(style='plain')
	# ax3.ticklabel_format(style='plain')
	# ax2.set_xticks([])

	# plt.suptitle("chromosome " + chromosomes[i])	# plt.show()
	tmp_repli = df_repli[df_repli["chrom"]==chromosomes[i]]
	sig_regions_repli = tmp_repli[tmp_repli["pval"] <= 0.05]
	for index,row in sig_regions_repli.iterrows():
		ax.axvspan(xmin=row["start"], xmax=row["stop"], facecolor="green", alpha=0.8)
	plt.savefig("".join(argv[1].split(".")[0:3]) +".windows.and.vlincs.TEST"+chromosomes[i]+".png", dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
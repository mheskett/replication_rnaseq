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

def get_age(x,age_dict):
	return age_dict[x] if x in age_dict.keys() else get_age(x[:-1],age_dict)

line_age = pd.read_csv("/Users/heskett/replication_rnaseq/annotation.files/warburton.line.age.ranking.txt",
						names=["line_name","evolutionary_age_ranking"],sep="\t")
all_ucsc_lines = pd.read_csv("/Users/heskett/replication_rnaseq/annotation.files/all.l1.elements.txt",header=None,sep="\t",names=["line_name"])

all_lines_df = pd.read_csv("/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed",sep="\t",names=["chrom","start","stop","name","score","strand"],
							header=None,dtype={"chrom":str,"start":int,"stop":int,"name":str,"score":int,"strand":str})

unique,counts = np.unique(all_lines_df["name"],return_counts=True)
num_l1s_each_type = dict(zip(unique,counts))

num_l1s_df = pd.DataFrame({"line_name":unique,"number_occurences":counts})
num_l1s_df.loc[:,"relative_frequency"] = num_l1s_df["number_occurences"] / num_l1s_df["number_occurences"].sum()

# print(num_l1s_df) ## use this to sample from the multiomodal dist


asars_file = "/Users/heskett/replication_rnaseq/scripts/gm12878.rep1.hg19Aligned.outgm12878.rep1.expressed.vlincs.skewed.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"
all_vlincs_file = "/Users/heskett/replication_rnaseq/scripts/gm12878.rep1.hg19Aligned.outgm12878.rep1.expressed.vlincs.all.bed"
asars = pybedtools.BedTool(asars_file)
lines = pybedtools.BedTool(lines_file)
all_vlincs = pybedtools.BedTool(all_vlincs_file)
df = asars.intersect(lines,wa=True,wb=True).\
					to_dataframe(names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
														"binom_pval","fdr_pval","fdr_reject","total_reads","skew","chrom_line",
														"start_line","stop_line","name_line","score_line","strand_line"],
														dtype={"chrom":str,"start":int,"stop":int,"name_lncrna":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
														"binom_pval":float,"fdr_pval":float,"fdr_reject":str,"total_reads":int,"skew":float,"chrom_line":str,
														"start_line":int,"stop_line":int,"name_line":str,"score_line":float,"strand_line":str})
df_all = all_vlincs.intersect(lines,wa=True,wb=True).\
					to_dataframe(names=["chrom","start","stop","name_lncrna","rpkm","strand","l1_fraction","hap1_reads","hap2_reads",
														"binom_pval","fdr_pval","fdr_reject","total_reads","skew","chrom_line",
														"start_line","stop_line","name_line","score_line","strand_line"],
														dtype={"chrom":str,"start":int,"stop":int,"name_lncrna":str,"rpkm":float,"strand":str,"l1_fraction":float,"hap1_reads":int,"hap2_reads":int,
														"binom_pval":float,"fdr_pval":float,"fdr_reject":str,"total_reads":int,"skew":float,"chrom_line":str,
														"start_line":int,"stop_line":int,"name_line":str,"score_line":float,"strand_line":str})
df_all.to_csv("gm12878.rep1.all.vlincs.intersect.lines.bed",sep="\t",header=None,index=None)
# df_not_skewed = df_all[(df_all["binom_pval"]>10**-6) | (df_all["skew"]<0.25)]

# age_dict = dict(zip(line_age.line_name, line_age.evolutionary_age_ranking))


# ###
# all_lines_df.loc[:,"line_age"] = all_lines_df.apply(lambda x: get_age(x["name"], age_dict),axis=1)
# all_lines_df.to_csv("ucsc.L1.filtered.age.hg19.bed",sep="\t",header=None,index=None)
# ###
# df.loc[:,"line_age"] = df.apply(lambda x: get_age(x["name_line"], age_dict),axis=1) # higher number is younger age
# df_not_skewed.loc[:,"line_age"] = df_not_skewed.apply(lambda x: get_age(x["name_line"], age_dict),axis=1) # higher number is younger age



# sns.set(font_scale=1.5)

# fig, ax = plt.subplots()
# sns.kdeplot(df.groupby("name_lncrna")["line_age"].mean(),shade=True,label="LINE_age_score_skewed")
# # sns.kdeplot(df_not_skewed.groupby("name_lncrna")["line_age"].mean(),shade=True,label="LINE_age_score_biallelic")

# # ax.set_ylim([0,0.19])
# # ax.set_xlim([-2,30])
# plt.show()
# plt.close()

###
### first experiment. are there more L1s of any kind (as fraction of total sequence) than 
### you would expect compared to equivolent length of intergenic sequence.
### this is like sampling from a jar with two colors of marbles.

### second experiment. for each specific LINE in each ASAR, does that ASAR contain more 
### of that specific LINE that you would expect. this is like sampling from a jar where
## there are as many colors as there are types of LINES. this can be independent
## of the above experiment.
## try this numpy.random.multinomial(n, pvals, size=None)
# n is number of lines in an asar, size is how many times you simulate, pvals is relative freq of lines in human genome

# Use R package EMT:multinomial.test




# df = df.loc[:,"line_age"] = df.apply(lambda x: age_dict[x["name_line"]],axis=1)
# print(df)
## some strategies
## 1. create "line age" score where you calculate the average (or whatever metric) evolutionary age of all the lines in an ASAR
## 2. check to see if there's more of any particular L1 subtype than you would expect based on random sampling of all of them
## based on (warburton et al plos comp bio)
## 3. Determine if ASAR L1PA2s have shared features such as deletions, substitutions etc so for this you could cluster all L1PA2s
## and label by within an ASAR or not within an ASAR
## 4. cluster ASARs based on L1 composition -- the relative frequency of each LINE can be the features.



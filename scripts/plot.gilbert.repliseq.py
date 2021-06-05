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
# hg19
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

df_repli_hap1 = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/ethan.repliseq.analysis/ebv_10_hap1_w50000.bg",
						sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","log2r"],
						dtype = {"chrom":str,"start":int,"stop":int,"log2r":float})
df_repli_hap2 = pd.read_csv("/Users/mike/replication_rnaseq/bouhassira_data/repliseq.dec.20/ethan.repliseq.analysis/ebv_10_hap2_w50000.bg",
						sep="\t",header=None,index_col=None,
						names=["chrom","start","stop","log2r"],
						dtype = {"chrom":str,"start":int,"stop":int,"log2r":float})
print(df_repli_hap1)
print(df_repli_hap2)

for i in range(len(chromosomes)):
	f,ax = plt.subplots(figsize=(12,2))

	plt.plot(df_repli_hap1[df_repli_hap1["chrom"]==chromosomes[i]]["start"],
		df_repli_hap1[df_repli_hap1["chrom"]==chromosomes[i]]["log2r"])
	plt.plot(df_repli_hap2[df_repli_hap2["chrom"]==chromosomes[i]]["start"],
		df_repli_hap2[df_repli_hap2["chrom"]==chromosomes[i]]["log2r"])
	plt.show()
	


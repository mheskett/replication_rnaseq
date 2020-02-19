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

df_early = pd.read_csv("/Users/heskett/replication_rnaseq/scripts/4e.combined.samtool.rmdup.50kb.coverage.bed",sep="\t",header=None,index_col=None,
	names=["chrom","start","stop","counts"])

df_late = pd.read_csv("/Users/heskett/replication_rnaseq/scripts/4l.combined.samtool.rmdup.50kb.coverage.bed",sep="\t",header=None,index_col=None,
	names=["chrom","start","stop","counts"])


chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
gray_chromosomes = ["1","3","5","7","9","11","13","15","17","19","21","X"]
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

f, ax = plt.subplots(1, len(chromosomes), sharex=False,
							sharey=False,
							figsize=(15,1),
							gridspec_kw={'width_ratios': ratios})
for i in range(len(chromosomes)):
	# if theres no data then just set the formatting and skip to next one
	ax[i].scatter(df_early[df_early["chrom"]==chromosomes[i]]["start"],
			np.log10(df_early[df_early["chrom"]==chromosomes[i]]["counts"]),
			s=6,
			lw=0.1,
			edgecolor="black",
			c= "blue",
			vmin=4,
			vmax=30,
			alpha=0.6)
	ax[i].scatter(df_late[df_late["chrom"]==chromosomes[i]]["start"],
			-np.log10(df_late[df_late["chrom"]==chromosomes[i]]["counts"]),
			s=6,
			lw=0.1,
			edgecolor="black",
			c= "yellow",
			vmin=4,
			vmax=30,
			alpha=0.6)

	ax[i].axhline(y=0,linestyle="--",c="black")
	ax[i].set_yticks([0,1,2,3,4])
	ax[i].set_xticks([])
	ax[i].margins(x=0,y=0)
	ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
	ax[i].axvline(x=int(centromere[chromosomes[i]]), linestyle = "--", lw = 0.5,color="black")
	ax[i].set_xlim([0,lengths[i]])
	ax[i].set_ylim()
	plt.subplots_adjust(wspace=0, hspace=0)
for i in range(len(chromosomes)):
	if chromosomes[i] in gray_chromosomes:
	    ax[i].axvspan(xmin=0, xmax=lengths[i], ymin=0, ymax=1,
	     alpha=0.2,facecolor="gray")
plt.show()
plt.close()
for i in range(len(chromosomes)):
	f, ax = plt.subplots(2,1,sharex=True)
	ax[0].scatter(df_early[df_early["chrom"]==chromosomes[i]]["start"],
			np.log10(df_early[df_early["chrom"]==chromosomes[i]]["counts"]),
			s=6,
			lw=0.1,
			edgecolor="black",
			c= "blue",
			vmin=4,
			vmax=30,
			alpha=0.6)
	ax[1].scatter(df_late[df_late["chrom"]==chromosomes[i]]["start"],
			np.log10(df_late[df_late["chrom"]==chromosomes[i]]["counts"]),
			s=6,
			lw=0.1,
			edgecolor="black",
			c= "yellow",
			vmin=4,
			vmax=30,
			alpha=0.6)
	ax[0].margins(x=0,y=0)
	ax[1].margins(x=0,y=0)
	ax[0].set_ylim([0,4])
	ax[1].set_ylim([0,4])
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.show()
	plt.close()

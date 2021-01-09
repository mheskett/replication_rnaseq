import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
"13","14","15","16","17","18","19","20","21","22","X"]
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
df_early = pd.read_csv("/Users/heskett/replication_rnaseq/data/4e.combined.counts.bed", names=["chrom","start","stop","counts"], sep="\t")
df_late = pd.read_csv("/Users/heskett/replication_rnaseq/data/4l.combined.counts.bed", names=["chrom","start","stop","counts"], sep="\t")


df_control = pd.read_csv("/Users/heskett/replication_rnaseq/repli_seq/RT_GM12878_B-Lymphocyte_Int61574576_hg19_nochr.bed",names=["chrom","start","stop","log2r"], sep="\t")
df_control2  = pd.read_csv("/Users/heskett/replication_rnaseq/repli_seq/RT_GM12878_B-Lymphocyte_Int72761980_hg19_nochr.bed",names=["chrom","start","stop","log2r"], sep="\t")

## filtering

df_combined = pd.concat([df_early,df_late["counts"]],axis=1)

df_combined.columns = ["chrom","start","stop","early_counts","late_counts"]
df_combined = df_combined.dropna(how="any",axis=0)
df_combined = df_combined[(df_combined["early_counts"] + df_combined["late_counts"] != 0)]
## smooth log2r values
## percent early seems to be less affected by batch...
df_combined.loc[:,"percent_early"] = df_combined["early_counts"] / (df_combined["early_counts"] + df_combined["late_counts"])
df_combined.loc[:,"log2r"] = np.log2 ( df_combined["early_counts"] /  df_combined["late_counts"] )

plt.rc("xtick",labelsize=15)
plt.rc("ytick",labelsize=20)
plt.rc('figure',titlesize=15)
result = []
for i in range(len(chromosomes)):
	df_chr = df_combined[df_combined["chrom"]==chromosomes[i]]
	df_chr.loc[:,"log2r_smoothed"] = sm.nonparametric.lowess(endog=df_chr["log2r"], 
																	exog = df_chr["start"],
																	return_sorted=False,
																	frac = 6 / len(df_chr.index)  ) # uses ~6 rows, which should be 300kb at 50kb windows
	result += [df_chr]

df_combined_smoothed = pd.concat(result, axis=0)

######################
for i in range(len(chromosomes)):

	f, ax = plt.subplots(figsize=(12,2))
	ax.plot(df_combined_smoothed[df_combined_smoothed["chrom"]==chromosomes[i]]["start"],
			df_combined_smoothed[df_combined_smoothed["chrom"]==chromosomes[i]]["log2r_smoothed"],c="blue",linewidth=1)
	ax.margins(x=0, y=0)
	ax.set_ylim([-4,4])
	ax.set_yticks([])
	ax.set_xticks([])
	ax.set_xlim([0, chromosome_length[chromosomes[i] ] ] )
	ax.axhline(y=0,linestyle="--",c="black")
	ax.ticklabel_format(style="plain")
	### this is the published data
	# ax2=ax.twinx()
	# ax2.plot(df_control[df_control["chrom"]=="X"]["start"],
	# 	df_control[df_control["chrom"]=="X"]["log2r"],c="orange")
	# plt.show()
	# plt.close()

	plt.savefig("/Users/heskett/replication_rnaseq/scripts/4x.repliseq."+chromosomes[i]+".png", 
		dpi=400, transparent=True, bbox_inches='tight', pad_inches = 0)
	plt.close()
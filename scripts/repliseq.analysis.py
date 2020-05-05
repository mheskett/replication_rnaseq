import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
"13","14","15","16","17","18","19","20","21","22","X"]

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
df_combined.loc[:,"percent_early"] = df_combined["early_counts"] / (df_combined["early_counts"] + df_combined["late_counts"])

result = []
for i in range(len(chromosomes)):
	df_chr = df_combined[df_combined["chrom"]==chromosomes[i]]
	df_chr.loc[:,"percent_early_smoothed"] = sm.nonparametric.lowess(endog=df_chr["percent_early"], 
																	exog = df_chr["start"],
																	return_sorted=False,
																	frac = 6 / len(df_chr.index)  ) # uses ~6 rows, which should be 300kb at 50kb windows
	result += [df_chr]

df_combined_smoothed = pd.concat(result, axis=0)

######################
f, ax = plt.subplots()
ax.plot(df_combined_smoothed[df_combined_smoothed["chrom"]=="X"]["start"],
		df_combined_smoothed[df_combined_smoothed["chrom"]=="X"]["percent_early_smoothed"],c="blue")
### this is the published data
ax2=ax.twinx()
ax2.plot(df_control[df_control["chrom"]=="X"]["start"],
	df_control[df_control["chrom"]=="X"]["log2r"],c="orange")
plt.show()
plt.close()
import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
from scipy.stats import norm

from sys import argv
import glob
import statsmodels.api as sm
from sklearn.cluster import KMeans
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/heskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
chromosome_length = {"chr1":249250621,
"chr2":243199373,
"chr3":198022430,
"chr4":191154276,
"chr5":180915260,
"chr6":171115067,
"chr7":159138663,
"chr8":146364022,
"chr9":141213431,
"chr10":135534747,
"chr11":135006516,
"chr12":133851895,
"chr13":115169878,
"chr14":107349540,
"chr15":102531392,
"chr16":90354753,
"chr17":81195210,
"chr18":78077248,
"chr19":59128983,
"chr20":63025520,
"chr21":48129895,
"chr22":51304566,
"chrX":155270560}

def intersect_tables(df1,df2):
    ### return all df1 rows that intersect df2 by >0bs
    ### run time is n_squared ish.....dont use this for big DFs. just for small ones
    ### three possibilities of overlap: 1) start is between start/stop. 2) stop is between start/stop. 3) start <= start AND stop >= stop
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    return result

def add_binom_pval(df):
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
                            row["hap1_counts"]+row["hap2_counts"],
                            p=0.5,
                            alternative="two-sided"), # v slow for some reason 
                            axis=1)
    switcherss = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = switcherss[1]
    df["fdr_reject"] = switcherss[0]
    return

def get_arms(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes)):
        # should be (p end, q end)
        arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
        cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
    return arm_dict
def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["hap1_counts"] >= x["hap2_counts"]:
        return x["hap1_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["hap2_counts"]  / x["total_reads"] + 0.5
    return
def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
		# should be (p end, q end)
		arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
		cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
	return arm_dict
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

def smooth_repli(df):
    ## returns the Y- values after smoothing
    p=[]
    if len(df[df["arm"]=="p"]) <= 10:
        frac_p = 1
    elif len(df[df["arm"]=="p"]) >10:
        frac_p= 6 / len(df[df["arm"]=="p"])

    if len(df[df["arm"]=="p"]) > 0:
        p = sm.nonparametric.lowess(endog=df[df["arm"]=="p"]["logr"], 
                exog=df[df["arm"]=="p"]["start"], 
                return_sorted=False, frac = frac_p )
    ###
    q=[]
    print(len(df[df["arm"]=="q"]) )
    if len(df[df["arm"]=="q"]) <= 10:
        frac_q = 1
    elif len(df[df["arm"]=="q"]) > 10:
        frac_q = 6 / len(df[df["arm"]=="q"])
    if len(df[df["arm"]=="p"]) > 0:
        q = sm.nonparametric.lowess(endog=df[df["arm"]=="q"]["logr"], 
            exog=df[df["arm"]=="q"]["start"], 
            return_sorted=False, frac = frac_q) 
    return p,q

def smooth_vector(x,y):
    y_smooth = []
    if len(x) <= 4:
        frac = 1
    elif len(x) >4:
        frac= 4 / len(x)
    if len(x) > 0:
        y_smooth = sm.nonparametric.lowess(endog=y, 
                exog=x, 
                return_sorted=False, frac = frac )
    return y_smooth

###########################################
###########################################
parser = argparse.ArgumentParser(description="make windows")

parser.add_argument("--maternal_rt",
   type=str,
   metavar="[maternal_rt_file]",
   required=False,
   help="input bed file to plot genome features distributions")
parser.add_argument("--paternal_rt",
   type=str,
   metavar="[paternal_rt_file]",
   required=True,
   help="full path to output results")

arguments = parser.parse_args()

arm_dict = get_arms(cytoband)
all_files_repli = [arguments.maternal_rt, arguments.paternal_rt]
filenames_repli=[os.path.basename(x)[0:9] for x in all_files_repli]
repli_li = []

print(filenames_repli)
maternal_df= pd.read_csv(arguments.maternal_rt,sep="\t",
                        names= ["chrom","start","stop","log2rt"],
                        dtype = {"chrom":str,"start":int,"stop":int,"log2rt":float})
maternal_df["sample"] = os.path.basename(arguments.maternal_rt)[0:9]

paternal_df= pd.read_csv(arguments.paternal_rt,sep="\t",
                        names= ["chrom","start","stop","log2rt"],
                        dtype = {"chrom":str,"start":int,"stop":int,"log2rt":float})
paternal_df["sample"] = os.path.basename(arguments.paternal_rt)[0:9]

## remove Y....
maternal_df = maternal_df[maternal_df["chrom"]!="chrY"]
paternal_df = paternal_df[paternal_df["chrom"]!="chrY"]

## need to functionalize this arm shit, thought i already did...
maternal_df["arm"] = maternal_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
paternal_df["arm"] = paternal_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)



print(maternal_df)
for chrom in chromosomes:
    plt.rc('xtick', labelsize=10) 
    plt.rc('ytick', labelsize=10) 
    f, ax = plt.subplots(1,1,figsize=(16,3))

    maternal = maternal_df[(maternal_df["chrom"]==chrom)]
    paternal = paternal_df[(paternal_df["chrom"]==chrom)]
    print(maternal)
    print(paternal)
    ax.axhline(y=0,lw=0.5,c="black")
    ax.plot(maternal[maternal["arm"]=="p"]["start"],
                smooth_vector(maternal[maternal["arm"]=="p"]["start"],maternal[maternal["arm"]=="p"]["log2rt"]),
            c="blue",lw=1)
    ax.plot(paternal[paternal["arm"]=="p"]["start"],
            smooth_vector(paternal[paternal["arm"]=="p"]["start"],paternal[paternal["arm"]=="p"]["log2rt"]),
            c="red",lw=1) ## -- line style is haplotype 2
   
    ax.plot(maternal[maternal["arm"]=="q"]["start"],
                smooth_vector(maternal[maternal["arm"]=="q"]["start"],maternal[maternal["arm"]=="q"]["log2rt"]),
            c="blue",lw=1)
    ax.plot(paternal[paternal["arm"]=="q"]["start"],
            smooth_vector(paternal[paternal["arm"]=="q"]["start"],paternal[paternal["arm"]=="q"]["log2rt"]),
            c="red",lw=1) ## -- line style is haplotype 2

    ax.set_ylim([-2.2,2.2])
    ax.set_yticks([-2,-1,0,1,2])
    ax.set_xlim([0,chromosome_length[chrom]])
    ax.set_xticks(np.linspace(0,chromosome_length[chrom],10))
    plt.suptitle(arguments.maternal_rt.removesuffix(".bedgraph")+"_"+arguments.paternal_rt.removesuffix(".bedgraph")+"_"+chrom)

    plt.savefig(arguments.maternal_rt.removesuffix(".bedgraph")+"_"+arguments.paternal_rt.removesuffix(".bedgraph")+"_"+chrom+".png",
        dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)    




exit()

repli_li.append(df_repli)
repli_df = pd.concat(repli_li)

paternal=[x for x in filenames_repli if "_p" in x]
maternal=[x for x in filenames_repli if "_m" in x]
paternal.sort()
maternal.sort()
print(paternal,maternal)




exit()
# repli_df.loc[:,"logr_hap1"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
# repli_df.loc[:,"logr_hap2"] = repli_df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
# repli_df.loc[:,"logr_diff_abs"] = abs(repli_df["logr_hap1"] - repli_df["logr_hap2"]) ## 
# repli_df.loc[:,"logr_diff_raw"] = repli_df["logr_hap1"] - repli_df["logr_hap2"] # positive if hap1 early, negative if hap2 early
# repli_df.loc[:,"logr"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+x["hap2_early"]+1) / (x["hap1_late"]+x["hap2_late"]+1)), axis=1 )
# sig_repli = np.percentile(a = repli_df[repli_df["chrom"]!="X"]["logr_diff_abs"], q = 95)
# repli_df["sig_repli"]=["True" if x > sig_repli else "False" for x in repli_df["logr_diff_raw"]]
# repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
# color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "Blue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
# repli_df["repli_color"] = color_vector
# asynchronous_regions = repli_df[repli_df["sig_repli"]=="True"]
#####
zscore = lambda x: (x - x.mean()) / x.std()

# for i in repli_df["samples"].unique():
#     mean=repli_df[(repli_df["chrom"]!="X")&(repli_df["sample"]==i)]["logr_diff_abs"].mean()
#     std=repli_df[(repli_df["chrom"]!="X")&(repli_df["sample"]==i)]["logr_diff_abs"].std()

### this is figure 4a

repli_df = repli_df[repli_df["chrom"]!="X"]
repli_df["logr_diff_abs_sample_zscore"] = repli_df.groupby("sample")["logr_diff_abs"].transform(zscore)
print(repli_df)
samples=repli_df["sample"].unique()

for i in range(len(samples)):
	for j in range(len(autosomes)):
		f,ax = plt.subplots(1,1,figsize=(10,2))
		df_chrom_p = repli_df[(repli_df["sample"]==samples[i])&(repli_df["chrom"]==autosomes[j])&(repli_df["arm"]=="p")]
		df_chrom_q = repli_df[(repli_df["sample"]==samples[i])&(repli_df["chrom"]==autosomes[j])&(repli_df["arm"]=="q")]
		print(df_chrom_p)
		print(df_chrom_q)
		max_val = max(repli_df[(repli_df["sample"]==samples[i])&(repli_df["chrom"]==autosomes[j])]["logr_diff_abs_sample_zscore"])
		min_val = min(repli_df[(repli_df["sample"]==samples[i])&(repli_df["chrom"]==autosomes[j])]["logr_diff_abs_sample_zscore"])
		ax.plot(df_chrom_p["start"],df_chrom_p["logr_diff_abs_sample_zscore"],c="black")
		ax.plot(df_chrom_q["start"],df_chrom_q["logr_diff_abs_sample_zscore"],c="black")

		for index,row in df_chrom_p[(df_chrom_p["logr_diff_abs_sample_zscore"]>=2.5)].iterrows():
			rect=Rectangle((row["start"],min_val),width=row["stop"]-row["start"],height=max_val+abs(min_val),
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax.add_patch(rect)		
		for index,row in df_chrom_q[(df_chrom_q["logr_diff_abs_sample_zscore"]>=2.5)].iterrows():
			rect=Rectangle((row["start"],min_val),width=row["stop"]-row["start"],height=max_val+abs(min_val),
				facecolor=row["repli_color"],alpha=0.6,fill="True")
			ax.add_patch(rect)
		ax.set_ylim([min_val,max_val])
		ax.set_xlim([0,chromosome_length[autosomes[j]]])
		ax.set_xticks(np.linspace(0,chromosome_length[autosomes[j]],16))
		# plt.show()
		plt.savefig(samples[i]+"as.repli."+str(autosomes[j])+".png",dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)

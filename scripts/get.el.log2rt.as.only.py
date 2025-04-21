import os
import re
import csv
import numpy as np
import pandas as pd
import pybedtools
import scipy.stats
import seaborn as sns
from scipy.stats import norm
from sys import argv
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


def library_size_normalization(df):
    new_df = df.copy()
    total_reads = new_df["hap1_counts"].sum() + new_df["hap2_counts"].sum()
    new_df["hap1_lsm"] = (new_df["hap1_counts"] / total_reads) * 10**6
    new_df["hap2_lsm"] = (new_df["hap2_counts"] / total_reads) * 10**6

    return new_df

######
### paternal is LEFT column aka column 4, maternal is right column aka column 5
### this is for paternal or maternal files which have one column
### have to figure out the system for hap1 vs hap2 when these files are going to be ALL hap1 or ALL hap2

if "_p_" in argv[1]:
    prefix="paternal"
elif "_m_" in argv[1]:
    prefix="maternal"
print("file detected: "+argv[1][0:12]+" "+prefix)


early_df= pd.read_csv(argv[1],sep="\t",
                        names= ["chrom","start","stop","paternal_coverage_early","maternal_coverage_early"],
                        dtype = {"chrom":str,"start":int,"stop":int}).set_index(["chrom","start","stop"])
#early_df["sample"] = os.path.basename(argv[1])[0:11]

late_df= pd.read_csv(argv[2],sep="\t",
                        names= ["chrom","start","stop","paternal_coverage_late","maternal_coverage_late"],
                        dtype = {"chrom":str,"start":int,"stop":int}).set_index(["chrom","start","stop"])
#late_df["sample"] = os.path.basename(argv[2])[0:11]

## hard coding naming convention
df=pd.concat([early_df,late_df],axis=1).reset_index().dropna(how="any", axis="index")

early_read_sum = df["paternal_coverage_early"].sum(axis="rows") + df["maternal_coverage_early"].sum(axis="rows")
df["paternal_early_lsm"] = df["paternal_coverage_early"] /  early_read_sum * 10**6
df["maternal_early_lsm"] = df["maternal_coverage_early"] / early_read_sum * 10**6

late_read_sum = df["paternal_coverage_late"].sum(axis="rows") + df["maternal_coverage_late"].sum(axis="rows")
df["paternal_late_lsm"] = df["paternal_coverage_late"] / late_read_sum * 10**6
df["maternal_late_lsm"] = df["maternal_coverage_late"] / late_read_sum * 10**6


df["log2rt"] = np.log2((df[prefix+"_early_lsm"]+1) / (df[prefix+"_late_lsm"]+1))

df=df[df["log2rt"]!=0]
df[["chrom","start","stop","log2rt"]].to_csv(argv[1][0:10]+"el.log2rt.as.only.bed",sep="\t",header=False,index=False)





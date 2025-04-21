import csv
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
# pybedtools.set_bedtools_path("/Users/michaelheskett/miniconda3/envs/for_bedtools/bin/bedtools/")
import matplotlib.patheffects as path_effects
import scipy.stats
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from scipy.signal import find_peaks
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Shadow
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
from scipy.stats import ttest_ind
import glob

def add_binom_pval(df):
    """
    example:
    >>> scipy.stats.binomtest(5,10,p=0.5)
    BinomTestResult(k=5, n=10, alternative='two-sided', statistic=0.5, pvalue=1.0)
    """
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binomtest(row["paternal_counts"],
                                row["paternal_counts"]+row["maternal_counts"],
                            p=0.5,
                            alternative="two-sided").pvalue, # v slow for some reason 
                            axis=1)
    results = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = results[1]
    df["fdr_reject"] = results[0]


def helper_func(x):
    if x["total_reads"]==0: # try this for filtering
        return 0
    elif x["paternal_counts"] >= x["maternal_counts"]:
        return x["paternal_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["maternal_counts"]  / x["total_reads"] + 0.5
    return

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
def get_arms(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes)):
        # should be (p end, q end)
        arm_dict[chromosomes[i]] = (cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("p"))]["stop"].max(),
        cytoband[(cytoband["chrom"]==chromosomes[i]) & (cytoband["arm"].str.contains("q"))]["stop"].max())
    return arm_dict

def get_arms_nochr(cytoband):
    ## given a data frame with genome elements, add the arm information to a new column
    arm_dict = {}
    for i in range(len(chromosomes_nochr)):
        # should be (p end, q end)
        arm_dict[chromosomes_nochr[i]] = (cytoband_nochr[(cytoband_nochr["chrom"]==chromosomes_nochr[i]) & (cytoband_nochr["arm"].str.contains("p"))]["stop"].max(),
        cytoband_nochr[(cytoband_nochr["chrom"]==chromosomes_nochr[i]) & (cytoband_nochr["arm"].str.contains("q"))]["stop"].max())
    return arm_dict


def sum_region_length(df):
    diffs = df["stop"] - df["start"]
    return diffs.sum()

### REMOVIONG 15 for BAD DATA!!!!
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X"]
autosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

chromosomes_nochr = ["1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X"]

arms = ["p","q"]
#### for arm level data to skip over centromeres                
cytoband = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.chr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])

cytoband_nochr = pd.read_table("/Users/michaelheskett/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
                            names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)
arm_dict_nochr=get_arms_nochr(cytoband_nochr)

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

chromosome_length_nochr = {"1":249250621,
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


chromosome_length_mouse = {"chr1":    195471971,
"chr2":    182113224,
"chrX":    171031299,
"chr3":    160039680,
"chr4":    156508116,
"chr5":    151834684,
"chr6":    149736546,
"chr7":   145441459,
"chr10":   130694993,
"chr8":    129401213,
"chr14":   124902244,
"chr9":   124595110,
"chr11":   122082543,
"chr13":   120421639,
"chr12":   120129022,
"chr15":   104043685,
"chr16":   98207768,
"chr17":   94987271,
"chrY" :   91744698,
"chr18":   90702639,
"chr19":   61431566
}

chromosomes_mouse = {"chr1",
"chr2",
"chrX",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr10",
"chr8",
"chr14",
"chr9",
"chr11",
"chr13",
"chr12",
"chr15",
"chr16",
"chr17",
"chrY",
"chr18",
"chr19"}


### peakfinder for RT STD DEV data 

### use segment.repliseq.py on RT STD DEV track to find peaks. 


df_repli = pd.read_csv("acp6.rt.txt",sep="\t",header=0,index_col=False,
	dtype={"chrom":str,"start":int,"stop":int})
regions=[]
df_repli = df_repli.sort_values(["chrom","start","stop"])
print(df_repli)
for chrom in chromosomes_nochr:
		
	for arm in arms:

		tmp = df_repli[(df_repli["chrom"]==chrom) & (df_repli["arm"]==arm)].reset_index(drop=True)
		print("arm", arm)
		print("tmp", tmp)
		if len(tmp)<=10:
		    continue

		tmp_smoothed = lowess(endog=tmp["acp6_std_dev_both_haps"],
		                     exog=tmp["start"],
		                      frac=10/len(tmp["start"]),return_sorted=False)
		#### fix these hard coded filters for the peak finding to deal with STD DEV
		#### 
		####
		peaks,properties = find_peaks(tmp_smoothed , width=0, height=0, distance=1, prominence=.1, threshold=0)
		# valleys,valley_properties = find_peaks(-tmp_smoothed,width=0, height=-20, distance=1, prominence=0.6, threshold=0)
		#####
		# print("valleys valley properties",valleys,valley_properties)

		for i in range(len(properties["left_ips"])):

			regions+=[ [ chrom, 
						tmp.loc[properties["left_ips"][i].round(),:]["start"], 
						tmp.loc[properties["right_ips"][i].round(),:]["start"], 
						arm, 
						properties["prominences"][i], 
						tmp.loc[properties["right_ips"][i].round(),:]["start"] - tmp.loc[properties["left_ips"][i].round(),:]["start"],  
						properties["peak_heights"][i], 
						"peak" ]   ]

	### warning. height parameter set to 0 seems to actually filter based on height somehow. bug? 
	### oh because youre using negative numbers it might fuck this up. find way to shift the curves up
	## so min is set to 0. can add the abs value of the min to each curve to shift to positive realm
	## but i think youre looking for prominence, which should be in an absolute unit
		print("regions",regions)
		f,ax = plt.subplots(figsize=(12,2))
		ax.plot(tmp["start"],tmp_smoothed,linewidth=2)

		ax.plot(tmp.loc[peaks,:]["start"],
		        tmp_smoothed[peaks], "x",markersize=5)

		# ax.plot(tmp.loc[valleys,:]["start"],
		#         tmp_smoothed[valleys], "x",c="green",markersize=5)

		plt.ylim([0,1.25])

		###
		for i in range(len(properties["width_heights"])):
		    rect = Rectangle(xy = (tmp.loc[properties["left_ips"][i].round(),:]["start"]   , -10),
		                        width= tmp.loc[properties["right_ips"][i].round(),:]["start"] - tmp.loc[properties["left_ips"][i].round(),:]["start"],
		                        height=20, fill=True, alpha=0.5,facecolor="gray")
		    ax.add_patch(rect)

		plt.suptitle(chrom+"_"+arm)
		plt.margins(0,0)
		# plt.show()
		plt.savefig("acp6."+chrom+arm+".peaks.pdf")
		plt.close()
result = pd.DataFrame(regions,columns=["chrom","start","stop","arm","peak_prominence","peak_width","peak_height","peak_valley"])

	# result.to_csv(arguments.repli_bed.rstrip(".bed")+"_peaks.bed",header=None,index=None,sep="\t")
	# result.to_csv(arguments.repli_bed.rstrip(".bed")+"_peaks.txt",header=True,index=None,sep="\t")


plt.figure()
plt.hist(result["peak_prominence"],bins=30)
# plt.xticks(list(range(0,16)))
plt.suptitle("peak prominence")
# plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_prominence.pdf")
plt.show()

plt.figure()
plt.hist(result["peak_width"],bins=30)
# plt.xticks(list(range(0,16)))
plt.suptitle("peak width")
# plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_width.pdf")
plt.show()

plt.figure()
plt.hist(result["peak_height"],bins=30)
# plt.xticks(list(range(0,16)))
plt.suptitle("peak height")
# plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_height.pdf")
plt.show()

print(result)

plt.scatter(result["peak_height"],result["peak_width"],s=15,edgecolor="black",lw=0.2)
plt.show()






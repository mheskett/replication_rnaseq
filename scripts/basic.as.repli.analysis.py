import csv
import os
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import matplotlib.patheffects as path_effects
import scipy.stats
from scipy.stats import norm
import statsmodels.api as sm
import statsmodels.stats.multitest as mt
import pickle
import statsmodels.api as sm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Shadow

def add_binom_pval(df):
    ## binomial p value is a good starting point, but not as stringent
    ## as finding outliers in the distribution of number of biased reads per 
    ## window. this function directly modifies the input data frame
    df["binom_pval"] = df.apply(lambda row: scipy.stats.binom_test(row["hap1_counts"],
                            row["hap1_counts"]+row["hap2_counts"],
                            p=0.5,
                            alternative="two-sided"), # v slow for some reason 
                            axis=1)
    switchers = mt.multipletests(pvals=df["binom_pval"], 
                                alpha=0.01,
                                method="fdr_bh")
    df["fdr_pval"] = switchers[1]
    df["fdr_reject"] = switchers[0]
    return
def calculate_skew(x):
    if x["total_reads"]==0:
        return 0
    elif x["hap1_counts"] >= x["hap2_counts"]:
        return x["hap1_counts"]  / x["total_reads"] - 0.5
    else:
        return -x["hap2_counts"]  / x["total_reads"] + 0.5
    return
def intersect_tables(df1,df2):
    ### return all df1 rows that intersect df2 by >0bs
    ### this is only for smaller data tables. for very large files, run separately. 
    a = pybedtools.BedTool.from_dataframe(df1)
    b = pybedtools.BedTool.from_dataframe(df2)
    ## for plotting add slop to the AS-RT region which is b. a is the lncrnas
    b = pybedtools.BedTool.slop(b,b=250000,g="human_g1k_v37.fasta.fai")
    result = a.intersect(b,wa=True,wb=True).to_dataframe(names=list(df1.columns) + [x+'1' for x in df2.columns])
    result["chrom"] = result["chrom"].astype(str)
    return result
def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column.
    ## this is required for plotting whole chromosome views
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
    necessary for comparing samples vs sample, but not necessary for comparing
    within one sample
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return (df_qn)

def smooth_repli(df):
    ## returns the Y- values after smoothing using lowess
    ## must smooth each chromosome arm separately
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
    ## smooth smaller fragments of data for plotting purposes
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
#################################################################
# Necessary variables for plotting and etc. These can be loaded from a
# separate specifications file, but are kept in-text here for easier
# development. 
chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
				"13","14","15","16","17","18","19","20","21","22","X"]
arms = ["p","q"]
cytoband = pd.read_table("/Users/mike/replication_rnaseq/scripts/cytoband.nochr.hg19.bed",sep="\t",
							names =["chrom","start","stop","arm","band"])
arm_dict = get_arms(cytoband)
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
############################
# load and process repli-seq files
############## 
all_files_repli = [""]

filenames_repli=[os.path.basename(x)[0:15] for x in all_files_repli]
repli_li = []
for i in range(len(all_files_repli)):
    df_repli = pd.read_csv(all_files_repli[i],sep="\t",
                        names= ["chrom","start","stop","hap1_early","hap2_early","hap1_late","hap2_late"],
                        dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts_plus":str,"hap2_counts_plus":str,"hap1_counts_minus":str,"hap2_counts_minus":str,
                        "hap1_early":str,"hap2_early":str,"hap1_late":str,"hap2_late":str})

    tmp = df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]].replace(".",0)
    tmp = tmp.astype(int)
    df_repli.loc[:,["hap1_early","hap2_early","hap1_late","hap2_late"]] = tmp
    df_repli = df_repli.set_index(["chrom","start","stop"])
    df_repli = df_repli[df_repli.sum(axis="columns")!=0]
    df_repli = df_repli.reset_index()
    df_repli["sample"] = filenames_repli[i]
    repli_li.append(df_repli)
repli_df = pd.concat(repli_li)

repli_df.loc[:,"logr_hap1"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+1) / (x["hap1_late"]+1)), axis=1 )
repli_df.loc[:,"logr_hap2"] = repli_df.apply(lambda x: np.log2((x["hap2_early"]+1) / (x["hap2_late"]+1)), axis=1 )
repli_df.loc[:,"logr_diff_abs"] = abs(repli_df["logr_hap1"] - repli_df["logr_hap2"]) ## 
repli_df.loc[:,"logr_diff_raw"] = repli_df["logr_hap1"] - repli_df["logr_hap2"] # positive if hap1 early, negative if hap2 early
repli_df.loc[:,"logr"] = repli_df.apply(lambda x: np.log2((x["hap1_early"]+x["hap2_early"]+1) / (x["hap1_late"]+x["hap2_late"]+1)), axis=1 )
repli_df["arm"] = repli_df.apply(lambda x: "q" if (x["stop"] > arm_dict[x["chrom"]][0]) & (x["stop"] <= arm_dict[x["chrom"]][1]) else "p", axis=1)
color_vector = ["Red" if (row["logr_hap1"] >= row["logr_hap2"]) else "mediumblue" for index,row in repli_df.iterrows() ] # red if hap1 early, blue if hap2 early
repli_df["repli_color"] = color_vector
#####
zscore = lambda x: (x - x.mean()) / x.std()
repli_df = repli_df[repli_df["chrom"]!="X"]
repli_df["logr_diff_abs_sample_zscore"] = repli_df.groupby("sample")["logr_diff_abs"].transform(zscore)
###################################
### now get all the TLs an DAE TLs in EB2 and EB10
vlinc_files=[""]
dfs = []
for i in range(len(vlinc_files)):
    df = pd.read_csv(vlinc_files[i],sep="\t",
                            names= ["chrom","start","stop","name","rpkm","strand","l1_fraction","hap1_counts","hap2_counts"],
                            dtype = {"chrom":str,"start":int,"stop":int,"hap1_counts":int,"hap2_counts":int})
    df["total_reads"] = df["hap1_counts"] + df["hap2_counts"]
    df["skew"] = df.apply(helper_func, axis = 1)
    df["sample"] = os.path.basename(vlinc_files[i])[0:9]
    add_binom_pval(df)
    dfs += [df]
df = pd.concat(dfs)
unique_genes = list(df["name"].drop_duplicates())
model = pickle.load(open("eb.variance.coding.model.sav", 'rb'))
df["significant_deviation"] = df.apply(lambda x: True if abs(x["hap1_counts"] - x["total_reads"]/2) >= model.predict(np.array([x["total_reads"]]).reshape(1,-1))*2.5 else False,
    axis=1)
df_significant_rows = df[df["significant_deviation"]==True]
df_nonsignificant_rows = df[df["significant_deviation"]==False]
df=df[df["total_reads"]>=10]
color_dict = {"bouha.4.a":"r","bouha.15.":"c","bouha.10.":"orange","bouha.3.a":"g",
"bouha.2.a":"b","bouha.13.":"green"}
df["color"]=[color_dict[x] for x in df["sample"]]
########
########
########
## dual plot of ASRT regions and DAE lncrna
## DAE and DART are independent normally distributed random variables (different distribution for each expression level) 
## just used simple outlier detection strategy for both 1-dimensional variables.
## a rigorous statistical procedure to calculate p-values involves consideration of
## the allele-specific read counts of ALL heterozygous SNPs within a lncRNA or
## replication timing window. because we have fully haplotype resolved genomes,
## that is, 1 haplotype block per chromosome, statistical power is increased
## by enumerating haplotype specific reads across large regions. 
## the null hypothesis for DAE lncRNA was generated using the expected variation
## of allele specific reads across ALL expressed regions at various expression levels.
## total read counts in repli-seq do not vary significantly across the genome
## analogously to expression levels.
## a normal distribution of DART was generated, and then right-tailed outliers were considered as
## DART. in future experiments where large sample numbers may be involved, 
## concordance or discordance across samples can be used as additional criteria
## to identify DART regions--the ones presented here are only preliminary and must 
## be validated. the goal here is to identify the strongest examples of DART within the low sample
## size set that we currently have. 
dae_dart = intersect_tables(df[(df["significant_deviation"]==True) & (df["chrom"]!="X")],
    repli_df[repli_df["logr_diff_abs_sample_zscore"]>=2])
sample_pairs = [("bouha.2.a","bouha.2.repli.5"),("bouha.10.","bouha.10.repli.")]
## hap1 lnc, hap1 early, hap1 lnc hap2 early, hap2 lnc hap1 early, hap2 lnc  hap2 early
counter = [] 
for i in range(len(sample_pairs)): 
    lnc_as = dae_dart[(dae_dart["sample"]==sample_pairs[i][0]) & (dae_dart["sample1"]==sample_pairs[i][1])] ## sample specific
    for index,row in lnc_as.drop_duplicates(["chrom","start","stop"]).iterrows():
        counts = [0,0]
        if row["skew"]>0:
            counts[0]=1
        if row["logr_diff_raw1"] >0:
            counts[1]=1
        counter += [counts]
        f,ax = plt.subplots(1,1,figsize=(2,2),sharex=False)
        start=row["start"]
        stop=row["stop"]
        chrom=str(row["chrom"])
        plt.suptitle(chrom+":"+sample_pairs[i][0])
        ax.tick_params(axis='x', labelsize= 4)
        ax.set_xlim([max(0,start-2000000),stop+2000000])
        ax.set_xticks(np.linspace(max(0,start-2000000),stop+2000000, 8))

        repli_tmp =repli_df[(repli_df["chrom"]==chrom)
            &(repli_df["sample"]==sample_pairs[i][1]) 
            & (repli_df["start"].between(start-2500000,stop+2500000))]
        ax.plot(repli_tmp["start"],
                smooth_vector(repli_tmp["start"],repli_tmp["logr_hap1"]),
                lw=1,zorder=1,c=color_dict[sample_pairs[i][0]])
        ax.plot(repli_tmp["start"],
            smooth_vector(repli_tmp["start"],repli_tmp["logr_hap2"]),
            linestyle="--",lw=1,zorder=2,c=color_dict[sample_pairs[i][0]])

        for index3, row3 in repli_tmp[repli_tmp["logr_diff_abs_sample_zscore"]>=2].iterrows():
            rect=Rectangle((row3["start"], -10), width=row3["stop"]-row3["start"], height=20,
                         facecolor="lightgray",alpha=1,fill=True)   
            ax.add_patch(rect)

        ax_lnc=ax.twinx()
        ax_lnc.set_ylim([-0.6,0.6])
        ax_lnc.set_yticks([-.5,-.25,0.,.25,.5])
        for index2, row2 in lnc_as[(lnc_as["chrom"]==chrom) & (lnc_as["start"]==start) & (lnc_as["stop"]==stop) ].iterrows():
            rect=Rectangle((row2["start"], row2["skew"]-.05), width=row2["stop"]-row2["start"], height=0.05,
                     facecolor=row2["color"], edgecolor="black",fill=True,lw=.5)
            shadow = Shadow(rect, 10000,-0.0015 )                                        
            ax_lnc.add_patch(shadow)
            ax_lnc.add_patch(rect)
        plt.rc('xtick', labelsize=4)
        plt.rc('ytick', labelsize=8) 
        ##########0 
        ax.set_ylim([min(min(repli_tmp["logr_hap1"]),min(repli_tmp["logr_hap2"])),
        max(max(repli_tmp["logr_hap1"]),max(repli_tmp["logr_hap2"]))])
        plt.savefig("vlinc.logr.diff."+sample_pairs[i][0]+str(chrom)+"."+str(start)+"."+str(stop)+".png",
            dpi=400,transparent=True, bbox_inches='tight', pad_inches = 0)

        plt.close()
    print("sample: ",sample_pairs[i][0])
    plt.show()
test = pd.DataFrame(counter,columns=["hap1_skew","hap1_early"]).groupby(["hap1_skew","hap1_early"]).size()
test.plot.pie(figsize=(2,2),fontsize=20)
plt.show()
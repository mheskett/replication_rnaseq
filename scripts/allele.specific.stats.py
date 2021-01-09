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
from sys import argv

df_plus = pd.read_csv(argv[1],sep="\t",names=["chrom","start","stop","hap1","hap2","count1","count2"],
	dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"count1":int,"count2":int})

df_minus = pd.read_csv(argv[2],sep="\t",names=["chrom","start","stop","hap1","hap2","count1","count2"],
dtype={"chrom":str,"start":int,"stop":int,"hap1":str,"hap2":str,"count1":int,"count2":int})

df_plus["strand"] = "+"
df_minus["strand"] = "-"

df = pd.concat([df_plus,df_minus])
print("non-zero read snps: ",len(df[df["count1"]+df["count2"]>0]))
df1 = df[df["count1"]+df["count2"]>0]
df1["sum"] = df1["count1"]+df1["count2"]
stats = (df1["count1"]+df1["count2"]).describe(percentiles=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,])
print("num non zero SNPs: ",len(df1))
plt.hist(df1[df1["sum"]<=50]["sum"],bins=20)
plt.figtext(0.5,0.2,stats.to_string()
	)
plt.suptitle(argv[1]+"\n"+"num covered het snps: "+str(len(df1)))
plt.savefig(argv[1]+".png")
plt.close()

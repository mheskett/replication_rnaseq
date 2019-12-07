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


df = pd.read_csv("/Users/heskett/replication_rnaseq/eqtl.data/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.egenes.significant.bed",sep="\t",index_col=None,header=None,
	names=["gene_chr","gene_start","gene_end","gene_name","num_var","strand","beta_shape1","beta_shape2","true_df", "pval_true_df","variant_id",
	"tss_distance","chr","pos","ref","alt","num_alt_per_site","rs_id_dbSNP147_GRCh37p13","minor_allele_samples","minor_allele_count","maf","ref_factor",
	 "pval_nominal","slope","slope_se","pval_perm", "pval_beta", "qval", "pval_nominal_threshold", "log2_aFC"])

print(df)
print(plt.hist(df["log2_aFC"],bins=50))
plt.show()
plt.close()
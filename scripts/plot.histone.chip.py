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
import scipy.stats
import seaborn as sns
import statsmodels.api as sm


df = pd.read_table("H3K4me3-human.all.nochr.allele.counts.haplotype.resolved.counts.bed",
	sep="\t",
	names=["chrom","start","stop","hap1","hap2","hap1_counts","hap2_counts"])






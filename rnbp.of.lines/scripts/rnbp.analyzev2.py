import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools
import scipy.stats
import seaborn as sns
import multiprocessing as mp


## for local use
eclip_file = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.hepg2.chr6.nochr.sorted.bed"
introns_file =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs_file = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines_file = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

################
eclip = pybedtools.BedTool(eclip_file) # change this to bed file of previously determined windows of interest
introns = pybedtools.BedTool(introns_file)
vlincs = pybedtools.BedTool(vlincs_file)
lines = pybedtools.BedTool(lines_file)

################
df_eclip = pd.read_table(eclip_file,
					names=["chrom", "start", "stop", "name","score", 
					"strand", "signalvalue","pvalue","qvalue", "peak"],
					 dtype={"chrom":str, "start":int, "stop":int,
					 "name":str, "score":float, "strand":str,
					 "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})
df_eclip = df_eclip[ (~df_eclip.name.str.contains("IDR",case=False)) ] # filter out replication experiments
eclip_experiments = list( df_eclip.name.unique() )

##############

def introns_overlap_rnbp(introns_file, rnbp_name , eclip_df):
	## should return a data frame of all introns and the counts of their overlaps with a RNBP
	eclip_object = eclip_df[ eclip_df["name"] == rnbp_name ].from_dataframe()
	introns.coverage(eclip_object, F=0.9).to_dataframe(names=["chrom", "start", "stop", "intron_name",
														"score", "strand", "peak_count", "num_bases_covered",
														"intron_length", "fraction_peak"],
														dtype={"chrom":str,"start":int,"stop":int,
														"intron_name":str, "score":str,"strand":int,
														"peak_count":int,"num_bases_covered":int, "intron_length":int, "fraction_peak":float})
	return 
def sample_introns(length, introns_file, rnbp_name):
	## intersect introns file with RNBP. This should only be done once per RNBP.
	## returns 
	return


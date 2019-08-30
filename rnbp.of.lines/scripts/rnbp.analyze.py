import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools


eclip = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.chr6.nochr.sorted.bed"
introns =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

a = pybedtools.BedTool(eclip) # change this to bed file of previously determined windows of interest
b = pybedtools.BedTool(introns)
c = pybedtools.BedTool(vlincs)
d = pybedtools.BedTool(lines)

lines_in_introns = d.intersect(b) # wb =true to know which intron a line is in
lines_in_vlincs = d.intersect(c) # wb = true to know which vlinc a line is in

#eclip_binding_lines_in_vlincs = a.coverage(lines_in_vlincs,counts=True)
eclip_binding_lines_in_introns = a.intersect(lines_in_introns).to_dataframe(names=["chrom", "start", "stop", "name", 
																		"score", "strand","signalvalue","pvalue","qvalue","peak"],
																		dtype={"chrom":str, "start":int, "stop":int,
																		"name":str, "score":float, "strand":str,
																		"signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})
print(eclip_binding_lines_in_introns)
### may need to flip this and do more manual counting.
eclip_binding_lines_in_vlincs = a.intersect(lines_in_vlincs,c=True).to_dataframe(names=["chrom", "start", "stop", "name", 
																		"score", "strand","signalvalue","pvalue","qvalue","peak","count"],
																		dtype={"chrom":str, "start":int, "stop":int,
																		"name":str, "score":float, "strand":str,
																		"signalvalue":float,"pvalue":float,"qvalue":int,"peak":int,"count":int})

## goal is to get binding sites per LINE for each protein for set of intronic lines and vlinc lines
## 
# print(eclip_binding_lines_in_vlincs.groupby)
# eclip_binding_lines_in_introns = eclip_binding_lines_in_introns[eclip_binding_lines_in_introns["count"]>0]
# eclip_binding_lines_in_vlincs = eclip_binding_lines_in_vlincs[eclip_binding_lines_in_vlincs["count"]>0]

# print(eclip_binding_lines_in_introns.groupby(["name"])["count"].mean())

# df = pd.read_table("/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.chr6.nochr.sorted.bed",
# 	sep="\t",names=["chrom","start","stop","name","score","strand","signal_value","pvalue","qvalue","peak"],
# 	index_col=None)

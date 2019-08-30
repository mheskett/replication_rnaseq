import os
import csv
import numpy as np
import pandas as pd
import argparse
import re
import matplotlib.pyplot as plt
import pybedtools

# for the cluster
 eclip = "/home/groups/Spellmandata/heskett/replication.rnaseq/rnbp.of.lines/big.downloadall.eclip.nochr.sorted.bed"
introns =  "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines = "/home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

## for local use
eclip = "/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.chr6.nochr.sorted.bed"
introns =  "/Users/heskett/replication_rnaseq/annotation.files/ucsc.introns.filtered.hg19.bed"
vlincs = "/Users/heskett/replication_rnaseq/annotation.files/mcaron.vlinc1541.vlinc2149.merged.final.hg19g1k.bed"
lines = "/Users/heskett/replication_rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed"

a = pybedtools.BedTool(eclip) # change this to bed file of previously determined windows of interest
b = pybedtools.BedTool(introns)
c = pybedtools.BedTool(vlincs)
d = pybedtools.BedTool(lines)

lines_in_introns = d.intersect(b,f=0.9) # wb =true to know which intron a line is in
lines_in_vlincs = d.intersect(c,f=0.9) # wb = true to know which vlinc a line is in

total_intronic_l1s = len(lines_in_introns)
total_vlinc_l1s = len(lines_in_vlincs)

#eclip_binding_lines_in_vlincs = a.coverage(lines_in_vlincs,counts=True)
eclip_binding_lines_in_introns = lines_in_introns.intersect(a,wb=True,F=0.9).to_dataframe(names=["chrom", "start", "stop", "name", 
																		 "score", "strand","chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip",
																		 "strand_eclip","signalvalue","pvalue","qvalue","peak"],
																		 dtype={"chrom":str, "start":int, "stop":int,
																		 "name":str, "score":float, "strand":str,
																		 "chrom_eclip":str,"start_eclip":int,"stop_eclip":int,"name_eclip":str,"score_eclip":int,
																		 "strand_eclip":str, "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})


eclip_binding_lines_in_vlincs = lines_in_vlincs.intersect(a,wb=True,F=0.9).to_dataframe(names=["chrom", "start", "stop", "name", 
																		 "score", "strand","chrom_eclip","start_eclip","stop_eclip","name_eclip","score_eclip",
																		 "strand_eclip","signalvalue","pvalue","qvalue","peak"],
																		 dtype={"chrom":str, "start":int, "stop":int,
																		 "name":str, "score":float, "strand":str,
																		 "chrom_eclip":str,"start_eclip":int,"stop_eclip":int,"name_eclip":str,"score_eclip":int,
																		 "strand_eclip":str, "signalvalue":float,"pvalue":float,"qvalue":int,"peak":int})

eclip_binding_lines_in_vlincs = eclip_binding_lines_in_vlincs[(~eclip_binding_lines_in_vlincs.name_eclip.str.contains("IDR",case=False)) 
																& (eclip_binding_lines_in_vlincs["pvalue"]>=2) 
																& (eclip_binding_lines_in_vlincs["signalvalue"]>=1)]

eclip_binding_lines_in_introns = eclip_binding_lines_in_introns[(~eclip_binding_lines_in_introns.name_eclip.str.contains("IDR",case=False)) 
																& (eclip_binding_lines_in_introns["pvalue"]>=2) 
																& (eclip_binding_lines_in_introns["signalvalue"]>=1)]



fraction_of_intronic_lines_bound = eclip_binding_lines_in_introns.groupby(["name_eclip"]).size() / total_intronic_l1s
fraction_of_vlinc_lines_bound = eclip_binding_lines_in_vlincs.groupby(["name_eclip"]).size() / total_vlinc_l1s




fraction_of_vlinc_lines_bound.divide(fraction_of_intronic_lines_bound).to_csv("vlinc_vs_intron_l1_rnbp_binding.txt",sep="\t")





## goal is to get binding sites per LINE for each protein for set of intronic lines and vlinc lines
## 
# print(eclip_binding_lines_in_vlincs.groupby)
# eclip_binding_lines_in_introns = eclip_binding_lines_in_introns[eclip_binding_lines_in_introns["count"]>0]
# eclip_binding_lines_in_vlincs = eclip_binding_lines_in_vlincs[eclip_binding_lines_in_vlincs["count"]>0]

# print(eclip_binding_lines_in_introns.groupby(["name"])["count"].mean())

# df = pd.read_table("/Users/heskett/replication_rnaseq/rnbp.of.lines/data/all.eclip.chr6.nochr.sorted.bed",
# 	sep="\t",names=["chrom","start","stop","name","score","strand","signal_value","pvalue","qvalue","peak"],
# 	index_col=None)

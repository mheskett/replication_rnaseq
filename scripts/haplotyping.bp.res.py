import os
import re
import csv
import numpy as np
import pandas as pd
import argparse
import re
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generates batch scripts for mutation viewing")

    parser.add_argument("--bed",
       type=str,
       metavar="[NA12878 bed file]",
       required=True,
       help="")
    parser.add_argument("--out_directory", #differentiate between out directory and snapshot directory
       type=str,
       metavar="[out directory]",
       required=True,
       help="full path to output results")


    arguments = parser.parse_args()

    with open(arguments.bed) as f:
        lines = f.readlines()
        data = [x.rstrip("\n").split("\t") for x in lines]
        f.close()
## format: 1       740737  740738  C       T,<NON_REF>     6,2,0   1       740737  740738  T       C
## where haplotyped platinum genome is on the right side
    result = []
    # Formatting and filtering
    for i in range(len(data)):
      alt_genotype = data[i][4].split(",")
      counts = list(map(int,data[i][5].split(",")))
      ref_genotype = data[i][3]

      if sum(counts) != 0:
        tmp =data[i][0:4] + [alt_genotype[0]] + counts[0:2] + data[i][6:11]
      if (tmp[4] == "<NON_REF>") and tmp[6]!=0: # removes ambiguous cases, likely indels and non-calls
        continue
        #print(tmp)
    # now do the switcheroo
      if (tmp[3] == tmp[10]) and (tmp[4] == tmp[11] or tmp[4] == "<NON_REF>"):
        result += [ [tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]] ] # no need to keep showing NA17878 columns here
      if (tmp[3] == tmp[11]) and (tmp[4] == tmp[10] or tmp[4] == "<NON_REF>"): # debug with line below
        result += [ [tmp[0],tmp[1],tmp[2],tmp[4],tmp[3],tmp[6],tmp[5]] ]#tmp[7],tmp[8],tmp[9],tmp[10],tmp[11]]]

    #print(result)
    with open(arguments.out_directory+os.path.basename(arguments.bed).rstrip(".bed")+".haplotypes.bed","w") as f:
      writer = csv.writer(f,delimiter="\t")
      writer.writerows(result)
      f.close()


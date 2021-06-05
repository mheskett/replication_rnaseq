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

    # data = [[x[0],int(x[1]),int(x[2]),x[3],x[4],int(x[5]),
    #           int(x[6]),x[7],int(x[8]),int(x[9]),x[10],x[11]] for x in data]

    # print(len(data))
    # data = [x if (x[3]==x[10] and x[4]==x[11]) else [ x[0],x[1],x[2],x[4],x[3],x[6],x[5],x[7],x[8],x[9],x[10],x[11] ] for x in data]
    # qc = [x for x in data if (x[3]==x[10] and x[4]==x[11])]
    # qc = [[x[0],x[1],x[2],x[3],x[4],x[5],x[6]] for x in qc]
    # print(len(data))
    # print(len(qc))

    with open(arguments.out_directory+os.path.basename(arguments.bed).rstrip(".bed")+".haplotype.resolved.bed","w") as file:
      writer = csv.writer(file,delimiter="\t")
      writer.writerows(qc)
      file.close()








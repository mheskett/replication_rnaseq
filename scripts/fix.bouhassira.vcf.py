import csv
from sys import argv
"""
The GT (genotype) field encodes allele values separated by either of / or |.
The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. 
For diploid calls examples could be 0/1, 1|0, or 1/2, etc. / indicates an unphased genotype, and | indicates a phased genotype. 
For phased genotypes, the allele to the left of the bar is haplotype 1, and the allele to the right of the bar is haplotype 2.
file looks like this
#CHROM  -1      POS     REF     ALT     X3_2
"""
with open(argv[1]) as infile:
	lines = infile.readlines()
	data = [x.rstrip().split("\t") for x in lines]

result = []
for i in range(len(data)):
	gt=data[i][5]
	chrom=data[i][0]
	pos1=data[i][1]
	pos2=data[i][2]
	ref=data[i][3]
	alt=data[i][4]
	if gt=="0|0":
		result+=[[chrom,pos1,pos2,ref,ref]]
	elif gt=="0|1":
		result+=[[chrom,pos1,pos2,ref,alt]]
	elif gt=="1|1":
		result+=[[chrom,pos1,pos2,alt,alt]]
	elif gt=="1|0":
		result+=[[chrom,pos1,pos2,alt,ref]]

with open(argv[1].rstrip(".vcf")+".aligned.haplotypes.bed","w") as f:
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(result)




### Create a proper bed file from ensemble GTF file
import re
import csv
result = []
with open("/Users/heskett/replication_rnaseq/annotation.files/Homo_sapiens.GRCh37.87.gtf",'r') as f:
		next(f) # skip header
		for line in f:
			tmp = re.split(r'[;,\t\s]\s*',line.rstrip()) ## splits by semi colon, comma, space, tab followed by any whitespace
			for i in range(len(tmp)):
				if tmp[i] == "gene_id":
					gene_id = tmp[i+1]
				if tmp[i] == "gene_name":
					gene_name = tmp[i+1]
			result += [[ tmp[0], tmp[3], tmp[4], gene_id.strip('\"'), gene_name.strip('\"') , tmp[6],tmp[2]]]
with open("/Users/heskett/replication_rnaseq/annotation.files/Homo_sapiens.GRCh37.87.bed","w") as out:
    d = csv.writer(out, delimiter="\t")
    d.writerows(result)


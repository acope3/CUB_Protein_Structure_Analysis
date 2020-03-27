import re

locus = re.compile(r"\[locus_tag=([A-Z0-9\-]+)\]")

with open("mod_scerevisiae.fasta") as fin,open("scer_id.csv",'w') as out:
	out.write("CDS,Locus\n")
	for line in fin:
		if line[0] == ">":
			line_spt = line.split()
			curr = line_spt[0]
			locus_id = locus.search(line)
			print locus_id
			if locus_id != None:
				out.write(curr[1:]+","+locus_id.group(1)+"\n")



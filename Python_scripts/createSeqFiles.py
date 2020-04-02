import sys

seq = {}

with open(sys.argv[1]) as fin:
	current_id = None
	pep = ''
	for line in fin:
		if line[0] == ">":
			if current_id != None:
				seq[current_id] = pep
			current_id = line.split()[0]
			pep = ''
		else:
			pep = pep + line 

	seq[current_id] = pep

for item in seq.keys():
	pep = seq.get(item)
	with open(sys.argv[2]+"/"+item[1:]+".fasta",'w') as out:
		out.write(item+"\n")
		out.write(pep)

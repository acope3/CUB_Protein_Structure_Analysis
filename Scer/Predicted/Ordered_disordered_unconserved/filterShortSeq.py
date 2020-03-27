import glob
import pandas as pd 


phi = pd.read_csv("Disordered/disordered_phi.csv",header=0,index_col=0)
with open("Disordered/disordered.fasta") as fin, open("Disordered_filt/disordered.fasta",'w') as out:
	prot = []
	seq = ''
	for line in fin:
		if line[0] == ">":
			if seq != '' and len(seq) >= 30:
				out.write(gene)
				out.write(seq+"\n")
				prot.append(gene[1:].strip())
			gene = line
			seq = ''
		else:
			seq += line.strip()
phi_filt = phi.loc[prot,:]
phi_filt.to_csv("Disordered_filt/disordered_phi.csv",sep=",",header=True,index=True,quoting=False)

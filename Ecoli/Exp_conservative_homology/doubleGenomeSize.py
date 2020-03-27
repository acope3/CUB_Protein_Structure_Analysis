import glob
import os
import pandas as pd 

old_dir = "Secondary_structures_begin_end/"
new_dir = "Secondary_structures_begin_end_x2/"

def doubleFile(old_dir,new_dir):
	folders = glob.glob(old_dir+"*")
	for f in folders:
		if not os.path.isfile(f):
			to_make = f.split("/")[-1]
			
			try:
				os.mkdir(new_dir+to_make)
			except:
				print "Directory already exists"
			fasta = glob.glob(f+"/*.fasta")
			phi = glob.glob(f+"/*.csv")
			seq_file = open(fasta[0])
			phi_table = pd.read_csv(phi[0],header=0)
			entire_seq = seq_file.read()
			seq_file.close()
			with open(new_dir+to_make+"/"+fasta[0].split("/")[-1],'w') as out:
				out.write(entire_seq)
				out.write(entire_seq)
			df = pd.concat([phi_table,phi_table],axis=0, join='outer', ignore_index=False, keys=None,levels=None, names=None, verify_integrity=False)
			df.to_csv(new_dir+to_make+"/"+phi[0].split("/")[-1],header=True,index=False,quoting=False)
	return None
doubleFile(old_dir,new_dir)
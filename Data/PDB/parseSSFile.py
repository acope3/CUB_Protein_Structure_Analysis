import pandas as pd
import re
import os
import glob
import sys
import shutil


def createPsiPredLikeFiles():
	with open("ecoli_ss_liberal.txt") as fin:
		seq =''
		for line in fin:
			if line[0] == ">":
				if seq != '':
					spt_seq = list(seq)
					df = pd.DataFrame(data = {"Index":xrange(1,len(seq)+1),"AA":["A" for i in xrange(1,len(seq)+1)],"Str":spt_seq},columns=["Index","AA","Str"])
					df.to_csv("Ecoli_psipred_format_liberal/"+prot_id+".ss2",header=False,index=False,sep =" ")
				seq = ''
				prot_id = line.strip()[1:]
			else:
				seq = seq+line.strip()
		spt_seq = list(seq)
		df = pd.DataFrame(data = {"Index":xrange(1,len(seq)+1),"AA":["A" for i in xrange(1,len(seq)+1)],"Str":spt_seq},columns=["Index","AA","Str"])
		df.to_csv("Ecoli_psipred_format_liberal/"+prot_id+".ss2",header=False,index=False,sep =" ")

def createPsiPredLikeDisorderFiles():
	with open("ecoli_ss_disorder_liberal.txt") as fin:
		seq =''
		for line in fin:
			if line[0] == ">":
				if seq != '':
					spt_seq = list(seq)
					df = pd.DataFrame(data = {"Index":xrange(1,len(seq)+1),"AA":["A" for i in xrange(1,len(seq)+1)],"Str":spt_seq},columns=["Index","AA","Str"])
					df.to_csv("Ecoli_psipred_format_liberal/"+prot_id+".disorder",header=False,index=False,sep =" ")
				seq = ''
				prot_id = line.strip()[1:]
			else:
				seq = seq+line.strip()
		spt_seq = list(seq)
		df = pd.DataFrame(data = {"Index":xrange(1,len(seq)+1),"AA":["A" for i in xrange(1,len(seq)+1)],"Str":spt_seq},columns=["Index","AA","Str"])
		df.to_csv("Ecoli_psipred_format_liberal/"+prot_id+".disorder",header=False,index=False,sep =" ")

def mergeSecStrDisorder():
	with open("ecoli_ss_liberal.txt") as fin:
		for line in fin:
			if line[0] == ">":
				prot_id = line.strip()[1:]
				ss2 = pd.read_csv("Ecoli_psipred_format_liberal/"+prot_id+".ss2",header=-1,sep=" ")
				dis = pd.read_csv("Ecoli_psipred_format_liberal/"+prot_id+".disorder",header=-1,sep=" ")
				new_ss2 = pd.concat([ss2,dis],axis=1,join="inner")
				new_ss2.to_csv("Ecoli_psipred_format_liberal/"+prot_id+".ss2",header=False,index=False,sep =" ")



def pullSecStructures():
	df = pd.read_csv("ecoli_liberal_pid_80_60_coverage.blast",sep="\t",header=0)
	#df.drop_duplicates(subset=["Subject"],keep="first",inplace=True)
	proteins = df.Query.values
	pdb_ids = df.Subject.values
	for i in range(len(pdb_ids)):
		pdb_ids[i] = ":".join(pdb_ids[i].split(":")[:2]+["secstr"])
	df.index = pdb_ids
	is_seq = False
	seq = ''
	prot_id = ''
	already = list(df.index)
	total = 0
	check=df["Subject"].value_counts()
	with open("pdb_secstr.txt") as fin,open("ecoli_ss_liberal.txt",'w') as out:
		for line in fin:
			if ">" in line:
				if seq != '' and is_seq:
					seq = seq.replace(" ","C")
					for j in prot_id:
						out.write(">"+j+"\n")
						out.write(seq)
				line_spt = line.split()
				if line_spt[0][1:] in pdb_ids:
					is_seq = True
					get_seq = df.loc[line_spt[0][1:],"Query"]
					if isinstance(get_seq,str):
						prot_id = [get_seq]
					else:
						prot_id = list(get_seq)
					total += len(prot_id)
					for j in prot_id:
						already.remove(line_spt[0][1:])
				else:
					is_seq = False
				seq = ''

			else:
				if is_seq:
					seq += line


def pullDisorder():
	df = pd.read_csv("ecoli_liberal_pid_80_60_coverage.blast",sep="\t",header=0)
	#df.drop_duplicates(subset=["Subject"],keep="first",inplace=True)
	proteins = df.Query.values
	pdb_ids = df.Subject.values
	for i in range(len(pdb_ids)):
		pdb_ids[i] = ":".join(pdb_ids[i].split(":")[:2]+["disorder"])

	df.index = pdb_ids
	is_seq = False
	seq = ''
	prot_id = ''
	already = list(df.index)
	total = 0
	check=df["Subject"].value_counts()
	with open("pdb_disorder.txt") as fin,open("ecoli_ss_disorder_liberal.txt",'w') as out:
		for line in fin:
			if ">" in line:
				if seq != '' and is_seq:
					#seq = seq.replace(" ","C")
					for j in prot_id:
						out.write(">"+j+"\n")
						out.write(seq)
				line_spt = line.split()
				if line_spt[0][1:] in pdb_ids:
					is_seq = True
					get_seq = df.loc[line_spt[0][1:],"Query"]
					if isinstance(get_seq,str):
						prot_id = [get_seq]
					else:
						prot_id = list(get_seq)
					total += len(prot_id)
					for j in prot_id:
						already.remove(line_spt[0][1:])
				else:
					is_seq = False
				seq = ''

			else:
				if is_seq:
					seq += line



def pullSequences():
	df = pd.read_csv("ecoli_liberal_pid_80_60_coverage.blast",sep="\t",header=0)
	df.set_index(keys="Query",inplace=True,drop=False)
	proteins = df.Query.values
	is_seq = False
	seq = ''
	prot_id = ''
	seq_len = {}
	with open("Ecoli_K12_ncbi_protein.fasta") as fin,open("ecoli_exper_evidence_liberal.faa",'w') as out:
		for line in fin:
			if ">" in line:
				if seq != '' and is_seq:
					seq_len[prot_id.split()[0][1:]] = len(re.sub(r"\s+", "", seq))
					out.write(prot_id)
					out.write(seq)
				line_spt = line.split()
				if line_spt[0][1:] in proteins:
					is_seq = True
				else:

					is_seq = False
				prot_id = line
				seq = ''
			else:
				if is_seq:
					seq += line
	# tmp = pd.DataFrame({"Gene_length":pd.Series(seq_len)})
	# print tmp
	# df = df.join(tmp)
	# df.to_csv("test.blast",sep="\t",header=True,index=False,quoting=False)



def createSecStrFile():
	with open("ss_dis.txt") as fin,open("pdb_secstr.txt",'w') as out:
		is_seq = True
		seq = ''
		pdb_id = ''
		for line in fin:
			if ">" in line:
				if "secstr" in line:
					if seq != '':
						out.write(pdb_id)
						out.write(seq)
					pdb_id =line
					is_seq = True
					seq = ''
				else:
					is_seq = False
			else:
				if is_seq:
					seq += line

def createSeqFasta():
	with open("ss_dis.txt") as fin,open("pdb_sequences.fasta",'w') as out:
		is_seq = True
		seq = ''
		pdb_id = ''
		for line in fin:
			if ">" in line:
				if "sequence" in line:
					if seq != '':
						out.write(pdb_id)
						out.write(seq)
					pdb_id =line
					is_seq = True
					seq = ''
				else:
					is_seq = False
			else:
				if is_seq:
					seq += line

def createDisorder():
	with open("ss_dis.txt") as fin,open("pdb_disorder.txt",'w') as out:
		is_seq = True
		seq = ''
		pdb_id = ''
		for line in fin:
			if ">" in line:
				if "disorder" in line:
					if seq != '':
						out.write(pdb_id)
						out.write(seq)
					pdb_id =line
					is_seq = True
					seq = ''
				else:
					is_seq = False
			else:
				if is_seq:
					seq += line



def comparePredictions(test_case):
	turn_as_coil = 0
	turn_as_helix = 0
	turn_as_sheet = 0
	total = 0
	with open("ecoli_ss_liberal.txt") as fin:
		exp_seq = ''
		for line in fin:
			if line[0] == ">":
				if exp_seq != '':
					limit = min([len(pred_ss),len(exp_seq)])
					for i in range(limit):
						if exp_seq[i] in test_case:
							if pred_ss[i] == "H":
								turn_as_helix += 1
							elif pred_ss[i] == "C":
								turn_as_coil +=1
							elif pred_ss[i] == "E":
								turn_as_sheet += 1
				exp_seq = ''
				prot_id = line[1:].strip()
				df = pd.read_csv("../Alpha_psipred/"+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=2)
				pred_ss= "".join(df.iloc[:,2].values)
			else:
				exp_seq += line.strip()
	total = float(turn_as_coil + turn_as_sheet + turn_as_helix)
	print test_case,"identified as alpha helix:",turn_as_helix/total
	print test_case,"identified as beta sheet:",turn_as_sheet/total
	print test_case,"identified as coil:",turn_as_coil/total



def filterPdbSSFile():
	df = pd.read_csv("ecoli_pdb.txt",sep=",",header=0)
	df_filter = df.loc[(df.loc[:, "Percent_positive"] >= 85)  & (df.loc[:,"Gaps"] == 0) & (((df.loc[:,"Subject_end"] - df.loc[:,"Subject_start"])/df.loc[:,"Query_length"] >= 0.6)),:]
	seq_similarity = []
	rows = []
	
	with open("bc-100.out") as fin:
		for line in fin:
			line_spt = line.strip().split()
			seq_similarity.append(line_spt)
		proteins = df_filter["Query"].unique()
		for p in proteins:
			tmp = df_filter.loc[df_filter.Query == p,:]
			if len(tmp.index) == 1:
				best_fitting = tmp.Subject.values
			else:
				structures = tmp.Subject.values
				new_structures = [None for i in structures]
				for j in xrange(len(structures)):
					new_structures[j] = structures[j].split(":sequence")[0]
					new_structures[j] = new_structures[j].replace(":","_")
				best_fitting = None
				for c in seq_similarity:
					if new_structures[0] in c:
						for struct in c:
							if struct in new_structures:
								best_fitting = struct.replace("_",":")+":sequence"
								break
					if best_fitting != None:
						break
			results = tmp.loc[tmp.loc[:,"Subject"] == best_fitting,:]
			if len(results.index) > 1:
				print results
				results = results.iloc[0,:]
		#print results.isnull().values.any(),p
			rows.append(results)
					# assert check == p

	rows = pd.concat(rows,ignore_index=True)
	rows.drop(axis=1,labels=0,inplace=True)
	rows.dropna(axis=0,how="all",inplace=True)

	rows.to_csv("ecoli_liberal_ppos_85_80_coverage.blast",sep="\t",header=True,index=False,quoting=False,columns=["Query","Subject","Query_length","Subject_length","Query_start","Query_end","Subject_start","Subject_end","Evalue","Bitscore","Length","Percent_identity","Percent_positive","Gaps"])

def conservativeFilterPdbSSFile():
	df = pd.read_csv("ecoli_pdb.blast",sep="\t",header=0)
	df_filter = df.loc[(df.loc[:, "Percent_identity"] >= 95) & (df.loc[:,"Gaps"] == 0) & (((df.loc[:,"Subject_end"] - df.loc[:,"Subject_start"])/df.loc[:,"Query_length"] >= 0.8)) ,:]
	seq_similarity = []
	rows = []
	with open("bc-100.out") as fin:
		for line in fin:
			line_spt = line.strip().split()
			seq_similarity.append(line_spt)
		proteins = df_filter["Query"].unique()
		for p in proteins:
			tmp = df_filter.loc[df_filter.Query == p,:]
			if len(tmp.index) == 1:
				best_fitting = tmp.Subject.values
			else:
				structures = tmp.Subject.values
				new_structures = [None for i in structures]
				for j in xrange(len(structures)):
					new_structures[j] = structures[j].split(":sequence")[0]
					new_structures[j] = new_structures[j].replace(":","_")
				best_fitting = None
				for c in seq_similarity:
					if new_structures[0] in c:
						for struct in c:
							if struct in new_structures:
								best_fitting = struct.replace("_",":")+":sequence"
								break
					if best_fitting != None:
						break
			results = tmp.loc[tmp.loc[:,"Subject"] == best_fitting,:]
			if len(results.index) > 1:
				print(results.index)
				results = results.iloc[0,:]

		#print results.isnull().values.any(),p
			rows.append(results)
					# assert check == p
	#print rows

	rows = pd.concat(rows,ignore_index=True)
	#print(rows.columns)
	#rows.drop(axis=1,labels=0,inplace=True)
	rows.dropna(axis=0,how="all",inplace=True)
	rows.to_csv("ecoli_conservative_pid_95_80_coverage.blast",sep="\t",header=True,index=False,quoting=False,columns=["Query","Subject","Query_length","Subject_length","Query_start","Query_end","Subject_start","Subject_end","Evalue","Bitscore","Length","Percent_identity","Percent_positive","Gaps"])


def liberalFilterPdbSSFile():
	df = pd.read_csv("ecoli_pdb.blast",sep="\t",header=0)
	df_filter = df.loc[(df.loc[:, "Percent_identity"] >= 80) & (df.loc[:,"Gaps"] == 0) & (((df.loc[:,"Subject_end"] - df.loc[:,"Subject_start"])/df.loc[:,"Query_length"] >= 0.6)) ,:]
	seq_similarity = []
	rows = []
	with open("bc-100.out") as fin:
		for line in fin:
			line_spt = line.strip().split()
			seq_similarity.append(line_spt)
		proteins = df_filter["Query"].unique()
		for p in proteins:
			tmp = df_filter.loc[df_filter.Query == p,:]
			if len(tmp.index) == 1:
				best_fitting = tmp.Subject.values
			else:
				structures = tmp.Subject.values
				new_structures = [None for i in structures]
				for j in xrange(len(structures)):
					new_structures[j] = structures[j].split(":sequence")[0]
					new_structures[j] = new_structures[j].replace(":","_")
				best_fitting = None
				for c in seq_similarity:
					if new_structures[0] in c:
						for struct in c:
							if struct in new_structures:
								best_fitting = struct.replace("_",":")+":sequence"
								break
					if best_fitting != None:
						break
			results = tmp.loc[tmp.loc[:,"Subject"] == best_fitting,:]
			if len(results.index) > 1:
				print(results.index)
				results = results.iloc[0,:]

		#print results.isnull().values.any(),p
			rows.append(results)
					# assert check == p
	#print rows

	rows = pd.concat(rows,ignore_index=True)
	#print(rows.columns)
	#rows.drop(axis=1,labels=0,inplace=True)
	rows.dropna(axis=0,how="all",inplace=True)
	rows.to_csv("ecoli_liberal_pid_80_60_coverage.blast",sep="\t",header=True,index=False,quoting=False,columns=["Query","Subject","Query_length","Subject_length","Query_start","Query_end","Subject_start","Subject_end","Evalue","Bitscore","Length","Percent_identity","Percent_positive","Gaps"])

def readPechmann():
	files = glob.glob("Pechmann_data/profiles_PDB/*.dat")
	id_map = pd.read_csv("ecolievisiae_id_mapping.tsv",sep="\t",header=0,index_col=2)
	data_available = []
	for f in files:
		locus = f.split("/")[2].split(".")[0]
		prot = id_map.loc[locus,"Protein"]
		if prot != None:
			if os.path.exists("Psipred_format_liberal/"+prot+".ss2"):
				shutil.copy("Psipred_format_liberal/"+prot+".ss2","Pechmann_proteins")
		df = pd.read_csv(f,delim_whitespace=True,header=None,names=["Index","Codon","AA","Structure","Score","Conserved"])
		df.replace("-","C",inplace=True)
		df = df.iloc[:,[0,2,3,4,5]]
		df.to_csv("Psipred_format_pechmann/"+prot+".ss2",header=False,index=False,quoting=False,sep=" ")




# createSecStrFile()
# createSeqFasta()
# createDisorder()
liberalFilterPdbSSFile()
pullSecStructures()
pullSequences()
pullDisorder()
createPsiPredLikeFiles()
createPsiPredLikeDisorderFiles()
mergeSecStrDisorder()
# comparePredictions(["H"])
# comparePredictions(["E"])
# comparePredictions(["C","G","B","I"])
# comparePredictions(["T","S"])
#readPechmann()
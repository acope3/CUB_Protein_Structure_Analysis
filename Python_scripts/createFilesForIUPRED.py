import numpy as np
import re
import argparse
import pandas as pd
import glob


class Genome:
	def __init__(self):
		self.fasta_file = ""
		self.genes = {}
		self.id_mapper = {}
		self.phi = None

	## Currently assumes RefSeq cds file format
	def readFasta(self,fasta_file):
		self.fasta_file = fasta_file
		gene_id_pat = re.compile(r"\[gene=[A-Z0-9]+\]")
		locus_id_pat = re.compile(r"\[locus_tag=[A-Z0-9\-]+\]")
		prot_id_pat = re.compile(r"\[protein_id=[A-Z0-9\.]+\]")
		with open(fasta_file) as fin:
			cds_id = ''
			seq = ''
			for line in fin:
				if line[0] == ">":
					if cds_id != '':
						self.genes[cds_id] = seq	
					cds_id = line.split()[0][1:]
					gene_id = gene_id_pat.search(line)
					locus_id = locus_id_pat.search(line)
					prot_id = prot_id_pat.search(line)
					self.id_mapper[cds_id] = [gene_id,locus_id,prot_id]
					seq =''
				else:
					seq += line.strip()
			self.genes[cds_id] = seq

	def getGenes(self):
		return self.genes

	def setGenes(self,new_genes):
		self.genes = new_genes

	def setGenesById(self,new_gene_id,new_gene_seq):
		self.genes[new_gene_id] = new_gene_seq

	def getGeneById(self,gene_id):
		return self.genes.get(gene_id)

	def getCodonsForSeq(self,gene_id):
		seq = genes[gene_id]
		codons = re.findall(r"[ATCGYNW]{3}",gene.strip())
		return codons


	def writeFasta(self,filename):
		with open(filename,'w') as output:
			for i in self.genes:
				seq = self.genes[i]
				output.write(">"+i+"\n")
				output.write(seq+"\n")

	def createPhiTable(self,df_phi):
		gene_ids = self.genes.keys()
		self.phi = df_phi.loc[gene_ids,:]

	def getPhiTable(self):
		return self.phi

	def setPhiTable(self,phi_table):
		self.phi = phi_table

	def writePhiTable(self,filename):
		self.phi.to_csv(filename,header=True,index=True,sep=",",quoting=False)

	def mergeWithOtherGenome(self,genome):
		new_genome = Genome()
		all_ids_self = self.genes.keys()
		all_ids_other = genome.getGenes().keys()
		union_ids = set(all_ids_self).union(set(all_ids_other))
		for i in union_ids:
			seq = self.genes.get(i)
			other_seq = genome.getGeneById(i)
			if seq != None:
				seq = seq.strip()[3:-3] # remove start and stop codon which were added, and newline character
			else:
				seq = ''
			if other_seq != None:
				other_seq = other_seq.strip()[3:-3]
			else:
				other_seq = ''
			merged_seq = "".join([seq,other_seq])
			new_genome.setGenesById(i,merged_seq)
		all_phi = self.phi.append(genome.getPhiTable())
		all_phi = all_phi.loc[~all_phi.index.duplicated(keep='first')]
		all_phi = all_phi.loc[new_genome.getGenes().keys(),:]
		new_genome.setPhiTable(all_phi)
		return new_genome


def parseIUPREDFiles(iupred):
	with open(iupred) as fin:
		f=None
		for line in fin:
			if line[0] == ">":
				if f != None:
					f.close()
				f = open("../Data/IUPRED2_ecoli/"+line.split()[0][1:].strip(),'w')
			elif line[0]!="#":
				f.write(line)
		f.close()
	return None

# genome = Genome()
# genome.readFasta("../mod_scerevisiae.fasta")
# genes = genome.getGenes()

# tmp_genes = {}
# current=1
# for i in genes:
# 	if len(tmp_genes) == 650:
# 		tmp = Genome()
# 		tmp.setGenes(tmp_genes)
# 		tmp.writeFasta("scer_"+str(current)+".fasta")
# 		tmp_genes={}
# 		current+=1
# 	seq = genes.get(i)
# 	tmp_genes[i] = seq

# tmp = Genome()
# tmp.setGenes(tmp_genes)
# tmp.writeFasta("scer_"+str(current)+".fasta")

parseIUPREDFiles("../Data/Ecoli_K12_ncbi_protein.result")


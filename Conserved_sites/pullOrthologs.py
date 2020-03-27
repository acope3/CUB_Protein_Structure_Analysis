import numpy as np
import re
import argparse
import pandas as pd
import glob

class Genome:
	def __init__(self):
		self.fasta_file = ""
		self.faa_file = ""
		self.genes = {}
		self.proteins = {}
		self.id_mapper = {}
		self.phi = None

	## Currently assumes RefSeq cds file format
	def readFNA(self,fasta_file):
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

	def readFAA(self,fasta_file):
		self.faa_file = fasta_file
		with open(fasta_file) as fin:
			prot_id = ''
			seq = ''
			for line in fin:
				if line[0] == ">":
					if prot_id != '':
						if seq[-1]=="*":
							seq=seq[0:-1]
						self.proteins[prot_id] = seq	
					prot_id = line.split()[0][1:]
					seq =''
				else:
					seq += line.strip()
			self.proteins[prot_id] = seq

	def getGenes(self):
		return self.genes

	def setGenes(self,new_genes):
		self.genes = new_genes

	def setGenesById(self,new_gene_id,new_gene_seq):
		self.genes[new_gene_id] = new_gene_seq

	def getGeneById(self,gene_id):
		return self.genes.get(gene_id)

	def getProteins(self):
		return self.proteins

	def getProteinById(self,prot_id):
		return self.proteins.get(prot_id)


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


def pullOrthologs():

	target_directory = "Orthologs_sequences_proteins_sacch/"

	scer = Genome()
	sbay = Genome()
	spar = Genome()
	smik = Genome()
	skud = Genome()
	ncast = Genome()
	scer.readFAA("scerevisiae_R64_protein.faa")
	sbay.readFAA("sbayanus.faa")
	spar.readFAA("sparadoxus.faa")
	smik.readFAA("smikatae.faa")
	skud.readFAA("skudriavzevii.faa")
	ncast.readFAA("ncastellii.faa")
	with open("sacch_ortho_table.tsv") as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split()
			scer_prot_id = line_spt[0]
			ncast_prot_id = line_spt[1]
			sbay_prot_id = line_spt[2]
			skud_prot_id = line_spt[3]
			smik_prot_id = line_spt[4]
			spar_prot_id = line_spt[5]
			scer_prot = scer.getProteinById(scer_prot_id)
			ncast_prot = ncast.getProteinById(ncast_prot_id)
			sbay_prot = sbay.getProteinById(sbay_prot_id)
			skud_prot = skud.getProteinById(skud_prot_id)
			smik_prot = smik.getProteinById(smik_prot_id)
			spar_prot = spar.getProteinById(spar_prot_id)
			with open(target_directory+scer_prot_id+".faa",'w') as out:
				out.write(">"+scer_prot_id+"\n")
				out.write(scer_prot)
				out.write("\n>"+ncast_prot_id+"\n")
				out.write(ncast_prot)
				out.write("\n>"+sbay_prot_id+"\n")
				out.write(sbay_prot)
				out.write("\n>"+skud_prot_id+"\n")
				out.write(skud_prot)
				out.write("\n>"+smik_prot_id+"\n")
				out.write(smik_prot)
				out.write("\n>"+spar_prot_id+"\n")
				out.write(spar_prot)

def pullConservedSites():
	alignments = glob.glob("Protein_alignments_sacch/*")
	for orthologs in alignments:
		alignment = Genome()
		alignment.readFAA(orthologs)
		proteins = orthologs.split("/")[1].split(".faa")
		prot_1 = proteins[0]
		genes = alignment.getProteins()
		scer_seq = genes.get(prot_1)
		len_alignment = float(len(scer_seq))
		num_alignment = len(genes)
		gene_ids = genes.keys()
		seq = []
		for i in gene_ids:
			if i != prot_1:
				seq.append(list(genes.get(i)))
		# missing_1 = scer_seq.count("-")
		# missing_2 = sbay_seq.count("-")
		# gap_1 = missing_1/len_alignment
		# gap_2 = missing_2/len_alignment
		# if ((gap_1+gap_2) <= 0.3):
		with open("Conserved_positions_sacch/"+prot_1,'w') as out,open("Unconserved_positions_sacch/"+prot_1,'w') as out_2:
			position = 0
			for i,j,k,l,m,n in zip(list(scer_seq),*seq):
				if i != "-":
					if i == j and i == k and i == l and i == m and i == n:
						out.write(str(position)+"\n")
					else:
						out_2.write(str(position)+"\n")
					position+=1

#pullOrthologs()
pullConservedSites()


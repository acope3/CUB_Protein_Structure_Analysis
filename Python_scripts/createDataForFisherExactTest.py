## Author: Alexander Cope
## Create files from simulated data (see R/simulateRegions.R) to repeat Pechmann and Frydman 2013 analysis of position-specific codon biases in helices.



import re
import pandas as pd
import numpy as np
import os
import sys
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

def pullStartHelixForFisherExact(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="../Data/Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",experimental=True):
    blast_results=None
    if experimental:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin,open(output_folder+"helix.csv",'w') as out:
        out.write("Gene,Pos,Codon\n")
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                value = re.findall(pat,line.strip())   
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",delim_whitespace=True,skiprows=num_rows,header=None)
                    if experimental:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except (IOError,pd.errors.EmptyDataError):
                    continue
                starth_pat = re.compile(r"[HGI]{5}[HGI]+")
                ss = "".join(df.iloc[:,2].values)
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                starth_transitions = starth_pat.finditer(ss)
                offset = qstart - sstart
                starth_locations = []
                for i in starth_transitions:
                    begin = i.start()+offset
                    last = i.end()+offset
                    pos = 1
                    for position,codon in zip(xrange(begin,last),value[begin:last]):
                        out.write(",".join([prot_id,str(pos),codon]))
                        out.write("\n")
                        pos+=1        

for i in xrange(1,101):
    genome = Genome()
    genome.readFasta("../Scer/Test_Fisher_Exact/"+str(i)+"/sim.fasta")
    genome.writeFasta("../Scer/Test_Fisher_Exact/"+str(i)+"/mod_sim.fasta")
    pullStartHelixForFisherExact("../Scer/Test_Fisher_Exact/"+str(i)+"/",genome_file="../Scer/Test_Fisher_Exact/"+str(i)+"/mod_sim.fasta",ss2_loc = "../Data/Exp_structure_data/Scer_psipred_format_conservative/")

pullPredScer("../Scer/Predicted/",ss2_loc="../Data/Alpha_psipred/")

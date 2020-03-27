import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
import os
import sys
import glob
import seaborn as sns


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



def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))



def testMovingAverage():
    pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    dEta = pd.read_csv("selection_mod_scerevisiae.csv",sep=",",header=0)
    with open("mod_scerevisiae.fasta") as fin:
        for line in fin:
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
            if line[0] != ">":
                df = pd.read_csv("Orig_psipred/"+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=2)
                df["Structure"] = np.where(df.iloc[:,2]!='C', 0, 1)
                movingAverage(line.strip(),dEta,df,prot_id,window=8,output="CUB_variation_graphs/"+prot_id+".pdf")

def movingAverage(gene_seq,dEta,structure,prot_id,window=8,output="CodonVariation.pdf"):
    pat = re.compile(r"[ATCGYNW]{3}")
    value = re.findall(pat,gene_seq.strip())
    df_seq = pd.DataFrame(data={"dEta":[0.0 for i in value[:-1]],"Codon":value[:-1]},)
    df_seq.index = df_seq.index
    df_seq = pd.merge(df_seq,dEta,how="left",left_on="Codon",right_on="Codon")
    df_seq = df_seq.fillna(0)
    df_seq["MovingAvg"] = df_seq.loc[:,"Posterior"].rolling(window=window).mean()
    df_seq["Structure"] = structure.loc[:,"Structure"]
    df_seq["MovingAvgStructure"] = df_seq.loc[:,"Structure"].rolling(window=window).mean()
    fig,ax1 = plt.subplots(1,1)
    ax1.plot(range(len(df_seq.index)),df_seq.loc[:,"MovingAvg"])
    ax2 = ax1.twinx()
    ax2.plot(range(len(df_seq.index)),df_seq.loc[:,"MovingAvgStructure"],"red")
    print "Made fig"
    ax1.set_title(prot_id)
    ax1.set_xlabel("Codon Position")
    ax1.set_ylabel("Moving Average (DeltaEta)")
    ax2.set_ylabel("Moving Average (Structure)")
    ax2.set_ylim(bottom=-0.05,top=1.05)
    plt.savefig(output)
    #fig.close()
    plt.close()
    #print df_seq
    #df_seq["MovingAvg"] = df_seq.rolling(window=window).mean()
    #print df_seq


def getCutoff(cutoff):
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    helix = []
    sheet = []
    coil = []
    with open("mod_scerevisiae.fasta") as fin:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
            if line[0] != ">":
                value = re.findall(pat,line.strip())
                try:
                    df = pd.read_csv("Orig_psipred/"+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=2)
                    h = df.loc[(df.iloc[:,2] == "H"),4].values
                    s = df.loc[(df.iloc[:,2] == "E"),5].values
                    c = df.loc[(df.iloc[:,2] == "C"),3].values
                    helix += list(h)
                    sheet += list(s)
                    coil += list(c)
                except IOError:
                    "File does not exist"
                    continue
    limit = {}
    limit["helix"] = np.percentile(helix,q=cutoff)
    limit["sheet"] = np.percentile(sheet,q=cutoff)
    limit["coil"] = np.percentile(coil,q=cutoff)
    return limit

def getSecondaryStructureConf():
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    helix = []
    sheet = []
    coil = []
    with open("mod_scerevisiae.fasta") as fin:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
            if line[0] != ">":
                value = re.findall(pat,line.strip())
                try:
                    df = pd.read_csv("Orig_psipred/"+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=2)
                    h = df.loc[(df.iloc[:,2] == "H"),4].values
                    s = df.loc[(df.iloc[:,2] == "E"),5].values
                    c = df.loc[(df.iloc[:,2] == "C"),3].values
                    helix += list(h)
                    sheet += list(s)
                    coil += list(c)
                except IOError:
                    "File does not exist"
                    continue

    fig, ax = plt.subplots(2,2,figsize=(12,12),sharex=False)
    #max_helix = max(helix)
    #helix.remove(max_helix)
    #max_coil = max(coil)
    #coil.remove(max_coil)

    dist = sns.distplot(helix,color="blue",ax=ax[0,0],axlabel="Confidence",kde=False)
    dist = sns.distplot(sheet,color="red",ax=ax[0,1],axlabel="Confidence",kde=False)
    dist = sns.distplot(coil,color="grey",ax=ax[1,0],axlabel="Confidence",kde=False)
    ax[0,0].set_title("Helix")
    ax[0,1].set_title("Sheet")
    ax[1,0].set_title("Coil")
    fig.suptitle("Distribution of Confidence for Secondary structure Predictions")

    dist.get_figure().savefig("confidence_score_distributions.pdf")
    return None


def getSecondaryStructureLengthDist(experimental=False,remove_disorder=False):
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    helix = []
    sheet = []
    coil = []
    turn = []
    qstart = 1
    qend = -1
    if experimental:
        blast_results = pd.read_csv("Exp_structure_data/ecoli_conservative_pid_95_80_coverage.blast",sep="\t",header=0,index_col=0)
        #blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    with open("mod_scerevisiae.fasta") as fin:
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
            if line[0] != ">" and qstart!=-1:
                value = re.findall(pat,line.strip())
                if experimental:
                    value = value[qstart:qend]
                else:
                #cds file goes to stop codon, but structural data does not, so exclude stop codon
                    value = value[0:-1]
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder:
                        df.loc[df.iloc[:,5] == "X",2] = "X"
                except IOError:
                    print "Could not find ",prot_id+".ss2"
                    continue
                helix_pat = re.compile(r"H+")
                sheet_pat = re.compile(r"E+")
                coil_pat = re.compile(r"[CBIG]+")
                turn_pat = re.compile(r"[TS]+")
                ss = "".join(df.iloc[:,2].values)
                length = len(ss)
                helix_results = helix_pat.findall(ss)
                sheet_results = sheet_pat.findall(ss)
                coil_results = coil_pat.findall(ss)
                turn_results = turn_pat.findall(ss)
                total_length = 0
                for i in helix_results:
                    if len(i) >= 100:
                        helix.append(100)
                    else:
                        helix.append(len(i))
                    total_length += len(i)
                for i in sheet_results:
                    if len(i) >= 100:
                        sheet.append(100)
                    else:
                        sheet.append(len(i))
                    total_length += len(i)
                for i in coil_results:
                    if len(i) >= 100:
                        coil.append(100)
                    else:
                        coil.append(len(i))
                    total_length += len(i)
                for i in turn_results:
                    if len(i) >= 100:
                        turn.append(100)
                    else:
                        turn.append(len(i))
                    total_length += len(i)

                assert total_length == len(ss)
    # fig, ax = plt.subplots(2,2,figsize=(12,12),sharex=False)
    print np.median(helix)
    print np.var(helix)

    print np.median(sheet)
    print np.var(sheet)

    print np.median(coil)
    print np.var(coil)

    print np.median(turn)
    print np.var(turn)
    #max_helix = max(helix)
    #helix.remove(max_helix)
    #max_coil = max(coil)
    #coil.remove(max_coil)
    # dist = sns.distplot(helix,bins=[0,10,20,30,40,50,60,70,80,90,100],color="blue",ax=ax[0,0],axlabel="Length",kde=False)
    # dist = sns.distplot(sheet,bins=[0,10,20,30,40,50,60,70,80,90,100],color="red",ax=ax[0,1],axlabel="Length",kde=False)
    # dist = sns.distplot(coil,bins=[0,10,20,30,40,50,60,70,80,90,100],color="grey",ax=ax[1,0],axlabel="Length",kde=False)
    # ax[0,0].set_title("Helix")
    # ax[0,0].set_xticklabels(["0","20","40","60","80",">=100"])
    # ax[0,1].set_title("Sheet")
    # ax[0,1].set_xticklabels(["0","20","40","60","80",">=100"])
    # ax[1,0].set_title("Coil")
    # ax[1,0].set_xticklabels(["0","20","40","60","80",">=100"])
    # fig.suptitle("Distribution of lengths of Secondary_structures")

    # dist.get_figure().savefig("test.pdf")
    return None


def structuresRepresented(output,genome_file = "mod_scerevisiae.fasta",ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,exclude_mito=False,only_chap=False,exclude_tm_sp=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
        blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] == 100.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    qstart = 0
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    with open(genome_file) as fin, \
    open(output,'w') as out:
        out.write("Protein\tCoil\tHelix\tSheet\tLength\n")
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                value = re.findall(pat,line.strip())
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    post_dict[position] = codon
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except IOError:
                    "File does not exist"
                    continue
                ss = "".join(df.iloc[:,2].values)
                ss=ss.replace("T","C")
                ss=ss.replace("S","C")
                ss=ss.replace("I","C")
                ss=ss.replace("G","C")
                ss=ss.replace("B","C")
                length = len(ss[sstart:send])
                ss = ss[sstart:send]
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                percent_coil = ss.count("C")/float(length)
                percent_helix = ss.count("H")/float(length)
                percent_sheet = ss.count("E")/float(length)
                out.write("\t".join([cds_id[1:],str(percent_coil),str(percent_helix),str(percent_sheet),str(length)]))
                out.write("\n")


def classifyCodon(position,ss):
    if (position == 0) or (position == len(ss)-1):
        return None
    else:
        pattern = ss[position-1:position+2]
        if pattern == "CHH" or pattern == "EHH":
            return "H1"
        elif pattern == "HHH":
            return "H2"
        elif pattern == "HHC" or pattern == "HHE":
            return "H3"
        elif pattern == "CEE" or pattern == "HEE":
            return "E1"
        elif pattern == "EEE":
            return "E2"
        elif pattern == "EEC" or pattern == "EEH":
            return "E3"
        elif pattern == "HCC" or pattern == "ECC":
            return "C1"
        elif pattern == "CCC":
            return "C2"
        elif pattern == "CCH" or pattern == "CCE":
            return "C3"
        else:
            return None

def createContingencyTable(output,genome_file = "mod_scerevisiae.fasta",ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,exclude_mito=False,only_chap=False,exclude_tm_sp=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
        #blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    qstart = 0
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    with open(genome_file) as fin, \
    open(output,'w') as out:
        out.write("Protein\tPosition\tCodon\tClass\n")
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                value = re.findall(pat,line.strip())
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    post_dict[position] = codon
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except IOError:
                    "File does not exist"
                    continue
                ss = "".join(df.iloc[:,2].values)
                ss=ss.replace("T","C")
                ss=ss.replace("S","C")
                ss=ss.replace("I","C")
                ss=ss.replace("G","C")
                ss=ss.replace("B","C")
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                seq = []
                for i,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if i < 100 or qend - i < 100:
                        continue
                    ##print ss[i-offset-1:i-offset+2]
                    classified = classifyCodon(i-offset,ss)
                    if classified != None:
                        out.write("\t".join([prot_id,str(i),codon,classified]))
                        out.write("\n")
    df = pd.read_csv(output,header=0,sep="\t")
    deta = pd.read_csv("scerevisiae_Selection.csv",header=0)
    df = df.merge(deta,on="Codon")
    df = df.drop(columns=["Position","Posterior","0.025%","0.975%"])
    df.to_csv(output,sep="\t",header=True,index=False,quoting=False)

                    


def structureDownstreamCat_simplified(output_folder,struct = "C",genome_file = "mod_scerevisiae.fasta",ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,genes_used=None,length_cutoff=1,exclude_mito=False,only_chap=False,exclude_tm_sp=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    type_dict = {"C":"coil","H":"helix","E":"sheet","T":"turn","X":"disordered"}
    equiv_types = {"C":["C","G","I","B"],"H":["H"],"E":["E"],"T":["S","T"],"X":["X"],"*":["C","G","B","I","T","S","H","E"]}
    if len(struct) == 1:
        if struct == "N":
            type = "nterminus"
        if struct == "X":
            type = "disordered"
    else:
        if struct == "X*":
            type = "downstream_disordered"
        else:
            type = type_dict.get(struct[0])+"_"+type_dict.get(struct[1])
    subfolder = type.capitalize()
    fasta_file = type+".fasta"
    if not os.path.exists(output_folder+"/"+subfolder):
        os.makedirs(output_folder+"/"+subfolder)
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
        #blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    qstart = 0
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    output = "/".join([output_folder,subfolder,fasta_file])
    with open(genome_file) as fin, \
    open(output,'w') as out:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                value = re.findall(pat,line.strip())
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    post_dict[position] = codon
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except IOError:
                    "File does not exist"
                    continue
                ss = "".join(df.iloc[:,2].values)
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                seq = []
                for i,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    sec_str = ss[i-offset]
                    if sec_str == "X":
                        continue
                    elif type == "nterminus":
                        if (i - 35) <  0:
                            seq.append(codon)
                        else:
                            continue
                    else:
                        if (i - 35) >= 0:
                            upstream_type = equiv_types.get(struct[0])
                            sec_type = equiv_types.get(struct[1])
                            upstream_str = ss[i-35-offset]
                            if (upstream_str in upstream_type) and (sec_str in sec_type):
                                seq.append(codon)
                seq = "".join(seq)
                if len(seq) >= length_cutoff:
                    writeToFile(seq,cds_id,out)
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[type]},columns=["Gene","Class"]))
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return(genes_used)



def writeToFile(seq,prot_id,file_var):
    if seq:
        file_var.write(prot_id+"\n")
        file_var.write(seq+"\n")



def createExpressionFiles(output_folder,genes_used,gene_expression_file="mod_scerevisiae_expression_wo_phi_allunique.csv"):
    files = glob.glob(output_folder+"/*/*.fasta")
    phi = pd.read_csv(gene_expression_file,header=0,index_col=0)
    for f in files:
        split_dir = f.split("/")
        subfolder = split_dir[-2]
        #print subfolder
        type_struct = split_dir[-1].split(".fasta")[0]
        struct_genes = genes_used.loc[genes_used.iloc[:,1] == type_struct,"Gene"]
        phi_struct = phi.loc[struct_genes,:]
        phi_struct.to_csv(output_folder+"/"+subfolder+"/"+type_struct+"_phi.csv",sep=",",header=True,index=True,quoting=False)
    return None

def filterShortSeq():
    with open("mod_scerevisiae.fasta") as fin,open("mod_scerevisiae_200_codons_minimum.fasta",'w') as out:
        seq_id = ""
        for line in fin:
            if line[0] == ">":
                seq_id = line
            else:
                if len(line.strip())/3 >= 200:
                    out.write(seq_id)
                    out.write(line)

def findTransitions(position,locations):
    for transition in locations:
        if ((transition[0]) <= position) and ((transition[1]) > position):
            return True
    return False

def findDownstream(position,locations,length,downstream_pos,transition_size):
    for transition in locations:
        ts = transition_size
        if transition[1] - transition[0] < transition_size:
            ts = transition[1] - transition[0]
        if ((transition[1] - 1 + downstream_pos) <= position) and ((transition[1] -1 + downstream_pos + ts) > position) and ((transition[1]-1+downstream_pos+ts) < length):
            return True
    return False


def getProtsForAnalysis(exclude_chap=False,only_chap=False,exclude_mito=True,exclude_tm_sp=False):
    prot_to_exclude = set([])
    prot_to_include = set([])
    if exclude_chap or only_chap:
        chap_df = pd.read_csv("chaperone_associated.tsv",sep="\t",header=-1)
        if exclude_chap:
            prot_to_exclude= prot_to_exclude.union(list(chap_df.iloc[:,0].values))
        else:
            prot_to_include= prot_to_include.union(list(chap_df.iloc[:,0].values))
    if exclude_tm_sp:
        tm_sp = pd.read_csv("yeast_tm_sp_prot.txt",sep="\n",header=-1)
        prot_to_exclude = prot_to_exclude.union(list(tm_sp.iloc[:,0].values))
    if exclude_mito:
        mito = pd.read_csv("Simulations_BF/mitochodrial_protein_exclude.txt",sep="\n",header=-1)
        prot_to_exclude = prot_to_exclude.union(list(mito.iloc[:,0].values))
    return prot_to_exclude,prot_to_include


def pullSecondaryOrder(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False,residues_to_keep=None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil_ord"):
        os.makedirs(output_folder+"/Coil_ord")
    if not os.path.exists(output_folder+"/Helix_ord"):
        os.makedirs(output_folder+"/Helix_ord")
    if not os.path.exists(output_folder+"/Sheet_ord"):
        os.makedirs(output_folder+"/Sheet_ord")
    if not os.path.exists(output_folder+"/Coil_dis"):
        os.makedirs(output_folder+"/Coil_dis")
    if not os.path.exists(output_folder+"/Helix_dis"):
        os.makedirs(output_folder+"/Helix_dis")
    if not os.path.exists(output_folder+"/Sheet_dis"):
        os.makedirs(output_folder+"/Sheet_dis")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
      #  blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Coil_ord/coil_ord.fasta",'w') as coil_ord,\
    open(output_folder+"/Sheet_ord/sheet_ord.fasta",'w') as sheet_ord,\
    open(output_folder+"/Helix_ord/helix_ord.fasta",'w') as helix_ord,\
    open(output_folder+"/Coil_dis/coil_dis.fasta",'w') as coil_dis,\
    open(output_folder+"/Sheet_dis/sheet_dis.fasta",'w') as sheet_dis,\
    open(output_folder+"/Helix_dis/helix_dis.fasta",'w') as helix_dis,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        if use_emp_seq:
                            sstart = qstart
                            send = qend
                        else:
                            sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                            send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=None,delim_whitespace=True)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                coilo_seq = ''
                helixo_seq = ''
                sheeto_seq = ''
                coild_seq = ''
                helixd_seq = ''
                sheetd_seq = ''
                nterm_seq = ''
                comp_seq = ''
                ss = "".join(df.iloc[:,2].values)
                order = "".join(df.iloc[:,6].values)
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                #offset = 0
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    nterm_class = False
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    if not nterm_class:
                        if ss[position-offset] == "H" and order[position-offset] == "O":
                            helixo_seq += codon
                            if len(helixo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_ord"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "H" and order[position-offset] == "D":
                            helixd_seq += codon
                            if len(helixd_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_dis"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "C" and order[position - offset]=="O":
                            coilo_seq += codon
                            if len(coilo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_ord"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "C" and order[position-offset] == "D":
                            coild_seq += codon
                            if len(coild_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_dis"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "E" and order[position-offset] == "O":
                            sheeto_seq += codon
                            if len(sheeto_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet_ord"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "E" and order[position-offset] == "D":
                            sheetd_seq += codon
                            if len(sheetd_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet_dis"]},columns=["Gene","Class"]))
                writeToFile(helixo_seq,cds_id,helix_ord)
                writeToFile(coilo_seq,cds_id,coil_ord)
                writeToFile(sheeto_seq,cds_id,sheet_ord)
                writeToFile(helixd_seq,cds_id,helix_dis)
                writeToFile(coild_seq,cds_id,coil_dis)
                writeToFile(sheetd_seq,cds_id,sheet_dis)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                assert len(comp_seq) == (len(nterm_seq) +len(helixo_seq) + len(sheeto_seq) + len(coilo_seq) +len(helixd_seq) + len(sheetd_seq) + len(coild_seq))
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used

def pullSecondaryBeginEnd(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,separate_cterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False,residues_to_keep=None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil"):
        os.makedirs(output_folder+"/Coil")
    if not os.path.exists(output_folder+"/Helix"):
        os.makedirs(output_folder+"/Helix")
    if not os.path.exists(output_folder+"/Sheet"):
        os.makedirs(output_folder+"/Sheet")
    if not os.path.exists(output_folder+"/Turn"):
        os.makedirs(output_folder+"/Turn")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    if not os.path.exists(output_folder+"/Cterminus"):
        os.makedirs(output_folder+"/Cterminus")
    if not os.path.exists(output_folder+"/End_coil"):
        os.makedirs(output_folder+"/End_coil")
    if not os.path.exists(output_folder+"/Start_coil"):
        os.makedirs(output_folder+"/Start_coil")
    if not os.path.exists(output_folder+"/Start_helix"):
        os.makedirs(output_folder+"/Start_helix")
    if not os.path.exists(output_folder+"/End_helix"):
        os.makedirs(output_folder+"/End_helix")
    if not os.path.exists(output_folder+"/Start_sheet"):
        os.makedirs(output_folder+"/Start_sheet")
    if not os.path.exists(output_folder+"/End_sheet"):
        os.makedirs(output_folder+"/End_sheet")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
      #  blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/End_coil/end_coil.fasta",'w') as endc,\
    open(output_folder+"/Start_coil/start_coil.fasta",'w') as startc,\
    open(output_folder+"/Start_helix/start_helix.fasta",'w') as starth,\
    open(output_folder+"/End_helix/end_helix.fasta",'w') as endh,\
    open(output_folder+"/Start_sheet/start_sheet.fasta",'w') as starte,\
    open(output_folder+"/End_sheet/end_sheet.fasta",'w') as ende,\
    open(output_folder+"/Coil/coil.fasta",'w') as coil,\
    open(output_folder+"/Sheet/sheet.fasta",'w') as sheet,\
    open(output_folder+"/Helix/helix.fasta",'w') as helix,\
    open(output_folder+"/Turn/turn.fasta",'w') as turn,\
    open(output_folder+"/Cterminus/cterminus.fasta",'w') as cterm,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        if use_emp_seq:
                            sstart = qstart
                            send = qend
                        else:
                            sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                            send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                startc_seq=''
                endc_seq=''
                starth_seq=''
                endh_seq=''
                starte_seq=''
                ende_seq=''
                coil_seq = ''
                helix_seq = ''
                sheet_seq = ''
                turn_seq = ''
                comp_seq = ''
                nterm_seq = ''
                cterm_seq = ''
                dis_seq = ''
                h_pat = re.compile(r"HHHHHHHH+")
                e_pat = re.compile(r"EEEEEEEE+")
                c_pat = re.compile(r"[CBGITS]{7}[CBGITS]+")
                ss = "".join(df.iloc[:,2].values)
                
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                #offset = 0
                h_transitions = h_pat.finditer(ss)
                e_transitions = e_pat.finditer(ss)
                c_transitions = c_pat.finditer(ss)

                startc_locations = []
                endc_locations = []
                starte_locations = []
                ende_locations = []
                starth_locations = []
                endh_locations = []
                for i in h_transitions:
                    endh_locations.append((i.end()-3+offset,i.end()+offset))
                    starth_locations.append((i.start()+offset,i.start()+3+offset))

                for i in e_transitions:
                    ende_locations.append((i.end()-3+offset,i.end()+offset))
                    starte_locations.append((i.start()+offset,i.start()+3+offset))
                for i in c_transitions:
                    endc_locations.append((i.end()-3+offset,i.end()+offset))
                    startc_locations.append((i.start()+offset,i.start()+3+offset))
                all_loc = [startc_locations,endc_locations,starth_locations,endh_locations,starte_locations,ende_locations]
                
                prefixes = ["start_coil","end_coil","start_helix","end_helix","start_sheet","end_sheet"]
        
                seq = [startc_seq,endc_seq,starth_seq,endh_seq,starte_seq,ende_seq]
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if ss[position - offset] == "X":
                        continue
                    if residues_to_keep != None:
                        if position not in keep:
                            continue
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    transit_class = False
                    nterm_class = False
                    cterm_class = False
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    elif position > len(value) - 35 and separate_cterm:
                        cterm_seq += codon
                        cterm_class = True
                        if len(cterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["cterminus"]},columns=["Gene","Class"]))

                    elif pull_transitions:
                        for i in xrange(len(all_loc)):
                            transit_class = findTransitions(position,all_loc[i])
                            if transit_class:
                                seq[i]+=codon
                                if len(seq[i]) == 3:
                                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[prefixes[i]]},columns=["Gene","Class"]))
                                break
                    if (not transit_class) and (not nterm_class) and (not cterm_class):
                        if ss[position-offset] == "H":
                            helix_seq += codon
                            if len(helix_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "C" or ss[position - offset]=="B" or ss[position - offset]=="G" or ss[position - offset]=="I":
                            coil_seq += codon
                            if len(coil_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "E":
                            sheet_seq += codon
                            if len(sheet_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "T" or ss[position - offset] == "S":
                            turn_seq += codon
                            if len(turn_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["turn"]},columns=["Gene","Class"]))
             
                writeToFile(seq[0],cds_id,startc)
                writeToFile(seq[1],cds_id,endc)
                writeToFile(seq[2],cds_id,starth)
                writeToFile(seq[3],cds_id,endh)
                writeToFile(seq[4],cds_id,starte)
                writeToFile(seq[5],cds_id,ende)
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(coil_seq,cds_id,coil)
                writeToFile(sheet_seq,cds_id,sheet)
                writeToFile(turn_seq,cds_id,turn)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                writeToFile(cterm_seq,cds_id,cterm)
                assert len(comp_seq) == (len(cterm_seq) + len(turn_seq) + len(nterm_seq) + len(seq[0]) + len(seq[1]) + len(seq[2]) + len(seq[3]) + len(seq[4]) + len(seq[5]) +len(helix_seq) + len(sheet_seq) + len(coil_seq))
                assert len(comp_seq)/3 == len(keep)
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used

def pullSecondaryDownstream(output_folder,genome_file = "mod_scerevisiae.fasta",ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,ds_helix=True,ds_sheet=True,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil"):
        os.makedirs(output_folder+"/Coil")
    if not os.path.exists(output_folder+"/Helix"):
        os.makedirs(output_folder+"/Helix")
    if not os.path.exists(output_folder+"/Sheet"):
        os.makedirs(output_folder+"/Sheet")
    if not os.path.exists(output_folder+"/Turn"):
        os.makedirs(output_folder+"/Turn")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    if not os.path.exists(output_folder+"/Downstream_helix"):
        os.makedirs(output_folder+"/Downstream_helix")
    if not os.path.exists(output_folder+"/Downstream_sheet"):
        os.makedirs(output_folder+"/Downstream_sheet")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
       # blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Downstream_helix/downstream_helix.fasta",'w') as dsh,\
    open(output_folder+"/Downstream_sheet/downstream_sheet.fasta",'w') as dse,\
    open(output_folder+"/Coil/coil.fasta",'w') as coil,\
    open(output_folder+"/Sheet/sheet.fasta",'w') as sheet,\
    open(output_folder+"/Helix/helix.fasta",'w') as helix,\
    open(output_folder+"/Turn/turn.fasta",'w') as turn,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:

            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except IOError:
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1 
                    sstart = 0
                    send = len(value) - 1
                dsh_seq=''
                dse_seq=''
                coil_seq = ''
                helix_seq = ''
                sheet_seq = ''
                turn_seq = ''
                comp_seq = ''
                nterm_seq = ''
                
                h_pat = re.compile(r"HHHHH+")
                e_pat = re.compile(r"EEEEE+")
                ss = "".join(df.iloc[:,2].values)
                
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                dsh_transitions = h_pat.finditer(ss)
                dse_transitions = e_pat.finditer(ss)

                dsh_locations = []
                dse_locations = []
                if ds_helix:
                    for i in dsh_transitions:
                        dsh_locations.append((i.end()-1+30+offset,i.end()-1+35+offset))
                if ds_sheet:
                    for i in dse_transitions:
                        dse_locations.append((i.end()-1+30+offset,i.end()-1+35+offset))
             
                all_loc = [dsh_locations,dse_locations]
                
                prefixes = ["downstream_helix","downstream_sheet"]
                
                seq = [dsh_seq,dse_seq]
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if ss[position- offset] == "X":
                        continue
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    transit_class = False
                    nterm_class = False
                    
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))

                    elif pull_transitions:
                        for i in xrange(len(all_loc)):
                            transit_class = findTransitions(position,all_loc[i])
                            if transit_class:
                                seq[i]+=codon
                                if len(seq[i]) == 3:
                                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[prefixes[i]]},columns=["Gene","Class"]))
                                break
                    if (not transit_class) and (not nterm_class):
                        if ss[position - offset] == "H":
                            helix_seq += codon
                            if len(helix_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "C" or ss[position- offset]=="B" or ss[position- offset]=="G" or ss[position- offset]=="I":
                            coil_seq += codon
                            if len(coil_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "E":
                            sheet_seq += codon
                            if len(sheet_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "T" or ss[position- offset] == "S":
                            turn_seq += codon
                            if len(turn_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["turn"]},columns=["Gene","Class"]))
                writeToFile(seq[0],cds_id,dsh)
                writeToFile(seq[1],cds_id,dse)
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(coil_seq,cds_id,coil)
                writeToFile(sheet_seq,cds_id,sheet)
                writeToFile(turn_seq,cds_id,turn)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)

                assert len(comp_seq) == (len(turn_seq) + len(nterm_seq) + len(seq[0]) + len(seq[1]) +len(helix_seq) + len(sheet_seq) + len(coil_seq))
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used


def pullDisordered(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False,residues_to_keep = None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Ordered"):
        os.makedirs(output_folder+"/Ordered")
    if not os.path.exists(output_folder+"/Disordered"):
        os.makedirs(output_folder+"/Disordered")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
       # blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Ordered/ordered.fasta",'w') as ordered,\
    open(output_folder+"/Disordered/disordered.fasta",'w') as disordered,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    print "check not only pulling chap"
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())   
                try:
                    df = pd.read_csv(ss2_loc+prot_id,sep="\t",header=None)
                    if residues_to_keep:
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    #print "Could not find ",residues_to_keep+prot_id
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1 
                    sstart = 0
                    send = len(value) - 1
                ordered_seq = ''
                disordered_seq = ''
                comp_seq = ''
                nterm_seq = ''
                df.loc[df.iloc[:,3] < 0.5,2] = "O"
                df.loc[df.iloc[:,3] >= 0.5,2] = "D"
                ss = "".join(df.iloc[:,2].values)

                ## Use qstart:qend here because structure files are based on s.cerevisiae sequences, not pdb sequences...
                ## probably want to use pdb sequences in future
                
                length = len(ss[qstart:qend])
             
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if residues_to_keep != None:
                        if position not in keep:
                            continue
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    nterm_class = False
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    if not nterm_class:
                        if ss[position] == "O":
                            ordered_seq += codon
                            if len(ordered_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["ordered"]},columns=["Gene","Class"]))
                        elif ss[position] == "D":
                            disordered_seq += codon
                            if len(disordered_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["disordered"]},columns=["Gene","Class"]))
                writeToFile(ordered_seq,cds_id,ordered)
                writeToFile(disordered_seq,cds_id,disordered)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                assert len(comp_seq) == (len(ordered_seq) + len(nterm_seq) + len(disordered_seq))
                assert len(comp_seq)/3 == len(keep)
    
    tmp = pd.concat(tmp,ignore_index=True)

    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used


def pullDisorderedDownstream(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Downstream_ordered"):
        os.makedirs(output_folder+"/Downstream_ordered")
    if not os.path.exists(output_folder+"/Downstream_disordered"):
        os.makedirs(output_folder+"/Downstream_disordered")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
       # blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Downstream_ordered/ordered.fasta",'w') as ordered,\
    open(output_folder+"/Downstream_disordered/disordered.fasta",'w') as disordered,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    print "check not only pulling chap"
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())   
                try:
                    df = pd.read_csv(ss2_loc+prot_id,sep="\t",header=None)
                except (IOError,pd.errors.EmptyDataError):
                    #print "Could not find ",residues_to_keep+prot_id
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1 
                    sstart = 0
                    send = len(value) - 1
                ordered_seq = ''
                disordered_seq = ''
                comp_seq = ''
                nterm_seq = ''
                df.loc[df.iloc[:,3] < 0.5,2] = "O"
                df.loc[df.iloc[:,3] >= 0.5,2] = "D"
                ss = "".join(df.iloc[:,2].values)

                ## Use qstart:qend here because structure files are based on s.cerevisiae sequences, not pdb sequences...
                ## probably want to use pdb sequences in future
                
                length = len(ss[qstart:qend])
             
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    nterm_class = False
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    if not nterm_class:
                        if ss[position - 35] == "O":
                            ordered_seq += codon
                            if len(ordered_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["ordered"]},columns=["Gene","Class"]))
                        elif ss[position-35] == "D":
                            disordered_seq += codon
                            if len(disordered_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["disordered"]},columns=["Gene","Class"]))
                writeToFile(ordered_seq,cds_id,ordered)
                writeToFile(disordered_seq,cds_id,disordered)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                assert len(comp_seq) == (len(ordered_seq) + len(nterm_seq) + len(disordered_seq))
       
    
    tmp = pd.concat(tmp,ignore_index=True)

    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used




def pullSecondaryStartHelix(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil"):
        os.makedirs(output_folder+"/Coil")
    if not os.path.exists(output_folder+"/Helix"):
        os.makedirs(output_folder+"/Helix")
    if not os.path.exists(output_folder+"/Sheet"):
        os.makedirs(output_folder+"/Sheet")
    if not os.path.exists(output_folder+"/Turn"):
        os.makedirs(output_folder+"/Turn")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    if not os.path.exists(output_folder+"/Start_helix"):
        os.makedirs(output_folder+"/Start_helix")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
       # blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Start_helix/start_helix.fasta",'w') as starth,\
    open(output_folder+"/Coil/coil.fasta",'w') as coil,\
    open(output_folder+"/Sheet/sheet.fasta",'w') as sheet,\
    open(output_folder+"/Helix/helix.fasta",'w') as helix,\
    open(output_folder+"/Turn/turn.fasta",'w') as turn,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:

            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    print "check not only pulling chap"
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())   
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",delim_whitespace=True,skiprows=num_rows,header=None)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except (IOError,pd.errors.EmptyDataError):
                    #print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1 
                    sstart = 0
                    send = len(value) - 1
                starth_seq = ''
                coil_seq = ''
                helix_seq = ''
                sheet_seq = ''
                turn_seq = ''
                comp_seq = ''
                nterm_seq = ''
                starth_pat = re.compile(r"HHHH+")
                ss = "".join(df.iloc[:,2].values)
                
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                starth_transitions = starth_pat.finditer(ss)
                offset = qstart - sstart
                starth_locations = []
                for i in starth_transitions:
                    starth_locations.append((i.start()+1+offset,i.start()+3+offset))
                all_loc = [starth_locations]  
                prefixes = ["start_helix"]
                in_coil = set([])
                seq = [starth_seq]
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if ss[position- offset] == "X":
                        continue
                    comp_seq+=codon
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    transit_class = False
                    nterm_class = False
                    
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))

                    elif pull_transitions:
                        for i in xrange(len(all_loc)):
                            transit_class = findTransitions(position,all_loc[i])
                            if transit_class:
                                seq[i]+=codon
                                if len(seq[i]) == 3:
                                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[prefixes[i]]},columns=["Gene","Class"]))
                                break
                    if (not transit_class) and (not nterm_class):
                        if ss[position - offset] == "H":
                            helix_seq += codon
                            if len(helix_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "C" or ss[position- offset]=="B" or ss[position- offset]=="G" or ss[position- offset]=="I":
                            coil_seq += codon
                            if len(coil_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "E":
                            sheet_seq += codon
                            if len(sheet_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet"]},columns=["Gene","Class"]))
                        elif ss[position- offset] == "T" or ss[position- offset] == "S":
                            turn_seq += codon
                            if len(turn_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["turn"]},columns=["Gene","Class"]))
    
                writeToFile(seq[0],cds_id,starth)
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(coil_seq,cds_id,coil)
                writeToFile(sheet_seq,cds_id,sheet)
                writeToFile(turn_seq,cds_id,turn)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                assert len(comp_seq) == (len(turn_seq) + len(nterm_seq) + len(seq[0]) +len(helix_seq) + len(sheet_seq) + len(coil_seq))
    
    tmp = pd.concat(tmp,ignore_index=True)

    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used


def pullCompleteSeq(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_chap=False,pull_transitions=False,offset = 0,genes_used = None,only_chap= False,exclude_tm_sp = False,exclude_mito=False,separate_nterm=False,nterm_region=35,experimental=False,remove_disorder_protein=False,remove_disorder_region=False,use_emp_seq=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Complete_seq"):
        os.makedirs(output_folder+"/Complete_seq")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_chap=exclude_chap,exclude_mito=exclude_mito,only_chap=only_chap,exclude_tm_sp=exclude_tm_sp)
    
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental or use_emp_seq:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
       # blast_results = blast_results.loc[(blast_results.loc[:,"Percent_identity"] >= 95.0) & (blast_results.loc[:,"Gaps"] == 0) & (blast_results.loc[:,"Query_length"] == blast_results.loc[:,"Subject_length"]),:]
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm,\
    open(output_folder+"/Complete_seq/complete_seq.fasta",'w') as comp:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
                if experimental or use_emp_seq:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_chap or exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                if only_chap:
                    if prot_id not in prot_to_include:
                        continue
                value = re.findall(pat,line.strip())
                
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=None,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_protein:
                        if len(df.loc[df.iloc[:,5] == "X",2].values)/float(len(df.iloc[:,2].values)) > 0.1:
                            continue
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                except IOError:
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1 
                    sstart = 0
                    send = len(value) - 1
                ss = "".join(df.iloc[:,2].values)
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                comp_seq = ''
                nterm_seq = ''
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if ss[position - offset] == "X":
                        continue
                    if position < nterm_region and separate_nterm:
                        nterm_seq += codon
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    else:
                        comp_seq+=codon
                        if len(comp_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used



def getProteinIdGeneName():
    pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    gene_name = re.compile(r"gene=([A-Z0-9\-_']{4,10})")
    with open("mod_scerevisiae.fasta") as fin,open("protein_gene_mapping.csv",'w') as out:
        out.write("Protein,Gene\n")
        for line in fin:
            if line[0] == ">":
                protein = pat_id.search(line)
                gene = gene_name.search(line)
                if protein != None and gene != None:
                    protein = protein.group(1)
                    gene = gene.group(1)
                    out.write(",".join([protein,gene]))
                    out.write("\n")

def getChaperoneAssociated():
    id_map = pd.read_csv("protein_gene_mapping.csv",header=0,index_col=1)
    proteins = {}
    with open("chaperone_interactions.tsv") as fin,open("chaperone_associated.tsv",'w') as out:
        fin.readline()
        for line in fin:
            line_spt = line.split()
            for gene in line_spt:
                gene_name = gene.split(",")[0].upper()
                try:
                    protein_name =id_map.loc[gene_name,"Protein"]
                    print protein_name
                    current = proteins.get(protein_name)
                    if current != None:
                        proteins[protein_name] = current + 1
                    else:
                        proteins[protein_name] = 1
                except KeyError:
                    continue
        for k in proteins.keys():
            chap_count = str(proteins.get(k))
            out.write("\t".join([k,chap_count]))
            out.write("\n")
    return None


def pullPredScer(location,ss2_loc="Alpha_psipred/"):

    # genes_used =  pullSecondaryStartHelix(location+"/Secondary_structures_pechmann/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures_pechmann/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")
    
    # genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures_emp_seq/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,use_emp_seq=True)
    # createExpressionFiles(location+"Secondary_structures_emp_seq/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures_unconserved/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,use_emp_seq=False,residues_to_keep="Conserved_sites/Unconserved_positions/")
    # createExpressionFiles(location+"Secondary_structures_unconserved/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")


    # genes_used =  pullDisordered(location+"Ordered_disordered/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Ordered_disordered/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    genes_used =  pullDisordered(location+"Ordered_disordered_conserved_sacch/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,residues_to_keep="Conserved_sites/Conserved_positions_sacch/")
    createExpressionFiles(location+"Ordered_disordered_conserved_sacch/",genes_used)

    genes_used =  pullDisordered(location+"Ordered_disordered_unconserved_sacch/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,residues_to_keep="Conserved_sites/Unconserved_positions_sacch/")
    createExpressionFiles(location+"Ordered_disordered_unconserved_sacch/",genes_used)
    
    # genes_used =  pullSecondaryOrder(location+"Secondary_structure_order/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structure_order/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")
     
    # genes_used =  pullDisorderedDownstream(location+"Downstream_ordered_disordered/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Downstream_ordered_disordered/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")
    

    # genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures_begin_end/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures_begin_end/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_structured/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Downstream_structured/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_helix/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,ds_sheet=False,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Downstream_helix/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_sheet/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,ds_helix=False,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Downstream_sheet/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")


    # genes_used =  pullCompleteSeq(location+"Complete_seq/",ss2_loc=ss2_loc,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=False,experimental=False)
    # createExpressionFiles(location+"Complete_seq/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullCompleteSeq(location+"Nterm/",ss2_loc=ss2_loc,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Nterm/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # groups = ["CC","CH","CE","HH","HC","HE","N","EE","EC","EH"]
    # #groups = ["CC","CH","CE","HH","HC","HE","N","EE","EC","EH"]
    # genes_used = None
    # for i in groups:
    #     genes_used = structureDownstreamCat_simplified(location+"Structure_Downstream_Combo/",ss2_loc=ss2_loc,genes_used=genes_used,struct=i,exclude_mito=True,only_chap=False,exclude_tm_sp=False,experimental=False)
    # createExpressionFiles(location+"Structure_Downstream_Combo/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

def pullEmpScer(location,ss2_loc="Alpha_psipred/",genome_file="mod_scerevisiae.fasta",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",remove_disorder_region=False):
    # genes_used =  pullSecondaryStartHelix(location+"/Secondary_structures_pechmann/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_pechmann/",genes_used)

    genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    createExpressionFiles(location+"Secondary_structures/",genes_used)

    # genes_used =  pullDisordered(location+"Ordered_disordered/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Ordered_disordered/",genes_used)


    # genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures_begin_end/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_begin_end/",genes_used)

    # genes_used =  pullSecondaryDownstream(location+"Downstream_structured/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_structured/",genes_used)

    # genes_used =  pullSecondaryDownstream(location+"Downstream_helix/",genome_file=genome_file,ss2_loc=ss2_loc,blast_file=blast_file,pull_transitions=True,exclude_mito=True,ds_sheet=False,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_helix/",genes_used)

    # genes_used =  pullSecondaryDownstream(location+"Downstream_sheet/",genome_file=genome_file,ss2_loc=ss2_loc,blast_file=blast_file,pull_transitions=True,exclude_mito=True,ds_helix=False,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_sheet/",genes_used)


    # genes_used =  pullCompleteSeq(location+"Complete_seq/",ss2_loc=ss2_loc,genome_file=genome_file,blast_file=blast_file,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=False,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Complete_seq/",genes_used)

    # genes_used =  pullCompleteSeq(location+"Nterm/",ss2_loc=ss2_loc,genome_file=genome_file,blast_file=blast_file,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Nterm/",genes_used)

    # groups = ["CC","CH","CE","CT","HH","HC","HE","HT","N","EE","EC","EH","ET","TT","TC","TE","TH"]
    # #groups = ["CC","CH","CE","HH","HC","HE","N","EE","EC","EH"]
    # genes_used = None
    # for i in groups:
    #     genes_used = structureDownstreamCat_simplified(location+"Structure_Downstream_Combo/",blast_file=blast_file,genome_file=genome_file,ss2_loc=ss2_loc,genes_used=genes_used,struct=i,exclude_mito=True,only_chap=False,exclude_tm_sp=False,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Structure_Downstream_Combo/",genes_used)

def pullEmpEcoli(location,ss2_loc="Alpha_psipred/",genome_file="Ecoli_K12_MG1655_main.fasta",remove_disorder_region=False,blast_file="Exp_structure_data/ecoli_conservative_pid_95_80_coverage.blast"):
    # genes_used =  pullSecondaryStartHelix(location+"/Secondary_structures_pechmann/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_pechmann/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures/",genome_file=genome_file,ss2_loc=ss2_loc,blast_file=blast_file,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    createExpressionFiles(location+"Secondary_structures/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")


    # genes_used =  pullSecondaryBeginEnd(location+"Secondary_structures_begin_end/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_begin_end/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_structured/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_structured/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_helix/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,ds_sheet=False,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_helix/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryDownstream(location+"Downstream_sheet/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,ds_helix=False,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Downstream_sheet/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")


    # genes_used =  pullCompleteSeq(location+"Complete_seq/",ss2_loc=ss2_loc,genome_file=genome_file,blast_file=blast_file,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=False,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Complete_seq/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullCompleteSeq(location+"Nterm/",ss2_loc=ss2_loc,genome_file=genome_file,blast_file=blast_file,pull_transitions=False,offset = 0,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Nterm/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")

    # groups = ["CC","CH","CE","CT","HH","HC","HE","HT","N","EE","EC","EH","ET","TT","TC","TE","TH"]
    # #groups = ["CC","CH","CE","HH","HC","HE","N","EE","EC","EH"]
    # genes_used = None
    # for i in groups:
    #     genes_used = structureDownstreamCat_simplified(location+"Structure_Downstream_Combo/",blast_file=blast_file,genome_file=genome_file,ss2_loc=ss2_loc,genes_used=genes_used,struct=i,exclude_mito=True,only_chap=False,exclude_tm_sp=False,experimental=True,remove_disorder_region=remove_disorder_region)
    # createExpressionFiles(location+"Structure_Downstream_Combo/",genes_used,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")



#pullEmpNoChapAssociated()
#pullEmpChapAssociated()
#pullEmpEcoli("Ecoli/Exp_liberal_homology/",ss2_loc = "Exp_structure_data/Ecoli_psipred_format_liberal/",blast_file="Exp_structure_data/ecoli_liberal_pid_80_60_coverage.blast",remove_disorder_region=False)
#pullEmpScer("Scer/Exp_conservative_homology_missing_data/",genome_file="mod_scerevisiae.fasta",ss2_loc = "Exp_structure_data/Scer_psipred_format_conservative/",remove_disorder_region=True,blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast")
pullPredScer("Scer/Predicted/",ss2_loc = "IUPRED2/IUPRED_Results/")
# for i in xrange(1,11):
#     genome = Genome()

#     genome.readFasta("Evaluate_MH_test/simulated_scer_exp_conservative_homology_seq_"+str(i)+".fasta")
#     genome.writeFasta("Evaluate_MH_test/mod_simulated_scer_exp_conservative_homology_seq_"+str(i)+".fasta")
#structuresRepresented("scer_protein_structure_percentage.tsv",genome_file="mod_scerevisiae.fasta",ss2_loc="Exp_structure_data/Scer_psipred_format_conservative/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",exclude_mito=True,experimental=True)
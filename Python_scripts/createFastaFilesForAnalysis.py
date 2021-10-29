## Author: Alexander Cope
## Code for parsing secondary structure files from PsiPred or from data downloaded from PDB.
import re
import pandas as pd
import numpy as np
import os
import sys
import glob
              

def mergePsiPredIUPRED():
    mapping = {"C":{"O":0,"D":0},"H":{"O":0,"D":0},"E":{"O":0,"D":0}}
    pat = re.compile("C+")
    order = glob.glob("../Data/IUPRED2_ecoli/*")
    for f in order:
        try:
            prot = f.split("/")[-1]
            ss2 = "../Data/Psipred_ecoli/"+prot+".ss2"
            df_ss = pd.read_csv(ss2,header=None,delim_whitespace=True,skiprows=2)
            df_ss.columns = ["Position","AA","SS","C","H","E"]
            df_order = pd.read_csv(f,header=None,sep="\t")
            assert len(df_ss.index) == len(df_order.index)
            scores = df_order.iloc[:,2].values
            anchor = df_order.iloc[:,3].values
            df_order["Order_Score"] = scores
            df_order.loc[df_order.loc[:,"Order_Score"] < 0.5,2] = "O"
            df_order.loc[df_order.loc[:,"Order_Score"] >= 0.5,2] = "D"
            ss = df_ss.iloc[:,2].values
            ss_string = "".join(ss)
            order = df_order.iloc[:,2].values
            df_ss = pd.concat([df_ss,pd.DataFrame({"Order":order,"Order_Score":scores})],axis=1)
            order_string = "".join(order)
            df_ss.to_csv("../Data/Psipred_iupred_ecoli/"+prot+".ss2",sep="\t",header=False,index=False,quoting=False)
        except IOError:
            continue

def writeToFile(seq,prot_id,file_var):
    if seq:
        file_var.write(prot_id+"\n")
        file_var.write(seq+"\n")



def createExpressionFiles(output_folder,genes_used,gene_expression_file="../mod_scerevisiae_expression_wo_phi_allunique.csv"):
    files = glob.glob(output_folder+"/*/*.fasta")
    phi = pd.read_csv(gene_expression_file,header=0,index_col=0)
    for f in files:
        split_dir = f.split("/")
        subfolder = split_dir[-2]
        type_struct = split_dir[-1].split(".fasta")[0]
        struct_genes = genes_used.loc[genes_used.iloc[:,1] == type_struct,"Gene"]
        phi_struct = phi.loc[struct_genes,:]
        phi_struct.to_csv(output_folder+"/"+subfolder+"/"+type_struct+"_phi.csv",sep=",",header=True,index=True,quoting=False)
    return None



def findTransitions(position,locations):
    for l in locations:
        if l == position:
            return True
    return False



def getProtsForAnalysis(exclude_mito=True):
    prot_to_exclude = set([])
    prot_to_include = set([])
    if exclude_mito:
        mito = pd.read_csv("../mitochodrial_protein_exclude.txt",sep="\n",header=-1)
        prot_to_exclude = prot_to_exclude.union(list(mito.iloc[:,0].values))
    return prot_to_exclude,prot_to_include


def pullSecondaryOrder(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",genes_used = None,exclude_mito=True,separate_nterm=True,residues_to_keep=None):
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
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
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
            if line[0] != ">" and qstart !=-1:
                if exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=None,delim_whitespace=True)
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
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
                        if ss[position] == "H" and order[position] == "O":
                            helixo_seq += codon
                            if len(helixo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_ord"]},columns=["Gene","Class"]))
                        elif ss[position] == "H" and order[position] == "D":
                            helixd_seq += codon
                            if len(helixd_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_dis"]},columns=["Gene","Class"]))
                        elif ss[position] == "C" and order[position]=="O":
                            coilo_seq += codon
                            if len(coilo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_ord"]},columns=["Gene","Class"]))
                        elif ss[position] == "C" and order[position] == "D":
                            coild_seq += codon
                            if len(coild_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_dis"]},columns=["Gene","Class"]))
                        elif ss[position] == "E" and order[position] == "O":
                            sheeto_seq += codon
                            if len(sheeto_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet_ord"]},columns=["Gene","Class"]))
                        elif ss[position] == "E" and order[position] == "D":
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


def pullSecondaryPaired(output_folder,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",genes_used = None,exclude_mito=True,separate_nterm=True,residues_to_keep=None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil_paired"):
        os.makedirs(output_folder+"/Coil_paired")
    if not os.path.exists(output_folder+"/Helix_paired"):
        os.makedirs(output_folder+"/Helix_paired")
    if not os.path.exists(output_folder+"/Sheet_paired"):
        os.makedirs(output_folder+"/Sheet_paired")
    if not os.path.exists(output_folder+"/Coil_unpaired"):
        os.makedirs(output_folder+"/Coil_unpaired")
    if not os.path.exists(output_folder+"/Helix_unpaired"):
        os.makedirs(output_folder+"/Helix_unpaired")
    if not os.path.exists(output_folder+"/Sheet_unpaired"):
        os.makedirs(output_folder+"/Sheet_unpaired")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    pat = re.compile(r"[ATCGYNW]{3}")
    #pat_id = re.compile(r"protein_id=(NP_[0-9]+\.[0-9])")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Coil_paired/coil_paired.fasta",'w') as coil_ord,\
    open(output_folder+"/Sheet_paired/sheet_paired.fasta",'w') as sheet_ord,\
    open(output_folder+"/Helix_paired/helix_paired.fasta",'w') as helix_ord,\
    open(output_folder+"/Coil_unpaired/coil_unpaired.fasta",'w') as coil_dis,\
    open(output_folder+"/Sheet_unpaired/sheet_unpaired.fasta",'w') as sheet_dis,\
    open(output_folder+"/Helix_unpaired/helix_unpaired.fasta",'w') as helix_dis,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
        for line in fin:
            post_dict = {}
            if line[0] == ">":
                prot_id = pat_id.search(line).group(1)
                cds_id = line.split()[0]
            if line[0] != ">" and qstart !=-1:
                if exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=None,delim_whitespace=True)
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
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
                order = df.iloc[:,8].values
                aa = "".join(df.iloc[:,1].values)
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    comp_seq+=codon
                    aa_pos = aa[position]
                    # if aa_pos == "L":
                    # 	paired_pat = re.compile(r"([\(\)]..)|(..[\(\)])")
                    # 	unpaired_pat = re.compile(r"\..\.")
                    # else:
                    # 	paired_pat = re.compile(r"..[\(\)]")
                    # 	unpaired_pat = re.compile(r"..\.")
                    paired_pat = re.compile(r"([\(\)]..)|(..[\(\)])|(.[\(\)].)")
                    unpaired_pat = re.compile(r"\.\.\.")
                    if len(comp_seq) == 3:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["complete_seq"]},columns=["Gene","Class"]))
                    nterm_class = False
                    if position < 35 and separate_nterm:
                        nterm_seq += codon
                        nterm_class = True
                        if len(nterm_seq) == 3:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                    if not nterm_class:
                        if ss[position] == "H" and (paired_pat.match(order[position]) != None):
                            helixo_seq += codon
                            if len(helixo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_paired"]},columns=["Gene","Class"]))
                        elif ss[position] == "H" and (unpaired_pat.match(order[position]) != None):
                            helixd_seq += codon
                            if len(helixd_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix_unpaired"]},columns=["Gene","Class"]))
                        elif ss[position] == "C" and (paired_pat.match(order[position]) != None):
                            coilo_seq += codon
                            if len(coilo_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_paired"]},columns=["Gene","Class"]))
                        elif ss[position] == "C" and (unpaired_pat.match(order[position]) != None):
                            coild_seq += codon
                            if len(coild_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil_unpaired"]},columns=["Gene","Class"]))
                        elif ss[position] == "E" and (paired_pat.match(order[position]) != None):
                            sheeto_seq += codon
                            if len(sheeto_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet_paired"]},columns=["Gene","Class"]))
                        elif ss[position] == "E" and (unpaired_pat.match(order[position]) != None):
                            sheetd_seq += codon
                            if len(sheetd_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet_unpaired"]},columns=["Gene","Class"]))
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


def pullSecondaryStructures(output_folder,minimum_length=0,maximum_length=10000000000000000,termini_size=2,genome_file = "../Data/Fasta/mod_scerevisiae.fasta" ,ss2_loc = "../Data/Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",pull_transitions=False,genes_used = None,exclude_mito=True,separate_nterm=True,experimental=False,remove_disorder_region=False,residues_to_keep=None):
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
    if not os.path.exists(output_folder+"/Core_sheet"):
        os.makedirs(output_folder+"/Core_sheet")
    if not os.path.exists(output_folder+"/Core_helix"):
        os.makedirs(output_folder+"/Core_helix")
    if not os.path.exists(output_folder+"/Core_coil"):
        os.makedirs(output_folder+"/Core_coil")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"protein_id=([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    end_gi = 0.0
    core_gi = 0.0
    start_gi = 0.0
    total_start =0.0
    total_end = 0.0
    total_core =0.0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/End_coil/end_coil.fasta",'w') as endc,\
    open(output_folder+"/Start_coil/start_coil.fasta",'w') as startc,\
    open(output_folder+"/Start_helix/start_helix.fasta",'w') as starth,\
    open(output_folder+"/End_helix/end_helix.fasta",'w') as endh,\
    open(output_folder+"/Start_sheet/start_sheet.fasta",'w') as starte,\
    open(output_folder+"/End_sheet/end_sheet.fasta",'w') as ende,\
    open(output_folder+"/Core_sheet/core_sheet.fasta",'w') as coree,\
    open(output_folder+"/Core_helix/core_helix.fasta",'w') as coreh,\
    open(output_folder+"/Core_coil/core_coil.fasta",'w') as corec,\
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
                if experimental:
                    try:
                        qstart = int(blast_results.loc[prot_id,"Query_start"]) - 1
                        qend = int(blast_results.loc[prot_id,"Query_end"])
                        sstart = int(blast_results.loc[prot_id,"Subject_start"]) -1
                        send = int(blast_results.loc[prot_id,"Subject_end"])
                    except KeyError:
                        qstart=-1
            if line[0] != ">" and qstart !=-1:
                if exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                value = re.findall(pat,line.strip())
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                   # print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not experimental:
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
                startc_seq=''
                corec_seq = ''
                endc_seq=''
                starth_seq=''
                coreh_seq = ''
                endh_seq=''
                starte_seq=''
                coree_seq = ''
                ende_seq=''
                coil_seq = ''
                helix_seq = ''
                sheet_seq = ''
                turn_seq = ''
                comp_seq = ''
                nterm_seq = ''
                dis_seq = ''
                min_length = minimum_length

                h_pat = re.compile(r"[HGI]+")
                e_pat = re.compile(r"[EB]+")
                c_pat = re.compile(r"[CTS]+")               
                ss = "".join(df.iloc[:,2].values)

                gene_df = pd.DataFrame({"Codon":value[qstart:qend],
                                    "Structure":df.iloc[sstart:send,2].values,
                                    "Position": xrange(qstart,qend),
                                    "Relative.Position":[0] * (qend-qstart),
                                    "Length":[0] * (qend-qstart)},index=xrange(qstart,qend))
                h_transitions = h_pat.finditer(ss)
                e_transitions = e_pat.finditer(ss)
                c_transitions = c_pat.finditer(ss)
                offset = qstart - sstart
                for i in h_transitions:
                    offset_start = i.start() + offset
                    offset_end = i.end() + offset -1 
                    gene_df.loc[offset_start:offset_end,"Relative.Position"] = gene_df.loc[offset_start:offset_end,"Position"] - offset_start
                    gene_df.loc[offset_start:offset_end,"Length"] = offset_end - offset_start + 1
                for i in e_transitions:
                    offset_start = i.start() + offset
                    offset_end = i.end() + offset -1 
                    gene_df.loc[offset_start:offset_end,"Relative.Position"] = gene_df.loc[offset_start:offset_end,"Position"] - offset_start
                    gene_df.loc[offset_start:offset_end,"Length"] = offset_end - offset_start + 1 
                for i in c_transitions:
                    offset_start = i.start() + offset
                    offset_end = i.end() + offset -1
                    gene_df.loc[offset_start:offset_end,"Relative.Position"] = gene_df.loc[offset_start:offset_end,"Position"] - offset_start
                    gene_df.loc[offset_start:offset_end,"Length"] = offset_end - offset_start + 1
                


                gene_df = gene_df.loc[(gene_df.loc[:,"Structure"] != "X"),:]
                if separate_nterm:
                    nterm_seq = "".join(gene_df.loc[gene_df.loc[:,"Position"] < 35,"Codon"].values)
                    if len(nterm_seq) != 0:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                
                    gene_df = gene_df.loc[gene_df.loc[:,"Position"] >= 35,:]
                
                helical_struct = ["H","G","I"]
                sheet_struct = ["E","B"]
                coil_struct = ["T","S","C"]

                

                gene_df = gene_df.loc[(gene_df.loc[:,"Length"] >= min_length) & (gene_df.loc[:,"Length"] <= maximum_length),:]

                helix_seq = "".join(gene_df.loc[(gene_df["Structure"].isin(helical_struct)) ,"Codon"].values)
                if len(helix_seq) != 0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                
                sheet_seq = "".join(gene_df.loc[(gene_df["Structure"].isin(sheet_struct)),"Codon"].values)
                if len(sheet_seq) !=0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet"]},columns=["Gene","Class"]))

                coil_seq = "".join(gene_df.loc[(gene_df.loc[:,"Structure"] =="C"),"Codon"].values)
                if len(coil_seq) != 0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil"]},columns=["Gene","Class"]))
                
                turn_seq = "".join(gene_df.loc[(gene_df.loc[:,"Structure"] =="T") | (gene_df.loc[:,"Structure"] =="S"),"Codon"].values)
                if len(turn_seq) !=0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["turn"]},columns=["Gene","Class"]))


                struct = [helical_struct,sheet_struct,coil_struct]
                start = ["start_helix","start_sheet","start_coil"]
                core = ["core_helix","core_sheet","core_coil"]
                end = ["end_helix","end_sheet","end_coil"]
                
                start_file = [starth,starte,startc]
                core_file = [coreh,coree,corec]
                end_file = [endh,ende,endc]


                if pull_transitions:
                    for i in xrange(len(struct)):

                        struct_df = gene_df.loc[(gene_df["Structure"].isin(struct[i])),:]
                        
                        start_seq = "".join(struct_df.loc[(struct_df.loc[:,"Relative.Position"] < termini_size),"Codon"].values)
                        if len(start_seq) != 0:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[start[i]]},columns=["Gene","Class"]))
                        
                        end_seq = "".join(struct_df.loc[((struct_df.loc[:,"Length"] - struct_df.loc[:,"Relative.Position"]) <= termini_size),"Codon"].values)
                        if len(end_seq) != 0:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[end[i]]},columns=["Gene","Class"]))

                        core_seq = "".join(struct_df.loc[(struct_df.loc[:,"Relative.Position"] >= termini_size) & ((struct_df.loc[:,"Length"] - struct_df.loc[:,"Relative.Position"]) > termini_size),"Codon"].values)
                        if len(core_seq) != 0:
                            tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":[core[i]]},columns=["Gene","Class"]))
                    
                        writeToFile(start_seq,cds_id,start_file[i])
                        writeToFile(end_seq,cds_id,end_file[i])
                        writeToFile(core_seq,cds_id,core_file[i])
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(coil_seq,cds_id,coil)
                writeToFile(sheet_seq,cds_id,sheet)
                writeToFile(turn_seq,cds_id,turn)
                writeToFile(nterm_seq,cds_id,nterm)

                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                
               
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used


def pullSecondarySplitHelixTypes(output_folder,minimum_length=4,maximum_length=10000000000000000,genome_file = "../Data/Fasta/mod_scerevisiae.fasta" ,ss2_loc = "../Data/Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",pull_transitions=False,genes_used = None,exclude_tm_sp = False,exclude_mito=True,separate_nterm=False,experimental=False,remove_disorder_region=False,residues_to_keep=None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Coil"):
        os.makedirs(output_folder+"/Coil")
    if not os.path.exists(output_folder+"/Helix"):
        os.makedirs(output_folder+"/Helix")
    if not os.path.exists(output_folder+"/GI_Helix"):
        os.makedirs(output_folder+"/GI_Helix")
    if not os.path.exists(output_folder+"/Sheet"):
        os.makedirs(output_folder+"/Sheet")
    if not os.path.exists(output_folder+"/Turn"):
        os.makedirs(output_folder+"/Turn")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")

    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito,exclude_tm_sp=exclude_tm_sp)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    end_gi = 0.0
    core_gi = 0.0
    start_gi = 0.0
    total_start =0.0
    total_end = 0.0
    total_core =0.0
    with open(genome_file) as fin, \
    open(output_folder+"/complete_seq.fasta",'w') as comp,\
    open(output_folder+"/Coil/coil.fasta",'w') as coil,\
    open(output_folder+"/Sheet/sheet.fasta",'w') as sheet,\
    open(output_folder+"/Helix/helix.fasta",'w') as helix,\
    open(output_folder+"/GI_Helix/gi_helix.fasta",'w') as gi_helix,\
    open(output_folder+"/Turn/turn.fasta",'w') as turn,\
    open(output_folder+"/Nterminus/nterminus.fasta",'w') as nterm:
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
                if exclude_tm_sp or exclude_mito:
                    if prot_id in prot_to_exclude:
                        print prot_id
                        continue
                value = re.findall(pat,line.strip())
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",header=-1,delim_whitespace=True,skiprows=num_rows)
                    if experimental and remove_disorder_region:
                        df.loc[df.iloc[:,5]=="X",2] = "X"
                    if residues_to_keep != None:          
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                   # print "Could not find ",ss2_loc+prot_id+".ss2"
                    continue
                if not (experimental or use_emp_seq):
                    qstart = 0
                    qend = len(value) - 1
                    sstart = 0
                    send = len(value) - 1
               
                coil_seq = ''
                helix_seq = ''
                gi_helix_seq = ''
                sheet_seq = ''
                turn_seq = ''
                comp_seq = ''
                nterm_seq = ''
         
    

                ss = "".join(df.iloc[:,2].values)
                
                length = len(ss[sstart:send])
                assert length == len(value[qstart:qend])
                offset = qstart - sstart
                #offset = 0
   
                for position,codon in zip(xrange(qstart,qend),value[qstart:qend]):
                    if ss[position - offset] == "X":
                        continue
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
    
                    if (not nterm_class):
                        if ss[position-offset] == "H":
                            helix_seq += codon
                            if len(helix_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                        if ss[position-offset] == "G" or ss[position - offset]=="I":
                            gi_helix_seq += codon
                            if len(gi_helix_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["gi_helix"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "C":
                            coil_seq += codon
                            if len(coil_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["coil"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "E" or ss[position-offset] == "B":
                            sheet_seq += codon
                            if len(sheet_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["sheet"]},columns=["Gene","Class"]))
                        elif ss[position-offset] == "T" or ss[position - offset] == "S":
                            turn_seq += codon
                            if len(turn_seq) == 3:
                                tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["turn"]},columns=["Gene","Class"]))
             
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(gi_helix_seq,cds_id,gi_helix)
                writeToFile(coil_seq,cds_id,coil)
                writeToFile(sheet_seq,cds_id,sheet)
                writeToFile(turn_seq,cds_id,turn)
                writeToFile(comp_seq,cds_id,comp)
                writeToFile(nterm_seq,cds_id,nterm)
                assert len(comp_seq) == (len(turn_seq) + len(nterm_seq) +len(helix_seq) +len(gi_helix_seq) + len(sheet_seq) + len(coil_seq))
                # if residues_to_keep:
                #     assert len(comp_seq)/3 == len(keep)
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    #print core_gi/total_core,start_gi/total_start,end_gi/total_end
    return genes_used






## Expects merged PsiPred + IUPRED2 output
def pullDisordered(output_folder,genome_file = "../Data/Fasta/mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",genes_used = None,exclude_mito=True,separate_nterm=True,residues_to_keep = None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Ordered"):
        os.makedirs(output_folder+"/Ordered")
    if not os.path.exists(output_folder+"/Disordered"):
        os.makedirs(output_folder+"/Disordered")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito)
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    pat = re.compile(r"[ATCGYNW]{3}")
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
    
            if line[0] != ">" and qstart !=-1:
                if exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                value = re.findall(pat,line.strip())   
                try:
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",sep="\t",header=None)
                    if residues_to_keep:
                        keep = pd.read_csv(residues_to_keep+"/"+prot_id,squeeze=True)
                except (IOError,pd.errors.EmptyDataError):
                    #print "Could not find ",residues_to_keep+prot_id
                    continue
                qstart = 0
                qend = len(value) - 1 
                sstart = 0
                send = len(value) - 1
                ordered_seq = ''
                disordered_seq = ''
                comp_seq = ''
                nterm_seq = ''
                #df.loc[df.iloc[:,3] < 0.5,2] = "O"
                #df.loc[df.iloc[:,3] >= 0.5,2] = "D"
                ss = "".join(df.iloc[:,6].values)

                ## Use qstart:qend here because structure files are based on s.cerevisiae sequences, not pdb sequences...
                ## probably want to use pdb sequences in future
                
                length = len(ss[qstart:qend])
             
                assert length == len(value[qstart:qend])

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
                if residues_to_keep:
                    assert len(comp_seq)/3 == len(keep)
    
    tmp = pd.concat(tmp,ignore_index=True)

    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used





def pullSecondaryStartHelix(output_folder,minimum_length=4,maximum_length=10000000000000000,genome_file = "mod_scerevisiae.fasta" ,ss2_loc = "Alpha_psipred/",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",genes_used = None,exclude_mito=True,separate_nterm=True,experimental=True,remove_disorder_region=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder+"/Helix"):
        os.makedirs(output_folder+"/Helix")
    if not os.path.exists(output_folder+"/Nterminus"):
        os.makedirs(output_folder+"/Nterminus")
    if not os.path.exists(output_folder+"/Start_helix"):
        os.makedirs(output_folder+"/Start_helix")
    if not os.path.exists(output_folder+"/Remainder_helix"):
        os.makedirs(output_folder+"/Remainder_helix")
    prot_to_exclude,prot_to_include = getProtsForAnalysis(exclude_mito=exclude_mito)
    
    if genes_used is None:
        genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    blast_results=None
    if experimental:
        blast_results = pd.read_csv(blast_file,sep="\t",header=0,index_col=0)
    pat = re.compile(r"[ATCGYNW]{3}")
    pat_id = re.compile(r"([YN]P_[0-9]+\.[0-9])")
    qstart = 0
    with open(genome_file) as fin, \
    open(output_folder+"/Start_helix/start_helix.fasta",'w') as starth,\
    open(output_folder+"/Remainder_helix/remainder_helix.fasta",'w') as remh,\
    open(output_folder+"/Helix/helix.fasta",'w') as helix, \
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
                if exclude_mito:
                    if prot_id in prot_to_exclude:
                        continue
                value = re.findall(pat,line.strip())   
                try:
                    if experimental:
                        num_rows = 0
                    else:
                        num_rows = 2
                    df = pd.read_csv(ss2_loc+prot_id+".ss2",delim_whitespace=True,skiprows=num_rows,header=None)
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
                remh_seq = ''
                helix_seq = ''
                nterm_seq = ''
                starth_pat = re.compile(r"[HGI]+")
                ss = "".join(df.iloc[:,2].values)
                
                min_length = minimum_length

                gene_df = pd.DataFrame({"Codon":value[qstart:qend],
                                    "Structure":df.iloc[sstart:send,2].values,
                                    "Position": xrange(qstart,qend),
                                    "Relative.Position":[0] * (qend-qstart),
                                    "Length":[0] * (qend-qstart)},index=xrange(qstart,qend))
                starth_transitions = starth_pat.finditer(ss)
         
                offset = qstart - sstart
                for i in starth_transitions:
                    offset_start = i.start() + offset
                    offset_end = i.end() + offset -1 
                    gene_df.loc[offset_start:offset_end,"Relative.Position"] = gene_df.loc[offset_start:offset_end,"Position"] - offset_start
                    gene_df.loc[offset_start:offset_end,"Length"] = offset_end - offset_start + 1


                gene_df = gene_df.loc[(gene_df.loc[:,"Structure"] != "X"),:]
                if separate_nterm:
                    nterm_seq = "".join(gene_df.loc[gene_df.loc[:,"Position"] < 35,"Codon"].values)
                    if len(nterm_seq) != 0:
                        tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["nterminus"]},columns=["Gene","Class"]))
                
                    gene_df = gene_df.loc[gene_df.loc[:,"Position"] >= 35,:]
                
                helical_struct = ["H","G","I"]

                gene_df = gene_df.loc[(gene_df.loc[:,"Length"] >= min_length) & (gene_df.loc[:,"Length"] <= maximum_length),:]

                helix_seq = "".join(gene_df.loc[(gene_df["Structure"].isin(helical_struct)) ,"Codon"].values)
                if len(helix_seq) != 0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["helix"]},columns=["Gene","Class"]))
                
                
                struct_df = gene_df.loc[(gene_df["Structure"].isin(helical_struct)),:]
                relative_position = struct_df.loc[:,"Relative.Position"]
                start_seq = "".join(struct_df.loc[(relative_position == 1) | (relative_position == 2) ,"Codon"].values)
                if len(start_seq) != 0:
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["start_helix"]},columns=["Gene","Class"]))
                remainder_seq = "".join(struct_df.loc[(relative_position != 1) & (relative_position != 2),"Codon"].values)
                if len(remainder_seq) != 0: 
                    tmp.append(pd.DataFrame(data={"Gene":[cds_id[1:]],"Class":["remainder_helix"]},columns=["Gene","Class"]))

        
                writeToFile(start_seq,cds_id,starth)
                writeToFile(remainder_seq,cds_id,remh)
                writeToFile(helix_seq,cds_id,helix)
                writeToFile(nterm_seq,cds_id,nterm)  
    tmp = pd.concat(tmp,ignore_index=True)

    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used




def pullPred(location,ss2_loc="../Data/Psipred/",genome_file = "mod_scerevisiae.fasta"):

    # genes_used =  pullSecondaryStartHelix(location+"/Secondary_structures_pechmann/",ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures_pechmann/",genes_used)#,gene_expression_file="Ecoli_K12_MG1655_main_phi.csv")
    
    # genes_used =  pullSecondaryStructures(location+"Secondary_structures/",genome_file = genome_file,ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures/",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullSecondaryStructures(location+"Secondary_structures/",genome_file = genome_file,ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures/",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")
    
    # genes_used =  pullSecondaryStructures(location+"Secondary_structures",genome_file = genome_file,ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,separate_nterm=True,experimental=False)
    # createExpressionFiles(location+"Secondary_structures/",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")

    
    # genes_used =  pullDisordered(location+"Ordered_disordered/",genome_file = genome_file,ss2_loc=ss2_loc,exclude_mito=True,separate_nterm=True)
    # createExpressionFiles(location+"Ordered_disordered/",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")

    # genes_used =  pullDisordered(location+"Ordered_disordered_conserved_sacch_no_ncast/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,residues_to_keep="../Conserved_sites/Conserved_positions_sacch_no_ncast/")
    # createExpressionFiles(location+"Ordered_disordered_conserved_sacch_no_ncast/",genes_used)

    # genes_used =  pullDisordered(location+"Ordered_disordered_unconserved_sacch_no_ncast/",ss2_loc=ss2_loc,pull_transitions=False,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,experimental=False,residues_to_keep="../Conserved_sites/Unconserved_positions_sacch_no_ncast/")
    # createExpressionFiles(location+"Ordered_disordered_unconserved_sacch_no_ncast/",genes_used)
    
    genes_used =  pullSecondaryPaired(location+"Secondary_structure_paired_any_nucleotide/",genome_file = genome_file,ss2_loc=ss2_loc,exclude_mito=True,separate_nterm=True)
    createExpressionFiles(location+"Secondary_structure_paired_any_nucleotide/",genes_used,gene_expression_file="../mod_scerevisiae_expression_wo_phi_allunique.csv")
     

    # genes_used =  pullSecondaryOrder(location+"Secondary_structure_order/",genome_file = genome_file,ss2_loc=ss2_loc,exclude_mito=True,separate_nterm=True)
    # createExpressionFiles(location+"Secondary_structure_order/",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")
     
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

    
def pullEmpScer(location,ss2_loc="Alpha_psipred/",genome_file="mod_scerevisiae.fasta",blast_file="Exp_structure_data/scer_conservative_pid_95_80_coverage.blast",remove_disorder_region=False,gene_expression_file="../mod_scerevisiae_expression_wo_phi_allunique.csv"):
    #genes_used =  pullSecondaryStartHelix(location+"/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/",minimum_length=6,genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,exclude_mito=True,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    #createExpressionFiles(location+"Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/",genes_used,gene_expression_file=gene_expression_file)
    
    # genes_used =  pullSecondarySplitHelixTypes(location+"Secondary_structures_helix_types/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_helix_types",genes_used)

    # genes_used =  pullSecondaryBeginEndMin(location+"Secondary_structures_begin_end_length_4_2_codon_for_termini/",minimum_length=4,maximum_length=4,genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_begin_end_length_4_2_codon_for_termini",genes_used)


    for i in xrange(6,8):
        genes_used =  pullSecondaryStructures(location+"Secondary_structures_begin_end_length_at_least_"+str(i)+"_2_codon_for_termini/",minimum_length=i,termini_size=2,genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
        createExpressionFiles(location+"Secondary_structures_begin_end_length_at_least_"+str(i)+"_2_codon_for_termini",genes_used,gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")
        


    # genes_used =  pullSecondaryBeginEndMin(location+"Secondary_structures_begin_end_exclude_less_than_10_2_codon_for_termini/",genome_file=genome_file,blast_file=blast_file,ss2_loc=ss2_loc,pull_transitions=True,exclude_mito=True,exclude_tm_sp=False,separate_nterm=True,remove_disorder_region=remove_disorder_region,experimental=True)
    # createExpressionFiles(location+"Secondary_structures_begin_end_exclude_less_than_10_2_codon_for_termini/",genes_used)

#mergePsiPredIUPRED()
pullPred("../Scer/Predicted/",genome_file="../Data/Fasta/mod_scerevisiae.fasta",ss2_loc="../Data/Psipred_iupred_mrna/")
#pullEmp("../Ecoli/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/",genome_file="../Data/Fasta/Ecoli_K12_MG1655_main.fasta",ss2_loc = "../Data/PDB/Ecoli_psipred_format_conservative/",remove_disorder_region=True,blast_file="../Data/PDB/ecoli_conservative_pid_95_80_coverage.blast",gene_expression_file="../Ecoli_K12_MG1655_main_phi.csv")

# for i in xrange(1,101):
#     genome = Genome()
#     genome.readFasta("../Scer/Test_Fisher_Exact/"+str(i)+"/sim.fasta")
#     genome.writeFasta("../Scer/Test_Fisher_Exact/"+str(i)+"/mod_sim.fasta")
#     pullStartHelixForFisherExact("../Scer/Test_Fisher_Exact/"+str(i)+"/",genome_file="../Scer/Test_Fisher_Exact/"+str(i)+"/mod_sim.fasta",ss2_loc = "../Data/Exp_structure_data/Scer_psipred_format_conservative/")

#pullPredScer("../Scer/Predicted/",ss2_loc="../Data/Alpha_psipred/")

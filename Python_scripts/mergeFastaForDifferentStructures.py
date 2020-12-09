## Author: Alexander Cope
## Code for merging structural regions into same category. Currently setup for performing analysis examining within genome 
## Always intended to integrate this into parseStructures.py so these grouping would be created after the different fasta files
## were created, but never got around to it. 

import pandas as pd
import glob
import os
import re


def mergeGenomes(list_to_merge,base_dir="./"):
    seq = {}
    filenames = []
    new_folder = "_".join(list_to_merge)
    if not os.path.exists(base_dir+"/"+new_folder):
        os.makedirs(base_dir+"/"+new_folder)
    for f in list_to_merge:
        filename = base_dir+"/" + f
        filenames.append(f)
        file = filename+"/"+f.lower()+".fasta"
        with open(file) as fin:
            for line in fin:
                if line[0] == ">":
                    accession = line
                else:
                    check = seq.get(accession)
                    if check == None:
                        seq[accession] = line.strip()
                    else:
                        prev_len = len(check)
                        seq[accession] = check+line.strip()
    genes_used = pd.DataFrame(data={"Gene":[],"Class":[]},columns=["Gene","Class"])
    tmp = []
    group = new_folder.lower()
    with open(base_dir+"/"+new_folder+"/"+group+".fasta",'w') as out:
        for key in seq.keys():
            out.write(key)
            out.write(seq.get(key))
            out.write("\n")
            tmp.append(pd.DataFrame(data={"Gene":[key[1:].strip()],"Class":[group]},columns=["Gene","Class"]))
    tmp = pd.concat(tmp,ignore_index=True)
    genes_used = pd.concat([genes_used,tmp],ignore_index=True)
    return genes_used,new_folder


def createExpressionFiles(output_folder,genes_used):
    #print output_folder
    files = glob.glob(output_folder+"/*.fasta")
    #print files
    phi = pd.read_csv("../mod_scerevisiae_expression_wo_phi_allunique.csv",header=0,index_col=0)
    for f in files:
        type_struct = f.split("/")[-1].split(".fasta")[0]
        struct_genes = genes_used.loc[genes_used.iloc[:,1] == type_struct,"Gene"]
        phi_struct = phi.loc[struct_genes,:]
        phi_struct.to_csv(output_folder+"/"+type_struct+"_phi.csv",sep=",",header=True,index=True,quoting=False)
    return None



genes_used,new_folder = mergeGenomes(["Turn","Coil"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)


genes_used,new_folder = mergeGenomes(["Start_coil","End_coil"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

genes_used,new_folder = mergeGenomes(["Start_sheet","End_sheet"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

genes_used,new_folder = mergeGenomes(["Start_helix","End_helix"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

genes_used,new_folder = mergeGenomes(["Start_coil","Core_coil","End_coil"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

genes_used,new_folder = mergeGenomes(["Start_sheet","Core_sheet","End_sheet"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

genes_used,new_folder = mergeGenomes(["Start_helix","Core_helix","End_helix"],base_dir="../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini")
createExpressionFiles("../Scer/Exp_conservative_homology_remove_X_effects_of_GI_in_H/Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini/"+new_folder,genes_used)

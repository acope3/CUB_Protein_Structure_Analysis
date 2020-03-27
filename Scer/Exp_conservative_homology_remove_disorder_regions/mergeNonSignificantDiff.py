import pandas as pd
import glob
from itertools import chain, combinations
import os
import re
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def mergeNonSigGenomes(list_to_merge,base_dir="./"):
    seq = {}
    filenames = []
    new_folder = "_".join(list_to_merge)
    if not os.path.exists(base_dir+"/"+new_folder):
        os.makedirs(base_dir+"/"+new_folder)
    for f in list_to_merge:
        #filename = f.split("/")[-2]
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
    phi = pd.read_csv("../../mod_scerevisiae_expression_wo_phi_allunique.csv",header=0,index_col=0)
    for f in files:
        type_struct = f.split("/")[-1].split(".fasta")[0]
        struct_genes = genes_used.loc[genes_used.iloc[:,1] == type_struct,"Gene"]
        phi_struct = phi.loc[struct_genes,:]
        phi_struct.to_csv(output_folder+"/"+type_struct+"_phi.csv",sep=",",header=True,index=True,quoting=False)
    return None


def splitSpan(string,sep="_",span=2):
    words = string.split(sep)
    return [sep.join(words[i:i+span]) for i in range(0,len(words),span)]


def testCheckMerge():
    assert checkCanMerge("Coil_coil_Helix_coil") == True
    assert checkCanMerge("Coil_coil_Coil_sheet") == True
    assert checkCanMerge("Helix_sheet_Sheet_coil") == False

    assert checkCanMerge("Helix_helix_Sheet_helix_Coil_helix") == True
    assert checkCanMerge("Sheet_helix_Sheet_sheet_Sheet_coil") == True
    assert checkCanMerge("Sheet_turn_Helix_sheet_Sheet_coil") == False
    assert checkCanMerge("Coil_helix_Coil_sheet_Coil_turn_Helix_coil") == False

    assert checkCanMerge("Helix_helix_Sheet_helix_Coil_helix_Turn_helix") == True
    assert checkCanMerge("Sheet_helix_Sheet_sheet_Sheet_coil_Sheet_turn") == True
    assert checkCanMerge("Sheet_sheet_Helix_sheet_Sheet_coil_Helix_coil") == True
    assert checkCanMerge("Sheet_turn_Helix_sheet_Sheet_coil_Helix_coil") == False
    assert checkCanMerge("Turn_turn_Helix_sheet_Sheet_coil_Helix_coil") == False

    assert checkCanMerge("Coil_helix_Coil_sheet_Coil_turn_Helix_helix_Helix_turn_Helix_sheet") == True
    assert checkCanMerge("Coil_helix_Coil_sheet_Helix_sheet_Helix_helix_Sheet_helix_Sheet_sheet") == True
    assert checkCanMerge("Coil_helix_Coil_sheet_Coil_turn_Helix_helix_Helix_turn_Helix_coil") == False
    assert checkCanMerge("Coil_helix_Coil_sheet_Coil_turn_Helix_helix_Helix_turn_Turn_coil") == False
    assert checkCanMerge("Coil_helix_Coil_sheet_Coil_turn_Helix_helix_Helix_turn_Helix_coil") == False


    assert checkCanMerge("Helix_helix_Sheet_helix_Coil_helix_Turn_helix_Helix_sheet_Sheet_sheet_Coil_sheet_Turn_sheet") == True
    assert checkCanMerge("Sheet_helix_Sheet_sheet_Sheet_coil_Sheet_turn_Helix_helix_Helix_sheet_Helix_coil_Helix_turn") == True
    assert checkCanMerge("Sheet_helix_Sheet_sheet_Sheet_coil_Turn_turn_Helix_helix_Helix_sheet_Helix_coil_Helix_turn") == False
    assert checkCanMerge("Sheet_helix_Sheet_sheet_Sheet_coil_Coil_turn_Helix_helix_Helix_sheet_Coil_coil_Helix_turn") == False


def checkCanMerge(i):
    words = i.split("_")
    total = len(words)/2
    asite_words = words[1::2]
    exitsite_words = words[::2]
    exitsite = set(words[::2])
    asite = set(words[1::2])
    if total == 2:
        if (len(exitsite) == 1 and len(asite) != 2) and (len(exitsite) != 2 and len(asite) == 1):
                return False
        if len(exitsite) == 2 and len(asite) == 2:
            return False
    elif total == 3:
        if len(exitsite) == 2 or len(asite) == 2:
            return False
        if (len(exitsite) == 3 and len(asite) != 1) or (len(exitsite) != 1 and len(asite) == 3):
            return False

    elif total == 4:
        if len(exitsite) == 3 or len(asite) == 3:
            return False
        if len(exitsite) == 2 and len(asite) == 2:
            good = True
            for j in exitsite:
                if exitsite_words.count(j) != 2:
                    good = False
            for j in asite:
                if asite_words.count(j) != 2:
                    good = False
            return good
        if not ((len(exitsite) == 1 and len(asite) == 4) or (len(exitsite) == 4 and len(asite) == 1)):
            return False
        if (len(exitsite) != 2 and len(asite) == 2) and (len(exitsite) == 2 and len(asite) != 2):
            return False
    elif total == 5 or total == 7 or total == 9 or total == 10 or total == 11 or total == 13 or total == 14 or total == 16:
        return False
    elif total == 6:
        if not ((len(exitsite) == 3 and len(asite) == 2) or (len(exitsite) == 2 and len(asite) == 3)):
            return False
    elif total == 8:
        if not ((len(exitsite) == 4 and len(asite) == 2) or (len(exitsite) == 2 and len(asite) == 4)):
            return False
    elif total == 12:
        if not ((len(exitsite) == 4 and len(asite) == 3) or (len(exitsite) == 3 and len(asite) == 4)):
            return False
    return True


# groups = set(["Coil_coil","Coil_helix","Coil_sheet","Helix_coil","Helix_helix","Helix_sheet","Sheet_coil","Sheet_helix","Sheet_sheet","Turn_coil","Turn_helix","Turn_sheet","Turn_turn","Coil_turn","Helix_turn","Sheet_turn"])
# # #testCheckMerge()
# combos = list(powerset(groups))
# #print checkCanMerge("Coil_helix_Coil_sheet_Helix_helix_Helix_turn")
# for i in combos:
#     if len(i) > 1:
#         combined = "_".join(i)
#         check = checkCanMerge(combined)
      
#         if check:
#             genes_used,new_folder = mergeNonSigGenomes(i,base_dir = "Structure_Downstream_Combo/")
#             if new_folder != None:
#                 createExpressionFiles("Structure_Downstream_Combo/"+new_folder,genes_used)

genes_used,new_folder = mergeNonSigGenomes(["Helix","Sheet"],base_dir="Secondary_structures")
createExpressionFiles("Secondary_structures/"+new_folder,genes_used)

genes_used,new_folder = mergeNonSigGenomes(["Turn","Coil"],base_dir="Secondary_structures")
createExpressionFiles("Secondary_structures/"+new_folder,genes_used)

# genes_used,new_folder = mergeNonSigGenomes(["Helix","Start_helix_sim_0.95_dEta"],base_dir="Secondary_structures_pechmann")
# createExpressionFiles("Secondary_structures_pechmann/"+new_folder,genes_used)

# genes_used,new_folder = mergeNonSigGenomes(["Helix","Start_helix_sim_inverse_dEta"],base_dir="Secondary_structures_pechmann")
# createExpressionFiles("Secondary_structures_pechmann/"+new_folder,genes_used)

# genes_used,new_folder = mergeNonSigGenomes(["Helix","Start_helix_sim_0.5_dEta"],base_dir="Secondary_structures_pechmann")
# createExpressionFiles("Secondary_structures_pechmann/"+new_folder,genes_used)


# genes_used,new_folder = mergeNonSigGenomes(["Helix","Start_helix_sim_0.9_dEta"],base_dir="Secondary_structures_pechmann")
# createExpressionFiles("Secondary_structures_pechmann/"+new_folder,genes_used)

#genes_used,new_folder = mergeNonSigGenomes(["Helix","Start_helix_sim_0.85_dEta"],base_dir="Secondary_structures_pechmann")
#createExpressionFiles("Secondary_structures_pechmann/"+new_folder,genes_used)


# genes_used,new_folder = mergeNonSigGenomes(["Helix","End_helix","Start_sheet","Sheet","End_sheet"],base_dir="Secondary_structures_begin_end")
# createExpressionFiles("Secondary_structures_begin_end/"+new_folder,genes_used)

# genes_used,new_folder = mergeNonSigGenomes(["Start_coil","Coil","End_coil","Turn"],base_dir="Secondary_structures_begin_end")
# createExpressionFiles("Secondary_structures_begin_end/"+new_folder,genes_used)

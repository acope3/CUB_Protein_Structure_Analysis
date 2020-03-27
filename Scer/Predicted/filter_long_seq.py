import pandas as pd


def filterShortSeq():
    with open("Ordered_disordered/Disordered/disordered.fasta") as fin,open("Ordered_disordered/Disordered_25_min/disordered_25_codons_min.fasta",'w') as out:
        seq_id = ""
        genes_to_keep = []
        index = 0
        for line in fin:
            if line[0] == ">":
                seq_id = line
            else:
                if len(line.strip())/3 >= 25:
                    out.write(seq_id)
                    out.write(line)
                    genes_to_keep.append(index)
                index+=1
    phi = pd.read_csv("Ordered_disordered/Disordered/disordered_phi.csv",header=0)
    phi_sub = phi.iloc[genes_to_keep,:]
    phi_sub.to_csv("Ordered_disordered/Disordered_25_min/disordered_25_codons_min_phi.csv",header=True,index=False,quoting=False)
    print phi_sub.iloc[:,1].mean()

filterShortSeq()
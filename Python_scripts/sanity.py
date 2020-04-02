

def parse(filename):
    seq_len = {}
    seq = {}
    with open(filename) as fin:
        for line in fin:
            if line[0] == ">":
                line_spt = line.split()
                gene = line_spt[0][1:]
            else:
                seq_len[gene] = len(line.strip())
                seq[gene] = line.strip()
                assert len(line.strip()) % 3 == 0
    return seq_len,seq



print "Checking..."
# file_1 = "Final_runs/Alpha/Structure_Downstream_Combo_fixed_cutoff_issue/Coil_helix_Helix_helix_Sheet_helix/coil_helix_helix_helix_sheet_helix.fasta"
# file_2 = "Final_runs/Alpha/Secondary_structures_no_mito_sep_nterm/Helix/helix.fasta"

file_1="Scer/Predicted/Ordered_disordered/complete_seq.fasta"
file_2 = "Scer/Predicted/Ordered_disordered/Ordered_Disordered_Nterminus/ordered_disordered_nterminus.fasta"

seq_len_1,seq_1 = parse(file_1)
seq_len_2,seq_2 = parse(file_2)

assert len(seq_1) == len(seq_2)
for i in seq_1.keys():
    len_1 = seq_len_1.get(i)
    len_2 = seq_len_2.get(i)
    seq_a = seq_1.get(i)
    seq_b = seq_2.get(i)
    assert len_1 == len_2
    print i,seq_a.count("A"),seq_b.count("A")
    assert seq_a.count("A") == seq_b.count("A")
    assert seq_a.count("B") == seq_b.count("B")
    assert seq_a.count("G") == seq_b.count("G")
    assert seq_a.count("T") == seq_b.count("T")
print "Done."

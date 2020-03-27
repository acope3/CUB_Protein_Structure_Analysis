import os
import shutil

def appendStartCodonForSim(genome,output):
	with open(genome) as fin,open(output,'w') as out:
		for line in fin:
			if line[0] == ">":
				out.write(line)
			else:
				seq = line.strip()
				if line[0:3] != "ATG":
					seq = "ATG"+seq
				if (line[-3:] != "TAA") and (line[-3:] != "TAG") and (line[-3:] != "TGA"):
					seq = seq + "TAA\n"
				else:
					seq = seq+"\n"
				out.write(seq)

	return None

def stripStopCodon(genome):
	with open(genome) as fin,open("tmp.fasta",'w') as out:
		current_gene = ''
		seq = ''
		for line in fin:
			if line[0] == ">":
				if current_gene != '':
					out.write(current_gene)
					out.write(seq[3:-3])
					out.write("\n")
				current_gene = line
				seq = ''
			else:
				seq = seq + line.strip()
		out.write(current_gene)
		out.write(seq[3:-3])

	os.remove(genome)
	os.rename("tmp.fasta",genome)
	return None


appendStartCodonForSim("Scer/Exp_conservative_homology/Complete_seq/Complete_seq/complete_seq.fasta","complete_seq_append_start_stop.fasta")

#stripStopCodon("Final_runs/Exp_conservative_homology/Secondary_structures_pechmann/Start_helix_sim_0.5_dEta/start_helix_sim_0.5_deta.fasta")
#stripStopCodon("Final_runs/Exp_conservative_homology/Secondary_structures_pechmann/Start_helix_sim_inverse_dEta/start_helix_sim_inverse_deta.fasta")
#stripStopCodon("Final_runs/Exp_conservative_homology/Secondary_structures_pechmann/Start_helix_sim_0.95_dEta/start_helix_sim_0.95_deta.fasta")
#stripStopCodon("Final_runs/Exp_conservative_homology/Secondary_structures_pechmann/Start_helix_sim_0.9_dEta/start_helix_sim_0.9_deta.fasta")
#stripStopCodon("Final_runs/Exp_conservative_homology/Secondary_structures_pechmann/Start_helix_sim_0.85_dEta/start_helix_sim_0.85_deta.fasta")
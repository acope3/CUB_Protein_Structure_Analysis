def main():
	with open("simulated_scer_exp_conservative_homology_seq.fasta",'r') as fin, open("mod_simulated_scer_exp_conservative_homology_seq.fasta",'w') as out:
		for line in fin:
			if line[0] == ">":
				out.write("\n")
				out.write(line)
			else:
				out.write(line.strip())	

				
	return 0

main()




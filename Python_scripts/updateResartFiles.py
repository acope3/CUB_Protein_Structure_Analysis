import sys
with open(sys.argv[1]) as fin,open(sys.argv[2],'w') as out:
    categories = None
    mutation_prior_sd = 0
    get_cat = False
    get_mpsd = False
    for line in fin:
        if line.strip() == ">numMutationCategories:":
            get_cat = True
            out.write(">numMutationCategories:\n")
        elif line.strip() == ">mutation_prior_sd:":
            get_mpsd = True
        elif get_cat:
            categories = int(line.strip())
            out.write(line)
            get_cat = False
        elif get_mpsd:
            mutation_prior_sd = line.strip()
            out.write(">mutation_prior_mean:\n***\n")
            for i in xrange(0,categories):
                for j in xrange(0,40):
                    out.write("0.0")
                    if j % 10 != 9:
                        out.write(" ")
                    else:
                        out.write("\n")
            out.write(">mutation_prior_sd:\n")
            for i in xrange(0,categories):
                out.write("***\n")
                for j in xrange(0,40):
                    out.write(mutation_prior_sd)
                    if j % 10 != 9:
                        out.write(" ")
                    else:
                        out.write("\n")
            get_mpsd = False
        else:
            out.write(line)

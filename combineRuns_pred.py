import re
import glob

alpha = re.compile(r"(Coil|Helix|Sheet)_helix_(Coil|Helix|Sheet)_helix_(Coil|Helix|Sheet)_helix")
beta = re.compile(r"(Coil|Helix|Sheet)_sheet_(Coil|Helix|Sheet)_sheet_(Coil|Helix|Sheet)_sheet")
coil = re.compile(r"(Coil|Helix|Sheet)_coil_(Coil|Helix|Sheet)_coil_(Coil|Helix|Sheet)_coil")
turn = re.compile(r"(Coil|Helix|Sheet)_turn_(Coil|Helix|Sheet)_turn_(Coil|Helix|Sheet)_turn")

ds_alpha = re.compile(r"Helix_(coil|helix|sheet)_Helix_(coil|helix|sheet)_Helix_(coil|helix|sheet)")
ds_beta = re.compile(r"Sheet_(coil|helix|sheet)_Sheet_(coil|helix|sheet)_Sheet_(coil|helix|sheet)")
ds_coil = re.compile(r"Coil_(coil|helix|sheet)_Coil_(coil|helix|sheet)_Coil_(coil|helix|sheet)")


asite_ch = re.compile(r"(Coil|Helix|Sheet)_(helix|coil)_(Coil|Helix|Sheet)_(helix|coil)_(Coil|Helix|Sheet)_(helix|coil)_(Coil|Helix|Sheet)_(helix|coil)_(Coil|Helix|Sheet)_(helix|coil)_(Coil|Helix|Sheet)_(helix|coil)")
asite_ce = re.compile(r"(Coil|Helix|Sheet)_(sheet|coil)_(Coil|Helix|Sheet)_(sheet|coil)_(Coil|Helix|Sheet)_(sheet|coil)_(Coil|Helix|Sheet)_(sheet|coil)_(Coil|Helix|Sheet)_(sheet|coil)_(Coil|Helix|Sheet)_(sheet|coil)")
asite_he = re.compile(r"(Coil|Helix|Sheet)_(helix|sheet)_(Coil|Helix|Sheet)_(helix|sheet)_(Coil|Helix|Sheet)_(helix|sheet)_(Coil|Helix|Sheet)_(helix|sheet)_(Coil|Helix|Sheet)_(helix|sheet)_(Coil|Helix|Sheet)_(helix|sheet)")

exit_ch = re.compile(r"(Coil|Helix)_(sheet|helix|coil)_(Coil|Helix)_(sheet|helix|coil)_(Coil|Helix)_(sheet|helix|coil)_(Coil|Helix)_(sheet|helix|coil)_(Coil|Helix)_(sheet|helix|coil)_(Coil|Helix)_(sheet|helix|coil)")
exit_ce = re.compile(r"(Coil|Sheet)_(sheet|coil|helix)_(Coil|Sheet)_(sheet|coil|helix)_(Coil|Sheet)_(sheet|coil|helix)_(Coil|Sheet)_(sheet|coil|helix)_(Coil|Sheet)_(sheet|coil|helix)_(Coil|Sheet)_(sheet|coil|helix)")
exit_eh = re.compile(r"(Helix|Sheet)_(helix|sheet|coil)_(Helix|Sheet)_(helix|sheet|coil)_(Helix|Sheet)_(helix|sheet|coil)_(Helix|Sheet)_(helix|sheet|coil)_(Helix|Sheet)_(helix|sheet|coil)_(Helix|Sheet)_(helix|sheet|coil)")



files = glob.glob("Final_runs/Exp_liberal_homology_pred_structures/Structure_Downstream_Combo/*")
print files
runs = {}
for f in files:
    run_name = f.split("/")[-1]
    result_a = alpha.match(run_name)
    result_b = beta.match(run_name)
    result_c = coil.match(run_name)

    result_dsa = ds_alpha.match(run_name)
    result_dsb = ds_beta.match(run_name)
    result_dsc = ds_coil.match(run_name)
    

    result_asite_ch = asite_ch.match(run_name)
    result_asite_ce = asite_ce.match(run_name)
    result_asite_he = asite_he.match(run_name)
    

    result_exit_ch = exit_ch.match(run_name)
    result_exit_ce = exit_ce.match(run_name)
    result_exit_eh = exit_eh.match(run_name)
    


    if result_a != None:
        value = result_a.group()
        if value[-1] != "_":
            runs["Helix"] = result_a.group()
    if result_b != None:
        value = result_b.group()
        if value[-1] != "_":
            runs["Sheet"] = result_b.group()
    if result_c != None:
        value = result_c.group()
        if value[-1] != "_":
            runs["Coil"] = result_c.group()
    

    if result_dsa != None:
        value = result_dsa.group()
        if value[-1] != "_":
            runs["Downstream_Helix"] = result_dsa.group()
    if result_dsb != None:
        value = result_dsb.group()
        if value[-1] != "_":
            runs["Downstream_Sheet"] = result_dsb.group()
    if result_dsc != None:
        value = result_dsc.group()
        if value[-1] != "_":
            runs["Downstream_Coil"] = result_dsc.group()
   

    if result_asite_ch != None:
        value = result_asite_ch.group()
        if value[-1] != "_":
            runs["Coil_Helix"] = result_asite_ch.group()
    if result_asite_ce != None:
        value = result_asite_ce.group()
        if value[-1] != "_":
            runs["Coil_Sheet"] = result_asite_ce.group()
    if result_asite_he != None:
        value = result_asite_he.group()
        if value[-1] != "_":
            runs["Helix_Sheet"] = result_asite_he.group()

    if result_exit_ch != None:
        value = result_exit_ch.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Helix"] = result_exit_ch.group()
    if result_exit_ce != None:
        value = result_exit_ce.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Sheet"] = result_exit_ce.group()
    if result_exit_eh != None:
        value = result_exit_eh.group()
        if value[-1] != "_":
            runs["Downstream_Sheet_Helix"] = result_exit_eh.group()




mergings = [["Helix","Sheet","Coil"],
            ["Helix","Coil_Sheet"],
            ["Sheet","Coil_Helix"],
            ["Coil","Helix_Sheet"],
            ["Downstream_Helix","Downstream_Sheet","Downstream_Coil"],
            ["Downstream_Helix","Downstream_Coil_Sheet"],
            ["Downstream_Sheet","Downstream_Coil_Helix"],
            ["Downstream_Coil","Downstream_Sheet_Helix"]]


print runs
fits_to_run = set(["Nterminus"])
with open("runs_to_compare_exp_proteins_predicted_structures_scer.txt","w") as out:
    for merge in mergings:
        out.write(",".join([runs[x] for x in merge]))
        for x in merge:
            fits_to_run.add(runs[x])
        out.write(",Nterminus\n")

with open("runs_for_exp_proteins_predicted_structures_scer.txt",'w') as out:
    for i in fits_to_run:
        out.write(i)
        out.write("\n")

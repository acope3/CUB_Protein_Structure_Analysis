import re
import glob

coil = re.compile(r"((Sheet|Coil|Helix|Turn)_coil_*){4}")
alpha = re.compile(r"((Sheet|Coil|Helix|Turn)_helix_*){4}")
beta = re.compile(r"((Sheet|Coil|Helix|Turn)_sheet_*){4}")
turn = re.compile(r"((Sheet|Coil|Helix|Turn)_turn_*){4}")

ds_coil = re.compile(r"(Coil_(sheet|coil|helix|turn)_*){4}")
ds_alpha = re.compile(r"(Helix_(sheet|coil|helix|turn)_*){4}")
ds_beta = re.compile(r"(Sheet_(sheet|coil|helix|turn)_*){4}")
ds_turn = re.compile(r"(Turn_(sheet|coil|helix|turn)_*){4}")

asite_ch = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|helix)_*){8}")
asite_ce = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|sheet)_*){8}")
asite_ct = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|turn)_*){8}")
asite_he = re.compile(r"((Sheet|Coil|Helix|Turn)_(helix|sheet)_*){8}")
asite_ht = re.compile(r"((Sheet|Coil|Helix|Turn)_(helix|turn)_*){8}")
asite_et = re.compile(r"((Sheet|Coil|Helix|Turn)_(sheet|turn)_*){8}")

exit_ch = re.compile(r"((Coil|Helix)_(sheet|coil|helix|turn)_*){8}")
exit_ce = re.compile(r"((Coil|Sheet)_(sheet|coil|helix|turn)_*){8}")
exit_ct = re.compile(r"((Coil|Turn)_(sheet|coil|helix|turn)_*){8}")
exit_eh = re.compile(r"((Sheet|Helix)_(sheet|coil|helix|turn)_*){8}")
exit_th = re.compile(r"((Turn|Helix)_(sheet|coil|helix|turn)_*){8}")
exit_et = re.compile(r"((Sheet|Turn)_(sheet|coil|helix|turn)_*){8}")

asite_che = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|helix|sheet)_*){12}")
asite_cht = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|helix|turn)_*){12}")
asite_cte = re.compile(r"((Sheet|Coil|Helix|Turn)_(coil|sheet|turn)_*){12}")
asite_hte = re.compile(r"((Sheet|Coil|Helix|Turn)_(helix|sheet|turn)_*){12}")


exit_che = re.compile(r"((Coil|Helix|Sheet)_(sheet|coil|helix|turn)_*){12}")
exit_cet = re.compile(r"((Coil|Sheet|Turn)_(sheet|coil|helix|turn)_*){12}")
exit_cth = re.compile(r"((Coil|Turn|Helix)_(sheet|coil|helix|turn)_*){12}")
exit_the = re.compile(r"((Turn|Helix|Sheet)_(sheet|coil|helix|turn)_*){12}")

files = glob.glob("Ecoli/Exp_conservative_homology/Structure_Downstream_Combo/*")
runs = {}
for f in files:
    run_name = f.split("/")[-1]
    result_a = alpha.match(run_name)
    result_b = beta.match(run_name)
    result_c = coil.match(run_name)
    result_t = turn.match(run_name)

    result_dsa = ds_alpha.match(run_name)
    result_dsb = ds_beta.match(run_name)
    result_dsc = ds_coil.match(run_name)
    result_dst = ds_turn.match(run_name)

    result_asite_ch = asite_ch.match(run_name)
    result_asite_ce = asite_ce.match(run_name)
    result_asite_he = asite_he.match(run_name)
    result_asite_ct = asite_ct.match(run_name)
    result_asite_ht = asite_ht.match(run_name)
    result_asite_et = asite_et.match(run_name)

    result_exit_ch = exit_ch.match(run_name)
    result_exit_ce = exit_ce.match(run_name)
    result_exit_eh = exit_eh.match(run_name)
    result_exit_ct = exit_ct.match(run_name)
    result_exit_th = exit_th.match(run_name)
    result_exit_et = exit_et.match(run_name)

    result_asite_che = asite_che.match(run_name)
    result_asite_cht = asite_cht.match(run_name)
    result_asite_cte = asite_cte.match(run_name)
    result_asite_hte = asite_hte.match(run_name)


    result_exit_che = exit_che.match(run_name)
    result_exit_cet = exit_cet.match(run_name)
    result_exit_cth = exit_cth.match(run_name)
    result_exit_the = exit_the.match(run_name)


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
    if result_t != None:
        value = result_t.group()
        if value[-1] != "_":
            runs["Turn"] = result_t.group()

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
    if result_dst != None:
        value = result_dst.group()
        if value[-1] != "_":
            runs["Downstream_Turn"] = result_dst.group()

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
    if result_asite_ct != None:
        value = result_asite_ct.group()
        if value[-1] != "_":
            runs["Coil_Turn"] = result_asite_ct.group()
    if result_asite_et != None:
        value = result_asite_et.group()
        if value[-1] != "_":
            runs["Turn_Sheet"] = result_asite_et.group()
    if result_asite_ht != None:
        value = result_asite_ht.group()
        if value[-1] != "_":
            runs["Helix_Turn"] = result_asite_ht.group()

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
    if result_exit_ct != None:
        value = result_exit_ct.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Turn"] = result_exit_ct.group()
    if result_exit_et != None:
        value = result_exit_et.group()
        if value[-1] != "_":
            runs["Downstream_Turn_Sheet"] = result_exit_et.group()
    if result_exit_th != None:
        value = result_exit_th.group()
        if value[-1] != "_":
            runs["Downstream_Helix_Turn"] = result_exit_th.group()

    if result_asite_che != None:
        value = result_asite_che.group()
        if value[-1] != "_":
            runs["Coil_Helix_Sheet"] = result_asite_che.group()
    if result_asite_cht != None:
        value = result_asite_cht.group()
        if value[-1] != "_":
            runs["Coil_Turn_Helix"] = result_asite_cht.group()
    if result_asite_cte != None:
        value = result_asite_cte.group()
        if value[-1] != "_":
            runs["Coil_Turn_Sheet"] = result_asite_cte.group()
    if result_asite_hte != None:
        value = result_asite_hte.group()
        if value[-1] != "_":
            runs["Helix_Turn_Sheet"] = result_asite_hte.group()


    if result_exit_che != None:
        value = result_exit_che.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Helix_Sheet"] = result_exit_che.group()
    if result_exit_cet != None:
        value = result_exit_cet.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Turn_Sheet"] = result_exit_cet.group()
    if result_exit_cth != None:
        value = result_exit_cth.group()
        if value[-1] != "_":
            runs["Downstream_Coil_Turn_Helix"] = result_exit_cth.group()
    if result_exit_the != None:
        value = result_exit_the.group()
        if value[-1] != "_":
            runs["Downstream_Helix_Turn_Sheet"] = result_exit_the.group()

mergings = [["Helix","Sheet","Coil","Turn"],
            ["Coil_Helix_Sheet","Turn"],
            ["Coil_Turn_Sheet","Helix"],
            ["Coil_Turn_Helix","Sheet"],
            ["Helix_Turn_Sheet","Coil"],
            ["Helix","Coil_Sheet","Turn"],
            ["Sheet","Coil_Helix","Turn"],
            ["Coil","Helix_Sheet","Turn"],
            ["Helix","Coil_Turn","Sheet"],
            ["Sheet","Helix_Turn","Coil"],
            ["Coil","Turn_Sheet","Helix"],
            ["Helix_Sheet","Coil_Turn"],
            ["Coil_Sheet","Helix_Turn"],
            ["Turn_Sheet","Coil_Helix"],
            ["Downstream_Helix","Downstream_Sheet","Downstream_Coil","Downstream_Turn"],
            ["Downstream_Coil_Helix_Sheet","Downstream_Turn"],
            ["Downstream_Coil_Turn_Sheet","Downstream_Helix"],
            ["Downstream_Coil_Turn_Helix","Downstream_Sheet"],
            ["Downstream_Helix_Turn_Sheet","Downstream_Coil"],
            ["Downstream_Helix","Downstream_Coil_Sheet","Downstream_Turn"],
            ["Downstream_Sheet","Downstream_Coil_Helix","Downstream_Turn"],
            ["Downstream_Coil","Downstream_Sheet_Helix","Downstream_Turn"],
            ["Downstream_Helix","Downstream_Coil_Turn","Downstream_Sheet"],
            ["Downstream_Sheet","Downstream_Helix_Turn","Downstream_Coil"],
            ["Downstream_Coil","Downstream_Turn_Sheet","Downstream_Helix"],
            ["Downstream_Sheet_Helix","Downstream_Coil_Turn"],
            ["Downstream_Coil_Sheet","Downstream_Helix_Turn"],
            ["Downstream_Turn_Sheet","Downstream_Coil_Helix"]]

fits_to_run = set(["Nterminus"])
with open("runs_to_compare_exp_conservative_ecoli.txt","w") as out:
    for merge in mergings:
        out.write(",".join([runs[x] for x in merge]))
        for x in merge:
            fits_to_run.add(runs[x])
        out.write(",Nterminus\n")

with open("runs_for_exp_conservative_ecoli.txt",'w') as out:
    for i in fits_to_run:
        out.write(i)
        out.write("\n")

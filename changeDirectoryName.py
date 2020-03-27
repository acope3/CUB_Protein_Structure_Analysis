import re
import glob
import os
import natsort

loc = glob.glob("Final_runs/Alpha/Results/Structure_Downstream_Combo/*/")
for dir in loc:
    runs = glob.glob(dir+"run_*/")
    ordered_runs = natsort.natsorted(runs)
    final_run = ordered_runs[-1]
    os.rename(final_run,dir+"final_run")

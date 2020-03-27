import re
import glob


files = glob.glob("Final_runs/Alpha/Structure_Downstream_Combo/*/")
struct = re.compile("[A-Z][a-z]+_(helix)|(sheet)_[A-Z][a-z]+_(sheet)|(helix)_[A-Z][a-z]+_(sheet)|(helix)_[A-Z][a-z]+_(sheet)|(helix)_[A-Z][a-z]+_(sheet)|(helix)_[A-Z][a-z]+_(sheet)|(helix)")
helix = re.compile("[A-Z][a-z]+_(helix)_[A-Z][a-z]+_(helix)_[A-Z][a-z]+_(helix)")
sheet = re.compile("[A-Z][a-z]+_(sheet)_[A-Z][a-z]+_(sheet)_[A-Z][a-z]+_(sheet)")
coil = re.compile("[A-Z][a-z]+_(coil)_[A-Z][a-z]+_(coil)_[A-Z][a-z]+_(coil)")

for f in files:
    category = f.split("/")[-1]
    check = category.split("_")
    if len(check) <= 6:

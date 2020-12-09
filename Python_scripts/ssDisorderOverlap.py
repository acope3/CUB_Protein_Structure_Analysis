import glob
import pandas as pd

files = glob.glob("../Data/Psipred_iupred/*.ss2")
c_o = 0
c_d = 0
h_o = 0
h_d = 0
e_o = 0
e_d = 0
for f in files:
	df = pd.read_csv(f,header=None,delim_whitespace=True)
	c = df.loc[df.iloc[:,2] == "C",:]
	h = df.loc[df.iloc[:,2] == "H",:]
	e = df.loc[df.iloc[:,2] == "E",:]

	c_o += len(c.loc[c.iloc[:,6] == "O",:].index)
	c_d += len(c.loc[c.iloc[:,6] == "D",:].index)

	h_o += len(h.loc[h.iloc[:,6] == "O",:].index)
	h_d += len(h.loc[h.iloc[:,6] == "D",:].index)

	e_o += len(e.loc[e.iloc[:,6] == "O",:].index)
	e_d += len(e.loc[e.iloc[:,6] == "D",:].index)


print c_o,c_d
print h_o,h_d
print e_o,e_d

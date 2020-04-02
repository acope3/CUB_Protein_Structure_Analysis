import glob
import pandas as pd

ss_files = glob.glob("Alpha_psipred/*.ss2")

freq = {"H":{"H":0,"E":0,"C":0},"E":{"H":0,"E":0,"C":0},"C":{"H":0,"E":0,"C":0}}

for ss in ss_files:
	df = df = pd.read_csv(ss,header=-1,delim_whitespace=True,skiprows=2)
	length = len(df.index)
	for index,row in df.iterrows():
		print index
		if (index + 35) >= length:
			break
		else:
			upstream = df.iloc[index,2]
			a_site = df.iloc[index+35,2]
			freq[upstream][a_site] += 1
print freq
##Result: {'H': {'H': 574518, 'C': 435852, 'E': 72119}, 'C': {'H': 449931, 'C': 770870, 'E': 142914}, 'E': {'H': 75036, 'C': 139112, 'E': 61246}}
##					 53%           40%          7%                  33%         56%           1%                27%            50%         22%
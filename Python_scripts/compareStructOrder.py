import re
import pandas as pd 
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind,mannwhitneyu,spearmanr
import numpy as np 


def compareHelixCoilAnchorOrdDisDist():

	files = glob.glob("Psipred_iupred_w_Anchor/*.ss2")
	helix_dis_score = []
	coil_dis_score = []
	for f in files:
		prot_struct = pd.read_csv(f,sep="\t",header=0)
		helix_dis = prot_struct.loc[(prot_struct.iloc[:,2] == "H") & (prot_struct.loc[:,"Order"] == "D") & (prot_struct.loc[:,"H"] > 0.75),:]
		coil_dis = prot_struct.loc[(prot_struct.iloc[:,2] == "C") & (prot_struct.loc[:,"Order"] == "D") & (prot_struct.loc[:,"C"] > 0.75),:]
		helix_dis_score = helix_dis_score + list(helix_dis.loc[:,"Anchor_Score"].values)
		coil_dis_score = coil_dis_score + list(coil_dis.loc[:,"Anchor_Score"].values)
	test_stat,p = mannwhitneyu(helix_dis_score,coil_dis_score,alternative="two-sided")
	print p
	print np.median(helix_dis_score),np.median(coil_dis_score)
	plt.hist(helix_dis_score,alpha=0.5,label="Disordered Helix",density=True)
	plt.hist(coil_dis_score,alpha=0.5,label="Disordered Coil",density=True)
	plt.legend(loc="upper left")
	plt.xlabel("Anchor Score")
	plt.ylabel("Density")
	plt.title("Comparison of Disordered Helix and Coil Anchor Scores (SS Score > 0.75)")
	plt.savefig("helix_coil_scores_anchor.pdf")

def compareHelixOrdDisDist():

	files = glob.glob("Psipred_iupred_w_Anchor/*.ss2")
	helix_ord_score = []
	helix_dis_score = []
	for f in files:
		prot_struct = pd.read_csv(f,sep="\t",header=0)
		helix_ord = prot_struct.loc[(prot_struct.iloc[:,2] == "H") & (prot_struct.loc[:,"Order"] == "O") ,:]
		helix_dis = prot_struct.loc[(prot_struct.iloc[:,2] == "H") & (prot_struct.loc[:,"Order"] == "D") ,:]
		helix_ord_score = helix_ord_score + list(helix_ord.loc[:,"Anchor_Score"].values)
		helix_dis_score = helix_dis_score + list(helix_dis.loc[:,"Anchor_Score"].values)
	test_stat,p = mannwhitneyu(helix_ord_score,helix_dis_score,alternative="two-sided")
	print p
	print np.median(helix_ord_score),np.median(helix_dis_score)
	plt.hist(helix_ord_score,alpha=0.5,label="Helix and Ordered",density=True)
	plt.hist(helix_dis_score,alpha=0.5,label="Helix and Disordered",density=True)
	plt.legend(loc="upper left")
	plt.xlabel("Anchor Score")
	plt.ylabel("Density")
	plt.title("Comparison of Helix Anchor Scores: Ordered vs. Disordered")
	plt.savefig("helix_scores_anchor.pdf")


def compareHelixOrdDisLength():
	files = glob.glob("Psipred_iupred/*.ss2")
	helix_length = []
	helix_dis_frac= []
	pat = re.compile("H+")
	for f in files:
		prot_struct = pd.read_csv(f,sep="\t",header=None)
		ss = "".join(prot_struct.iloc[:,2].values)
		order = "".join(prot_struct.iloc[:,6].values)
		helix = pat.finditer(ss)
		for h in helix:
			begin = h.start()
			end = h.end()
			sub = order[begin:end]
			total_dis = float(sub.count("D"))/len(sub)
			helix_length.append(len(sub))
			helix_dis_frac.append(total_dis)
	rho,p = spearmanr(helix_length,helix_dis_frac)
	print(rho,p)
	plt.scatter(np.log10(helix_length),helix_dis_frac,alpha=0.5)
	plt.xlabel("Log10(Helix Length)")
	plt.ylabel("Fraction Disordered")
	plt.title("Fraction of disordered regions in helices as function of length")
	plt.text(s=r'$\rho$ = '+str(round(rho,3)),x=2.5,y=0.9)
	plt.savefig("helix_length_disordered_prop.pdf")



def mergePsiPredIUPRED():
	mapping = {"C":{"O":0,"D":0},"H":{"O":0,"D":0},"E":{"O":0,"D":0}}
	pat = re.compile("C+")
	order = glob.glob("IUPRED2/IUPRED_w_Anchor/*")
	for f in order:
		prot = f.split("/")[-1]
		ss2 = "Alpha_psipred/"+prot+".ss2"
		df_ss = pd.read_csv(ss2,header=None,delim_whitespace=True,skiprows=2)
		df_ss.columns = ["Position","AA","SS","C","H","E"]
		df_order = pd.read_csv(f,header=None,sep="\t")
		assert len(df_ss.index) == len(df_order.index)
		scores = df_order.iloc[:,2].values
		anchor = df_order.iloc[:,3].values
		df_order["Order_Score"] = scores
		df_order.loc[df_order.loc[:,"Order_Score"] < 0.5,2] = "O"
		df_order.loc[df_order.loc[:,"Order_Score"] >= 0.5,2] = "D"
		ss = df_ss.iloc[:,2].values
		ss_string = "".join(ss)
		order = df_order.iloc[:,2].values
		df_ss = pd.concat([df_ss,pd.DataFrame({"Order":order,"Order_Score":scores,"Anchor_Score":anchor})],axis=1)
		order_string = "".join(order)
		df_ss.to_csv("Psipred_iupred_w_Anchor/"+prot+".ss2",sep="\t",header=True,index=False,quoting=False)
	# for i,j in zip(ss,order):
	# 	mapping[i][j] += 1

	# coils = pat.finditer(ss_string)

	# for i in coils:
	# 	subset = order_string[i.start():i.end()]
	# 	ord_count = subset.count("O")
	# 	dis_count = subset.count("D")
	# 	assert (ord_count + dis_count) == len(subset)
	# 	coil_info = coil_info.append({"Protein":prot,"Length":len(subset),"Order":ord_count,"Disorder":dis_count},ignore_index=True)
##coil_info.to_csv("coil_order_disorder.tsv",sep="\t",header=True,index=False,quoting=False)
##print mapping
#mergePsiPredIUPRED()
#compareHelixOrdDisDist()
compareHelixCoilAnchorOrdDisDist()
#compareHelixOrdDisLength()
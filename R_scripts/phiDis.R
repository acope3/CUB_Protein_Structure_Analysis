library(ggplot2)
library(ggpubr)
# phi <- read.table("../mod_scerevisiae_expression_wo_phi_allunique.csv",sep=",",header=T,stringsAsFactors=F)
# dis <- read.table("../scer_disordered_freq.tsv",sep="\t",header=T,stringsAsFactors=F)

phi <- read.table("../Ecoli_K12_MG1655_main_phi.csv",sep=",",header=T,stringsAsFactors=F)
dis <- read.table("../ecoli_disordered_freq.tsv",sep="\t",header=T,stringsAsFactors=F)


phi[,1] <- regmatches(phi[,1], regexpr("[N|Y]P_[0-9]+.[0-9]", phi[,1])) 

phi.dis <- merge(phi,dis,by.x="Gene_id",by.y="Protein")

p <- ggplot(phi.dis,aes(x=Disordered,y=Phi)) + geom_point(alpha=0.7) +
	 xlab("Disordered Frequency") +
	 ylab(expression(phi)) +
	 ggtitle("Comparing Relationship Between Phi\nand Disordered Frequency\nE. coli") +
	 stat_cor(method="spearman",label.sep="\n",label.x.npc="center",label.y.npc="top") + 
	 cowplot::theme_cowplot()

ggsave("../ecoli_compare_phi_dis.pdf")
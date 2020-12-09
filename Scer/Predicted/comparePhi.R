
ordered.phi <- read.table("Ordered_disordered/Ordered_25_min/ordered_25_codons_min_phi.csv",sep=",",header=T,stringsAsFactors=F)
disordered.phi <- read.table("Ordered_disordered/Disordered_25_min/disordered_25_codons_min_phi.csv",sep=",",header=T,stringsAsFactors=F)

ord.est <- read.table("Results/Ordered_est_phi_25_min/run_1/Parameter_est/gene_expression.txt",sep=",",header=T,stringsAsFactors=F)
dis.est <- read.table("Results/Disordered_est_phi_25_min/run_1/Parameter_est/gene_expression.txt",sep=",",header=T,stringsAsFactors=F)

ordered.phi[,"Ord.phi"] <- ord.est[,1]
disordered.phi[,"Dis.phi"] <- dis.est[,1]
ordered.phi[,"Log.Ord.phi"] <- log(ord.est[,1])
disordered.phi[,"Log.Dis.phi"] <- log(dis.est[,1])

phi <- merge(ordered.phi,disordered.phi,by=c("Gene_id","Phi"),all=T,sort=F)
phi[,"GeomMean"] <- exp(rowMeans(phi[,c("Log.Ord.phi","Log.Dis.phi")],na.rm=T))

pdf("compare_phi_comp_seq_ord_dis_split_25_min.pdf")
plot(phi$Phi,phi$GeomMean)
dev.off()


yas <- read.table("~/CUB_structure_project/yassour_w_cds_id.csv",sep=",",header=T,stringsAsFactors=F)
yas <- merge(phi,yas,by.x="Gene_id",by.y="ORF")

pdf("compare_to_yassour.pdf")
for (i in 8:ncol(yas))
{
	x <- which(is.finite(log10(yas[,i])))
	plot(log10(yas[x,i]),log10(yas[x,"Phi"]),main=paste("Yassour et. al:\n",colnames(yas)[i]),xlab=expression("Log"[10]*"[mRNA]"),ylab=expression("Log"[10]*phi~"(Full Sequence)"))
	xmax <- max(log10(yas[x,i]),na.rm=T)
	xmin <- min(log10(yas[x,i]),na.rm=T)
	ymax <- max(log10(yas[x,"Phi"]),na.rm=T)
	ymin <- min(log10(yas[x,"Phi"]),na.rm=T)
	height <-  ymax- ymin 
	width <-  xmax- xmin
	corr <- round(cor(log10(yas[x,i]),log10(yas[x,"Phi"]),use="pairwise.complete.obs"),4)
	text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))

	plot(log10(yas[,i]),log10(yas$GeomMean),main=paste("Yassour et. al:\n",colnames(yas)[i]),xlab=expression("Log"[10]*"[mRNA]"),ylab=expression("Log"[10]*phi~"(Geometric Mean, Ordered and Disordered)"))
	ymax <- max(log10(yas[x,"GeomMean"]),na.rm=T)
	ymin <- min(log10(yas[x,"GeomMean"]),na.rm=T)
	height <-  ymax- ymin 
	width <-  xmax- xmin
	corr <- round(cor(log10(yas[x,i]),log10(yas[x,"GeomMean"]),use="pairwise.complete.obs"),4)
	text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))
}
dev.off()



pdf("compare_to_weinberg.pdf")

weinberg <- read.table("../../weinberg_rpf.csv",sep=",",header=T,stringsAsFactors=F)
weinberg <- merge(phi,weinberg,by.x="Gene_id",by.y="CDS")

x <- which(weinberg[,"RPF"]!=0)
plot(log10(weinberg[x,"RPF"]),log10(weinberg[x,"Phi"]),main="Weinberg et. al:\nProtein Production Rates",xlab=expression("Log"[10]*"Protein Production Rates"),ylab=expression("Log"[10]*phi~"(Full Sequence)"))
xmax <- max(log10(weinberg[x,"RPF"]),na.rm=T)
xmin <- min(log10(weinberg[x,"RPF"]),na.rm=T)
ymax <- max(log10(weinberg[x,"Phi"]),na.rm=T)
ymin <- min(log10(weinberg[x,"Phi"]),na.rm=T)
height <-  ymax- ymin 
width <-  xmax- xmin
corr <- round(cor(log10(weinberg[x,"RPF"]),log10(weinberg[x,"Phi"]),use="pairwise.complete.obs"),4)
text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))


plot(log10(weinberg[x,"RPF"]),log10(weinberg[x,"GeomMean"]),main="Weinberg et. al:\nProtein Production Rates",xlab=expression("Log"[10]*"Protein Production Rates"),ylab=expression("Log"[10]*phi~"(Geometric Mean, Ordered and Disordered)"))
xmax <- max(log10(weinberg[x,"RPF"]),na.rm=T)
xmin <- min(log10(weinberg[x,"RPF"]),na.rm=T)
ymax <- max(log10(weinberg[x,"GeomMean"]),na.rm=T)
ymin <- min(log10(weinberg[x,"GeomMean"]),na.rm=T)
height <-  ymax- ymin 
width <-  xmax- xmin
corr <- round(cor(log10(weinberg[x,"RPF"]),log10(weinberg[x,"GeomMean"]),use="pairwise.complete.obs"),4)
text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))

weinberg <- read.table("../../weinberg_rna.csv",sep=",",header=T,stringsAsFactors=F)
weinberg <- merge(phi,weinberg,by.x="Gene_id",by.y="CDS")
x <- which(weinberg[,"RPF"] !=0 )
plot(log10(weinberg[x,"RPF"]),log10(weinberg[x,"Phi"]),main="Weinberg et. al:\nmRNA Abundances",xlab=expression("Log"[10]*"[mRNA]"),ylab=expression("Log"[10]*phi~"(Full Sequence)"))
xmax <- max(log10(weinberg[x,"RPF"]),na.rm=T)
xmin <- min(log10(weinberg[x,"RPF"]),na.rm=T)
ymax <- max(log10(weinberg[x,"Phi"]),na.rm=T)
ymin <- min(log10(weinberg[x,"Phi"]),na.rm=T)
height <-  ymax- ymin 
width <-  xmax- xmin
corr <- round(cor(log10(weinberg[x,"RPF"]),log10(weinberg[x,"Phi"]),use="pairwise.complete.obs"),4)
text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))


plot(log10(weinberg[x,"RPF"]),log10(weinberg[x,"GeomMean"]),main="Weinberg et. al:\nmRNA Abundances",xlab=expression("Log"[10]*"[mRNA]"),ylab=expression("Log"[10]*phi~"(Geometric Mean, Ordered and Disordered)"))
xmax <- max(log10(weinberg[x,"RPF"]),na.rm=T)
xmin <- min(log10(weinberg[x,"RPF"]),na.rm=T)
ymax <- max(log10(weinberg[x,"GeomMean"]),na.rm=T)
ymin <- min(log10(weinberg[x,"GeomMean"]),na.rm=T)
height <-  ymax- ymin 
width <-  xmax- xmin
corr <- round(cor(log10(weinberg[x,"RPF"]),log10(weinberg[x,"GeomMean"]),use="pairwise.complete.obs"),4)
text(x=xmin+0.1*width,y=ymax-height*0.1,label=bquote(rho~"="~.(corr)))

dev.off()
library(ggplot2)
library(cowplot)

pech <- read.table("../pechman_codon_opt_table.tsv",sep="\t",header=T,stringsAsFactors=F)
filter <- read.table("../Data/pechmann_prot_id.txt",sep="",header=F,stringsAsFactors=F)

odds.ratio.pech <- c(0.916,1.281,1.111,0.796)
or.plots <- list(length=4)
p.plots <- list(length=4)
binwidth <- c(0.02,0.05,0.05,0.03)
for (j in 1:4)
{
	p.value <- data.frame(Position=rep(1,100),P=numeric(length=100))
	odds <- data.frame(Position=rep(1,100),Odds=numeric(length=100))
	for (i in 1:100)
	{
		helix <- read.table(file.path("..","Scer","Test_Fisher_Exact",i,"helix.csv"),sep=",",header=T,stringsAsFactors=F)
		helix <- helix[which(helix$Gene %in% filter[,1]),]
		helix.opt <- merge(helix,pech,by.x="Codon",by.y="codon")
		helix.opt[,"Target.pos"] <- ifelse(helix.opt$Pos == j,as.character(j),"other")
		cont.table <- table(helix.opt$Target.pos,helix.opt$Scer)
		fisher <- fisher.test(cont.table)
		p.value[i,"P"] <- fisher$p.value
		odds[i,"Odds"] <- fisher$estimate
	}

	# pdf(paste0("alpha_helix_p_value_odd_ratio_distribution_simulated_data_",j,".pdf"))
	# hist(p.value,xlab="p-value",main=paste0("Distribution of p-values:\nHelix Relative Position ",j),breaks=20)
	# hist(odds,xlab="Odds Ratio",main=paste0("Distribution of odds ratios:\nHelix Relative Position ",j),breaks=20)
	# abline(v=odds.ratio.pech[j],lty=2,col="red")
	#dev.off()
	if (max(odds$Odds) < 1)
	{
		range.hist <- c(min(odds$Odds)-0.05,1)
	} else if (min(odds$Odds) > 1){
		range.hist <- c(1,max(odds$Odds)+0.05)
	} else {
		range.hist <- range(odds$Odds) + c(-0.05,0.05)
	}
	p <- ggplot(odds,aes(x=Odds)) + geom_histogram(binwidth=binwidth[j],color="black", fill="cyan") + ggtitle(paste("Position",j)) + xlab("Odds Ratio") + ylab("Counts") + xlim(range.hist) + theme_cowplot() + theme(legend.position="none")
	p <- p + geom_vline(xintercept=odds.ratio.pech[j],linetype="dashed", color = "red")
	or.plots[[j]] <- p
}

or.plots.comb <- plot_grid(or.plots[[1]],or.plots[[2]],or.plots[[3]],or.plots[[4]],labels=c("A","B","C","D"),nrow=2,ncol=2)
ggsave2("odds_ratio_plot.pdf",or.plots.comb)

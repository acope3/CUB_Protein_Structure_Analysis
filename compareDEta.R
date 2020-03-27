library(ggplot2)
library(viridis)
library(cowplot)
library(ggnewscale)

compareRankings <- function(runs.to.compare)
{
  dfs <- lapply(runs.to.compare,FUN=read.csv,header=T,stringsAsFactors=F)
  aa <- aminoAcids()
  for (a in aa)
  {
    if (a == "M" || a == "W" || a == "X") next
    for (i in 1:length(dfs))
    {
      index <- which(dfs[[i]][,"AA"] == a)
      posterior <- dfs[[i]][index,"Posterior"]
      rank.post <- rank(posterior)
      dfs[[i]][index,"Posterior"] <- rank.post
    }
  }
  print(dfs)
  agree <- which(dfs[[1]][,"Posterior"] == dfs[[2]][,"Posterior"] & dfs[[2]][,"Posterior"] == dfs[[3]][,"Posterior"])
  print(agree)
  return (length(agree)/nrow(dfs[[1]]))
}

compareAcrossCodon <- function(runs.to.compare,legend.info,colors=NULL,legend.title="Legend",...)
{
  if (length(runs.to.compare) != length(legend.info)) stop("runs.to.compare and legend.info should be same length")
  if (is.null(colors))
  {
    colors <- c("black","red","blue","green","purple")
  }
  cur.aa <- "A"
  labels <- c()
  splits <- c()
  
  placeholder <- data.frame(AA="",Codon="",Posterior=NA,Std.Dev=NA,X0.025.=NA,X0.975.=NA)
  dfs <- lapply(runs.to.compare,FUN=read.csv,header=T,stringsAsFactors=F)
  df.1 <- dfs[[1]]
  for (i in 1:nrow(df.1))
  {
    aa <- df.1[i,"AA"]
    codon <- df.1[i,"Codon"]
    if (aa != cur.aa)
    {
      splits <- c(splits,i)
      labels <- c(labels,"")
      cur.aa <- aa
    }
    labels <- c(labels,codon)
  }
  max.element <- max(unlist(lapply(dfs,FUN=function(x){max(x[,6])})))
  min.element <- min(unlist(lapply(dfs,FUN=function(x){min(x[,5])})))
  split.line <- which(labels == "")
  for (j in 1:length(dfs))
  {
    placeholder <- data.frame(AA="",Codon="",Posterior=NA,Std.Dev=NA,X0.025.=NA,X0.975.=NA)
    for (i in split.line)
    {
      dfs[[j]] <- rbind(dfs[[j]][1:i-1,],placeholder,dfs[[j]][i:nrow(dfs[[j]]),])
    }
  }
  ylim <- c(min.element - 0.05,max.element+0.05)
  
  par(mar=c(5.1,3.8,5.1,5.1))
  aa.spots <- c(0,split.line)+diff(c(0,split.line))/2
  plot(x=NULL,y=NULL,xlim=c(1,length(labels)),ylim=ylim,xlab="Codon",ylab=expression(Delta*eta),xaxt='n',...)
  for (i in 1:length(dfs))
  {
    points(x=1:length(labels),dfs[[i]][,3],col=colors[i])
    arrows(x0=1:length(labels),y0=dfs[[i]][,5],y1=dfs[[i]][,6],angle=90,length=0.05,col=colors[i],code=3)
  }
  abline(h=0,lty=2)
  axis(side=1,at=1:length(labels),labels=labels,las=2,cex.axis=0.75)
  axis(side=3,at=aa.spots,labels=unique(df.1[which(df.1[,"AA"]!=""),"AA"]),cex.axis=0.75)
  abline(v = split.line,lty=2,col="lightblue")
  legend(x = "topleft",inset=c(0.01,0.03),xpd=T,legend = legend.info,col = colors,pch=rep(1,length(dfs)),box.lty=0,title=legend.title)
  
}
getSignificantCodons<-function(dfs,category.names)
{
  numCat <- length(category.names)
  for (i in 1:numCat)
  {
    dfs[[i]][,"Significance"] <- rep("Not Significant",nrow(dfs[[i]])) 
  }
  for (i in 1:(numCat-1))
  {
    for (j in 2:numCat)
    {
      if (category.names[i] == category.names[j]) next
      df.1 <- dfs[[i]]
      df.2 <- dfs[[j]]
      sig <- which((df.1[,5] > df.2[,6]) | (df.2[,5] > df.1[,6]))
      for (k in sig)
      {
        if (df.1[k,"Significance"] == "Not Significant")
        {
          df.1[k,"Significance"] <- paste(category.names[i],category.names[j],sep="/")
        } else{
          df.1[k,"Significance"] <- paste(df.1[k,"Significance"],paste(category.names[i],category.names[j],sep="/"),sep=",")
        }
        
        if (df.2[k,"Significance"] == "Not Significant")
        {
          df.2[k,"Significance"] <- paste(category.names[i],category.names[j],sep="/")
        } else{
          df.2[k,"Significance"] <- paste(df.2[k,"Significance"],paste(category.names[i],category.names[j],sep="/"),sep=",")
        }
      }
      dfs[[i]] <- df.1
      dfs[[j]] <- df.2
    }
  }
  return(dfs)
}


highlightSignificant<-function(df)
{
  color.df <- data.frame(Color=c("blue","green","red","purple","purple","orange","orange","cyan","cyan",rep("brown",6),"blue"),row.names=c("Coil/Helix","Coil/Sheet","Helix/Sheet",
                                                                                                                                           "Coil/Helix,Coil/Sheet","Coil/Sheet,Coil/Helix","Coil/Helix,Helix/Sheet","Helix/Sheet,Coil/Helix","Coil/Sheet,Helix/Sheet","Helix/Sheet,Coil/Sheet",
                                                                                                                                           "Coil/Helix,Coil/Sheet,Helix/Sheet","Coil/Helix,Helix/Sheet,Coil/Sheet","Coil/Sheet,Coil/Helix,Helix/Sheet","Helix/Sheet,Coil/Helix,Coil/Sheet","Coil/Sheet,Helix/Sheet,Coil/Helix","Helix/Sheet,Coil/Sheet,Coil/Helix","Ordered/Disordered"),stringsAsFactors=F)
  # color.map <- c("blue","red","green","orange","purple")
  # significant.types <- unique(df$Significance)
  # significant.types <- significant.types[which(significant.types != "Not Significant")]
  # color.df <- data.frame(Color=color.map[1:length(significant.types)],row.names=significant.types,stringsAsFactors=F)
  df["Significance.Color"] <- rep("black",nrow(df))
  for (i in 1:nrow(df))
  {
    if (df[i,"Significance"] != "Not Significant")
    {
      tmp.color <- color.df[df[i,"Significance"],"Color"]
      df[i,"Significance.Color"] <- tmp.color
    }
  }
  return(df)
  
  
}

compareDeltaEta <- function(runs.to.compare,legend.info,ref.category=1,filter.aa=NULL,main="Differences in Selection",legend.title="Secondary Structure",legend.pos = c(0.1,0.85),...)
{
  dfs <- lapply(runs.to.compare,FUN=read.csv,header=T,stringsAsFactors=F)
  for (i in 1:length(dfs))
  {
    tmp <- dfs[[i]]
    tmp<- tmp[which(tmp$Posterior != 0),] 
    dfs[[i]] <- tmp
  }
  dfs <- getSignificantCodons(dfs,legend.info)
  ref.df <- dfs[[ref.category]]
  for (i in 1:length(dfs))
  {
    tmp <- dfs[[i]]
    tmp[,c("Posterior","X0.025.","X0.975.")] <- tmp[,c("Posterior","X0.025.","X0.975.")] - ref.df[,c("Posterior")] 
    dfs[[i]] <- tmp
  }
  for (i in 1:length(dfs))
  {
    dfs[[i]][,"Category"] <- legend.info[i]
  }
  
  df.test <- do.call("rbind",dfs)
  if (is.null(filter.aa)==F)
  {
    df.test <- df.test[which(df.test$AA %in% filter.aa),]
  }
  df.test[,"AA.Codon"] <- paste(df.test[,"AA"],df.test[,"Codon"],sep="/")
  for (i in 1:nrow(df.test))
  {
    aa.codon <- df.test[i,"AA.Codon"]
    rows <- which(df.test$AA.Codon == aa.codon)
    significance <- df.test[rows,"Significance"]
    tmp <- c()
    for (j in significance)
    {
      x<-unlist(strsplit(significance,split=","))
      tmp <- c(tmp,x)
    }
    tmp <- unique(tmp[which(tmp != "Not Significant")])
    if (length(tmp) > 0)
    {
      tmp <- paste(tmp,collapse=',')
      df.test[rows,"Significance"] <- tmp
    }else{
      df.test[rows,"Significance"] <- "Not Significant"
    }
  }
  codons <- df.test$Codon
  nuc <- strsplit(codons,split="")
  ct <- unlist(lapply(nuc,function(x){
    if (x[3] == "C" || x[3] == "T") return("CT")
    else return("AG")
    }))
  df.test["Third.Nucleotide"] <- ct
  df.test <- highlightSignificant(df.test)
  df.test[which(df.test$Significance != "Not Significant"),"Significance.Color"] <- "red" 
  uniqueInitials <- df.test$AA
  initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))
  p <- ggplot(df.test,aes(x=Codon,y=Posterior,label=AA))  + geom_hline(yintercept = 0,linetype="dashed")
  
  p <- p + geom_errorbar(aes(ymin=X0.025., ymax=X0.975.,linetype=Category),position=position_dodge(0.5),width=0,...) #
  #p <- p + geom_pointrange(aes(ymin=X0.025., ymax=X0.975.,color=Category,shape=uniqueInitials),position=position_dodge(0.5),...)
  #p <- p + scale_shape_manual(values=initialShapes)
  p <- p + scale_x_discrete(limits=unique(df.test[,"Codon"]))
  p <- p + scale_linetype_manual(values=c("solid","dotted","longdash"))
 # p <- p + scale_color_viridis(discrete=TRUE)+ theme_cowplot()
  
  p <- p + theme_cowplot()
  p <- p + labs(x="Codon",y=bquote(Delta*eta~"-"~Delta*eta[.(legend.info[ref.category])]),color=legend.title)
  
  
  p <- p + new_scale_color()
  
  p <- p + geom_text(aes(color=Third.Nucleotide),position=position_dodge2(0.5))
  p <- p + scale_shape_manual(values=initialShapes)
  p <- p + labs(color="3rd Base Nucleotide")
  p <- p  + scale_color_manual(values=c("blue","orange"))
  
  p <- p +ggtitle(main)
  p <- p + theme(axis.text.x = element_text(angle = 90,colour=df.test$Significance.Color))
  # p <- p + theme(axis.text.x = element_text(angle = 90,colour=df.test$Significance.Color),
  #                legend.position = legend.pos,
  #                legend.text = element_text(size=10,face="bold"),
  #                legend.title = element_text(size=12,face="bold"),
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + theme(legend.position = legend.pos,axis.title=element_text(size = 8,face="bold"),axis.text=element_text(size=10,face = "bold",color = "black"),plot.title = element_text(hjust = 0.5,size=10,face="bold",color="black"))
  #p <- p + theme(plot.title = element_text(hjust = 0.5,size=8,face="bold",color="black"),legend.position = legend.pos,)
  
  
  return(p)
}

runs.to.compare = c("Scer/Predicted/Results/Secondary_structures/Coil/final_run/Parameter_est/coil_Selection.csv",
                     "Scer/Predicted/Results/Secondary_structures/Helix/final_run/Parameter_est/helix_Selection.csv",
                    "Scer/Predicted/Results/Secondary_structures/Sheet/final_run/Parameter_est/sheet_Selection.csv")
legend.info=c("Coil","Helix","Sheet")
ref.category <- 1
p1<-compareDeltaEta(runs.to.compare = runs.to.compare,
                    legend.info = legend.info,
                    ref.category = ref.category,
                    filter.aa=c("C","D","E","F","H","N","K","Q","Y","Z"),
                    main=expression("2-Codon Amino Acids"),
                    legend.pos=c(0.1,0.85),size=0.25)
p1a <- p1+theme(legend.position='none')

ggsave(plot=p1,filename="Images_for_mike/two_codon_aa_dEta_comparison_example.pdf",device="pdf")
# p1a <- p1a + theme(axis.title.x = element_blank())

# p2<-compareDeltaEta(runs.to.compare = runs.to.compare,
#                     legend.info = legend.info,
#                     ref.category = ref.category,
#                     filter.aa=c("N","Q","S","Z","T"),
#                     main=expression("Polar Uncharged"),
#                     legend.pos='none',size=0.25)

# ggsave(plot=p2,filename="Images_for_mike/polar_uncharged_aa.pdf",device="pdf")
# p2 <- p2 + theme(axis.title.x = element_blank())

# p3 <- compareDeltaEta(runs.to.compare = runs.to.compare,
#                       legend.info = legend.info,
#                       ref.category = ref.category,
#                       filter.aa=c("A","F","I","L","V","Y"),
#                       main=expression("Hydrophobic"),
#                       legend.pos='none',size=0.25)
# ggsave(plot=p3,filename="Images_for_mike/hydro_aa.pdf",device="pdf")
# p3 <- p3 + theme(axis.title.x = element_blank())

# p4 <- compareDeltaEta(runs.to.compare = runs.to.compare,
#                       legend.info = legend.info,
#                       ref.category = ref.category,
#                       filter.aa=c("C","G","P"),
#                       main=expression("Other"),
#                       legend.pos='none',size=0.25)
# ggsave(plot=p4,filename="Images_for_mike/other_aa.pdf",device="pdf")


# cat.legend <- get_legend(p1)

# p5 <- plot_grid(p1a,p2,p3,p4,labels = c('A', 'B','C','D'), label_size = 12,align='v',nrow=4,axis='l')

# dummy.df <- data.frame(Significance=c("Not Significant","Coil is different from Helix",
#                             "Coil is different from Sheet",
#                             "Helix is different from Sheet",
#                             "Coil is different from Helix and Sheet",
#                             "Coil is different from Helix,\nand Helix is different from Sheet",
#                             "Coil is different from Sheet,\nand Sheet is different from Helix",
#                             "All structures are different"),Y=1:8,stringsAsFactors=F)
# cols <- c("Coil is different from Helix"="blue",
#                             "Coil is different from Sheet"="green",
#                             "Helix is different from Sheet"="red",
#                             "Coil is different from Helix and Sheet"="purple",
#                             "Coil is different from Helix,\nand Helix is different from Sheet"="orange",
#                             "Coil is different from Sheet,\nand Sheet is different from Helix"="cyan",
#                             "All structures are different"="brown","Not Significant"="black")
# # print(cols)
# dummy.df$Significance <- factor(dummy.df$Significance, levels = c("Not Significant","Coil is different from Helix",
#                             "Coil is different from Sheet",
#                             "Helix is different from Sheet",
#                             "Coil is different from Helix and Sheet",
#                             "Coil is different from Helix,\nand Helix is different from Sheet",
#                             "Coil is different from Sheet,\nand Sheet is different from Helix",
#                             "All structures are different"))
# dummy.leg.ord <- levels(dummy.df$Significance)
# dummyplot <- (ggplot(data=dummy.df,aes(x=Significance,y=Y,fill=Significance)) 
#                     + geom_bar(stat="identity") 
#                     + scale_fill_manual(values=cols)
#                     + theme(legend.title=element_text(size=8,face="bold"), legend.text = element_text(size=6,face="bold")))
# dummy.legend <- get_legend(dummyplot)
# blank_p <- plot_spacer() + theme_void()
# p6 <- plot_grid(cat.legend,dummy.legend,blank_p,nrow=2,ncol=1)
# p7 <- plot_grid(p5,p6,nrow=1,ncol=2,rel_widths=c(2.5,1))
# ggsave2(plot=p7,filename="Images_for_mike/ss_lattice.pdf",device="pdf",height=10,width=9)







# runs.to.compare <- c("Scer/Predicted/Results/Ordered_disordered_unconserved/Ordered/final_run/Parameter_est/ordered_Selection.csv",
#                     "Scer/Predicted/Results/Ordered_disordered_unconserved/Disordered/final_run/Parameter_est/disordered_Selection.csv")
# legend.info<-c("Ordered","Disordered")
# ref.category <- 2
# p1<-compareDeltaEta(runs.to.compare = runs.to.compare,
#                     legend.info = legend.info,
#                     ref.category = ref.category,
#                     filter.aa=c("D","E","H","K","R"),
#                     main=expression(atop("Overall Orderedness "*Delta*eta*":","Charged Amino Acids")),
#                     legend.pos=c(0.2,0.85))

# ggsave(plot=p1,filename="Images_for_mike/charged_aa_order_unconserved.pdf",device="pdf")

# p2<-compareDeltaEta(runs.to.compare = runs.to.compare,
#                     legend.info = legend.info,
#                     ref.category = ref.category,
#                     filter.aa=c("N","Q","S","Z","T"),
#                     main=expression(atop("Overall Orderedness "*Delta*eta*":","Polar Uncharged Amino Acids")),
#                     legend.pos=c(0.2,0.85))

# ggsave(plot=p2,filename="Images_for_mike/polar_charged_aa_order_unconserved.pdf",device="pdf")

# p3 <- compareDeltaEta(runs.to.compare = runs.to.compare,
#                       legend.info = legend.info,
#                       ref.category = ref.category,
#                       filter.aa=c("A","F","I","L","V","Y"),
#                       main=expression(atop("Overall Orderedness "*Delta*eta*":","Hydrophobic Amino Acids")),
#                       legend.pos=c(0.2,0.7),size=0.5,fatten=3)
# ggsave(plot=p3,filename="Images_for_mike/hydro_aa_order_unconserved.pdf",device="pdf")
# p4 <- compareDeltaEta(runs.to.compare = runs.to.compare,
#                       legend.info = legend.info,
#                       ref.category = ref.category,
#                       filter.aa=c("C","G","P"),
#                       main=expression(atop("Overall Orderedness "*Delta*eta*":","Other Amino Acids")),
#                       legend.pos=c(0.8,0.85))
# ggsave(plot=p4,filename="Images_for_mike/other_aa_order_unconserved.pdf",device="pdf")

#pdf("scer_compare_liberal_homology_deta_helix_termini.pdf",width=10)
# 
# compareAcrossCodon(runs.to.compare =c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/Turn_Coil/final_run/Parameter_est/turn_coil_Selection.csv",
#                                       "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/Helix_Sheet/final_run/Parameter_est/helix_sheet_Selection.csv"),
#                   legend.info = c("Turn/Coil","Helix/Sheet"),main=expression(Delta*eta*" Across Codons: Experimental Secondary Structures"))
# 
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Secondary_structures/Coil/final_run/Parameter_est/coil_Selection.csv",
#                                       "Scer/Predicted/Results/Secondary_structures/Helix/final_run/Parameter_est/helix_Selection.csv",
#                                       "Scer/Predicted/Results/Secondary_structures/Sheet/final_run/Parameter_est/sheet_Selection.csv"),
#                   legend.info=c("Coil","Helix","Sheet"),main=expression(Delta*eta*" Across Codons: Predicted Secondary Structures"),legend.title="Secondary Structure")

# compareAcrossCodon(runs.to.compare = c("Scer/Exp_liberal_homology/Results/Secondary_structures/Coil/final_run/Parameter_est/coil_Selection.csv",
#                                        "Scer/Exp_liberal_homology/Results/Secondary_structures/Turn/final_run/Parameter_est/turn_Selection.csv",
#                                        "Scer/Exp_liberal_homology/Results/Secondary_structures/Helix/final_run/Parameter_est/helix_Selection.csv",
#                                       "Scer/Exp_liberal_homology/Results/Secondary_structures/Sheet/final_run/Parameter_est/sheet_Selection.csv"),
#                    legend.info=c("Coil","Turn","Helix","Sheet"),main=expression(Delta*eta*" Across Codons: Experimental Secondary Structures"))

# pdf("Images_for_mike/ordered_disordered_conserved_vs_nonconserved_sites_deta_compare_across_codons.pdf",width=10)
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Ordered_disordered_conserved/Disordered/final_run/Parameter_est/disordered_Selection.csv",
#                                        "Scer/Predicted/Results/Ordered_disordered_unconserved/Disordered/final_run/Parameter_est/disordered_Selection.csv"),
#                    legend.info=c("Disordered (Conserved)","Disordered (Unconserved)"),main=expression(Delta*eta*" Across Codons: Conserved vs. Unconserved sites in Disordered Regions"),legend.title="Structure")
# 
# 
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Ordered_disordered_conserved/Ordered/final_run/Parameter_est/ordered_Selection.csv",
#                                        "Scer/Predicted/Results/Ordered_disordered_unconserved/Ordered/final_run/Parameter_est/ordered_Selection.csv"),
#                    legend.info=c("Ordered (Conserved)","Ordered (Unconserved)"),main=expression(Delta*eta*" Across Codons: Conserved vs. Unconserved sites in Ordered Regions"),legend.title="Structure")
# 
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Ordered_disordered_conserved/Ordered/final_run/Parameter_est/ordered_Selection.csv",
#                                        "Scer/Predicted/Results/Ordered_disordered_conserved/Disordered/final_run/Parameter_est/disordered_Selection.csv"),
#                    legend.info=c("Ordered (Conserved)","Disordered (Conserved)"),main=expression(Delta*eta*" Across Codons at Conserved Residues: Ordered vs. Disordered"),legend.title="Structure")
# 
# 
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Ordered_disordered_unconserved/Ordered/final_run/Parameter_est/ordered_Selection.csv",
#                                        "Scer/Predicted/Results/Ordered_disordered_unconserved/Disordered/final_run/Parameter_est/disordered_Selection.csv"),
#                    legend.info=c("Ordered (Non-conserved)","Disordered (Non-conserved)"),main=expression(Delta*eta*" Across Codons at Non-conserved Residues: Ordered vs. Disordered"),legend.title="Structure")
# 
# compareAcrossCodon(runs.to.compare = c("Scer/Predicted/Results/Ordered_disordered_conserved/Ordered/final_run/Parameter_est/ordered_Selection.csv",
#                                        "Scer/Predicted/Results/Ordered_disordered_conserved/Disordered/final_run/Parameter_est/disordered_Selection.csv",
#                                        
#
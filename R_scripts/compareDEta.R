## A plot for comparing \Delta\Etas across runs
## Note the input for these functions should be generated using the rescaleDEta.R script.

library(AnaCoDa)
library(deming)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(ggrepel)


colors <- c('"Not"~"Significant"'="black",
            "A"="blue",
            "D"="skyblue",
            "F"="purple",
            "G"="green",
            "K"="greenyellow",
            "L"="green4",
            "P"="yellow",
            "R"="orange",
            "S[4]"="red",
            "V"="cyan",
            "S[2]"="brown")



getSignificantCodons<-function(df.1,df.2)
{
  df.1[,"Significance"] <- rep("Not Significant",nrow(df.1)) 
  df.2[,"Significance"] <- rep("Not Significant",nrow(df.2)) 
  sig <- which((df.1[,5] > df.2[,6]) | (df.2[,5] > df.1[,6]))
  df.1[sig,"Significance"] <- "Significant"
  df.2[sig,"Significance"] <- "Significant"
  return(list(df.1,df.2))
}



demingRegression <- function(directory.1,directory.2,file.1,file.2,include.AA=c())
{
 
  sel.1 <- read.table(paste0(directory.1,"Parameter_est/",file.1),sep=",",header=TRUE,stringsAsFactors=F)
  sel.2 <- read.table(paste0(directory.2,"Parameter_est/",file.2),sep=",",header=TRUE,stringsAsFactors=F) 
  rownames(sel.1) <- sel.1[,2]
  rownames(sel.2) <- sel.2[,2]
  if (length(include.AA) != 0)
  {
    sel.1 <- sel.1[c(which(sel.1$AA %in% include.AA)),]
    sel.2 <- sel.2[c(which(sel.2$AA %in% include.AA)),]
  } 
  sd.1 <- sel.1$Std.Dev
  sd.2 <- sel.2$Std.Dev
  dfs <- getSignificantCodons(sel.1,sel.2)
  sel.1 <- dfs[[1]]
  sel.2 <- dfs[[2]]
  df <- plyr::join(sel.1,sel.2,by=c("AA","Codon","Significance"))
  colnames(df)[8] <- "Posterior.2"
  colnames(df)[9] <- "Std.Dev.2"
  colnames(df)[10] <- "X0.025.2"
  colnames(df)[11] <- "X0.975.2"

  ## Note that we regress through the origin (i.e. the y-intercept is 0), which is the reason for the "+ 0"
  ## This is justified because by rescaling \Delta\Eta by the mean, the average value for \Delta\Eta is 0 
  reg <- deming(Posterior.2 ~ Posterior+0,data=df,xstd=sd.1,ystd=sd.2)
  b1 <- reg$coefficients[2]
  ci.b1 <- reg$ci[2,]
  p.val <- pnorm(abs((b1-1)/sqrt(reg$variance[2,2])),lower.tail = F)*2
  return(list(df=df,Slope=b1,Slope.CI=ci.b1,SD.1=sd.1,SD.2=sd.2,p.val=p.val,std.err = sqrt(reg$variance[2,2])))
}

changeCodon <- function(codon)
{
  nucleotides <- unlist(strsplit(codon,""))
  if (nucleotides[3] == "A")
  {
    codon <- deparse(paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|G",collapse=""))
  } else if (nucleotides[3] == "C"){
    codon <- deparse(paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|T",collapse=""))
  }
  return(codon)
}



plotDeming<-function(data,b1,b0,categories=c("X","Y"),ci = F,bounds=NULL,file="title",title="Regression",range.xy = NULL)
{
  if (!is.null(bounds))
  {
    conf.int.1 <- bounds[1]
    conf.int.2 <- bounds[2]
    l <- data.frame(s=c(b1,conf.int.2,conf.int.1, 1.0),ic=c(b0,0.0,0.0,0.0),Line=c("Deming Slope","97.5% CI","2.5% CI","1:1 Line"),stringsAsFactors = F)
    l$Line <- factor(l$Line,levels=c("Deming Slope","97.5% CI","2.5% CI","1:1 Line"))
    levels(l$Line) <- c(paste0("Deming Slope: ",round(b1,3)),paste0("95% CI: ",round(conf.int.1,3),"-",round(conf.int.2,3)),paste0("95% CI: ",round(conf.int.1,3),"-",round(conf.int.2,3)),"1:1 Line")
    legend.colors <- c("black","black","red")
    lines.reg <- c("solid","dashed","solid")
    
  } else{
    l <- data.frame(s=c(b1,1.0),ic=c(b0,0.0),Line=c("Model II Regression","y = x"))
    legend.colors <- c("black","black")
    lines.reg <- c("solid","dashed")
  }
  p <- ggplot(data,aes(Posterior,Posterior.2))
  p <-(p + geom_abline(data=l,mapping=aes(slope=s,intercept=ic,linetype=Line,color=Line))
       + scale_color_manual(values=legend.colors,guide = guide_legend(order = 1))
       + scale_linetype_manual(values=lines.reg,guide = guide_legend(order = 1)))
  p <- p + labs(linetype="Deming Regression",color="Deming Regression")
  
  if (ci)
  {
    if (is.null(range.xy))
    {
      range.xy <- range(c(data[,5:6],data[,10:11]),na.rm=T)
    }
    xlim <- range.xy
    ylim <- range.xy
    p <- p + scale_x_continuous(limits = range.xy+c(-0.05,0.05)) + scale_y_continuous(limits = range.xy+c(-0.05,0.05))
    p <- (p + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),color="black") 
          + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),color="black"))
  } else{
    range.xy <- range(c(data[,3],data[,7]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  }
  p<-(p + new_scale_color()
       + geom_point(size=5,alpha=0.6)
       + labs(x=bquote(.(categories[1])~"("*Delta*eta*")"),y=bquote(.(categories[2])~"("*Delta*eta*")"))
       + ggtitle(label=title))
  p <- p + guides(color = T,linetype=F)
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.2,label=deparse(cor.exp),parse=T,fontface="bold",size=8)
  p <- (p + theme_bw()
       + theme(axis.title=element_text(face="bold",size=18),axis.text=element_text(face="bold",size=12))
       + theme(axis.line = element_line(colour = "black"))
       + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
       + theme(legend.position = c(0.25,0.8),legend.text=element_text(face="bold",size=14),legend.title=element_text(face="bold",size=14),legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
       + theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")))

  ggsave(filename = file,plot=p,device="pdf",width = 7,height = 7)
  return(p)
}


plotDetaWODeming<-function(data,p.val,categories=c("X","Y"),ci = F,file="title",title="Regression",range.xy = NULL)
{
  data[which(data$AA == "S"),"AA"] <- "S[4]"
  data[which(data$AA == "Z"),"AA"] <- "S[2]"
  data[,"Codon"] <- unlist(lapply(data$Codon,changeCodon))
 
  l <- data.frame(s=c(1.0),ic=c(0.0),Line=c("1:1 Line"),stringsAsFactors = F)
  l$Line <- factor(l$Line,levels=c("1:1 Line"))
  levels(l$Line) <- c("1:1 Line")
  legend.colors <- c("red")
  lines.reg <- c("solid")
    
  
  
  uniqueInitials <- paste(data$AA,data$Codon,sep=":")
  
  sig <- which(data$Significance == "Significant")
  #uniqueInitials[sig] <- paste0(uniqueInitials[sig],'*"*"',sep="")
  uniqueInitials[-sig] <- NA
  data[-sig,"AA"] <- '"Not"~"Significant"'
  p <- ggplot(data,aes(x=Posterior,y=Posterior.2))
  p <-(p + geom_abline(data=l,mapping=aes(slope=s,intercept=ic,linetype=Line,color=Line),show.legend=F)
       + scale_colour_manual(values=legend.colors,guide = guide_legend(order = 1),drop=F)
       + scale_linetype_manual(values=lines.reg,guide = guide_legend(order = 1)))
  if (ci)
  {
    if (is.null(range.xy))
    {
      range.xy <- range(c(data[,5:6],data[,10:11]),na.rm=T)
    }
    if (range.xy[1] > 0)
    {
      range.xy[1] <- -0.05
    }
    if (range.xy[2] < 0)
    {
      range.xy[2] <- 0.05
    }
    xlim <- range.xy
    ylim <- range.xy
    p <- p + scale_x_continuous(limits = range.xy+c(-0.005,0.005)) + scale_y_continuous(limits = range.xy+c(-0.005,0.005))
    
  } else{
    range.xy <- range(c(data[,3],data[,7]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  }
  
  p <- p + geom_hline(yintercept=0,color="red",linetype="dashed") + geom_vline(xintercept=0,color="red",linetype="dashed")
  p<-(p + new_scale_color()
      + geom_point(aes(color=AA),size=6,alpha=0.6)
      #+ geom_point(size=4,alpha=0.6)
      + geom_text_repel(label=uniqueInitials,size=6,box.padding = 1,point.padding=0.1,parse=T,segment.alpha=0.7,fontface="bold")
      #+ labs(x=bquote("Selection("~Delta*eta*"):"~.(categories[1])),y=bquote("Selection("~Delta*eta*"):"~.(categories[2])),color="3rd Position") 
      #+ labs(x=bquote(atop(.(categories[1]),atop("Selection Coefficients: G or T","for A/C "%<-%" Increasing Selection "%->%" for G/T"))),y=bquote(atop(.(categories[2]),atop("Selection Coefficients: G or T","for A/C"%<-%" Increasing Selection "%->%"for G/T"))),color="",parse=T)
      + labs(x=bquote(atop(.(categories[1])~"("*Delta*eta["AC|GT"]*")","for A/C "%<-%" Increasing Selection "%->%" for G/T")),y=bquote(atop(.(categories[2])~"("*Delta*eta["AC|GT"]*")","for A/C"%<-%" Increasing Selection "%->%"for G/T")),color="",parse=T)
      
      #+ scale_colour_manual(values=c("Pur-Pur"="blue","Pyr-Pyr"="orange","Pur-Pyr"="green","Pyr-Pur"="yellow"),guide = guide_legend(order = 2),drop=F)
      #+ scale_colour_discrete(breaks = c(unique(data$AA[sig]),labels=lapply(unique(data$AA[sig]),function(x)parse(text=x)),guide=guide_legend(nrow=ifelse(length(unique(data$AA[sig])) > 5,5,length(unique(data$AA[sig]))),ncol=ifelse(length(unique(data$AA[sig])) > 5,3,1)))
      #+ scale_colour_discrete(breaks = unique(data$AA[sig]),labels=lapply(unique(data$AA[sig]),function(x)parse(text=x)),guide=guide_legend(nrow=ifelse(length(unique(data$AA)) > 5,5,length(unique(data$AA))),ncol=ifelse(length(unique(data$AA)) > 5,3,1)))
      + scale_colour_manual(values=colors[unique(data$AA)],
        breaks=names(colors[unique(data$AA[sig])]),
        labels=lapply(names(colors[unique(data$AA[sig])]),function(x)parse(text=x)),
        guide=guide_legend(nrow=ifelse(length(unique(data$AA)) > 5,5,length(unique(data$AA))),ncol=ifelse(length(unique(data$AA)) > 5,3,1)))
      
      + ggtitle(label=bquote(.(title)))
      )

  p <- (p + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),color="black") 
          + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),color="black"))
  #p <- p + guides(color = guide_legend(override.aes = list(label="",size=1)))
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.2,label=deparse(cor.exp),parse=T,fontface="bold",size=6)
  #p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.10,label="*p < 0.05",parse=F,size=8)
  p <- (p + theme_bw()
        + theme(axis.title=element_text(face="bold",size=16),axis.text=element_text(face="bold",size=12))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.15,0.85),legend.text=element_text(face="bold",size=16),legend.title=element_text(face="bold",size=14),legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
        + theme(plot.title = element_text(hjust = 0.5,size=18,face="bold")))
  ggsave(filename = file,plot=p,device="pdf",width = 7,height = 7)
  return(p)
 
}



plotDemingRegressionSecondaryStructuresAndOrdered <- function(head.directory)
{
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p5 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Coil","Disordered Coil"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil (Structured vs. Disordered)",file="../Images_for_dissertation/scer_coil_ordered_vs_disordered.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p6 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Helix","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nHelix (Structured vs. Disordered)",file="../Images_for_dissertation/scer_helix_ordered_vs_disordered.pdf")


  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p7 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil vs. Helix (Structured Regions)",file="../Images_for_dissertation/scer_ordered_coil_vs_ordered_helix.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p8 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil vs. Sheet (Structured Regions)",file="../Images_for_dissertation/scer_ordered_coil_vs_ordered_sheet.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p9 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nHelix vs. Sheet (Structured Regions)",file="../Images_for_dissertation/scer_ordered_helix_vs_ordered_sheet.pdf")


  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p10 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Disordered Helix",file="../Images_for_dissertation/scer_disordered_helix_vs_disordered_coil_disordered.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p11 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Sheet",file="../Images_for_dissertation/scer_disordered_coil_vs_sheet.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p12 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Helix vs. Sheet",file="../Images_for_dissertation/scer_disordered helix_vs_sheet.pdf")

}


plotDemingRegressionOrdered <- function(head.directory)
{
  output <- demingRegression(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p4 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nStructured vs. Disordered Regions",file="../Images_for_dissertation/scer_ordered_vs_disorderd.pdf")

  ## Ordered and Disordered conserved vs non-conserved
  output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p.cons <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured vs. Disordered (Conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_vs_disorderd_conserved_sacch.pdf"))

  output <- demingRegression(file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p.non <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured vs Disordered (Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_vs_disorderd_nonconserved_sacch.pdf"))


  output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p.ord <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Residues","Non-conserved Residues"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured Regions (Conserved vs Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_conserved_vs_nonconserved_sacch.pdf"))

  output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p.dis <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Residues","Non-conserved Residues"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nDisordered Regions (Conserved vs Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_disorderd_conserved_vs_nonconserved_sacch.pdf"))


}


plotDemingRegressionSecondaryStructures <- function(head.directory)
{


  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p1 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_pres/scer_coil_vs_helix_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p2 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_pres/scer_coil_vs_sheet_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p3 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_pres/scer_helix_vs_sheet_already_normalized.pdf")

  
  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Coil","final_run/"),file.path(head.directory,"Secondary_structures_conserved_no_ncast","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p4 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_pres/scer_cons_coil_vs_helix_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Coil","final_run/"),file.path(head.directory,"Secondary_structures_conserved_no_ncast","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p5 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_pres/scer_cons_coil_vs_sheet_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Helix","final_run/"),file.path(head.directory,"Secondary_structures_conserved_no_ncast","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p6 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_pres/scer_cons_helix_vs_sheet_already_normalized.pdf")


  output <- demingRegression(file.path(head.directory,"Secondary_structures_variable_no_ncast","Coil","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p7 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_pres/scer_sim_var_coil_vs_helix_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_variable_no_ncast","Coil","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p8 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_pres/scer_sim_var_coil_vs_sheet_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_variable_no_ncast","Helix","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p9 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_pres/scer_sim_var_helix_vs_sheet_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Coil","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Coil","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p4 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Sites","Variable Sites"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil",file="../Images_for_pres/scer_coil_cons_var_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Sheet","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p5 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Sites","Variable Sites"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Sheet",file="../Images_for_pres/scer_sheet_cons_var_already_normalized.pdf")

  output <- demingRegression(file.path(head.directory,"Secondary_structures_conserved_no_ncast","Helix","final_run/"),file.path(head.directory,"Secondary_structures_variable_no_ncast","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  p6 <- plotDeming(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Sites","Variable Sites"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix",file="../Images_for_pres/scer_helix_cons_var_already_normalized.pdf")


  ## Return these plots so they can be combined into one plot using cowplot package
  return(list(p1,p2,p3))

}




plotDEtaByAA <- function(head.directory,target.directory,aa.groups=list(c("N","Q","S","T","C","Z"), 
                                                       c("D","E","H","K","R"),
                                                       c("A","F","G","I","L","P","V","Y")),
                                        groups=c("Uncharged polar","Charaged","Hydrophobic"))
{

  plots <- vector(mode="list",length=16)
  for (i in 1:16)
  {
    plots[[i]] <- vector(mode="list",length=3)
  }
  for (i in 1:length(aa.groups))
  {
    codon.group <- groups[i]
    ## Structures
    output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Coil","Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_coil_vs_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[1]][[i]] <- p
    

    output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[2]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[3]][[i]] <- p
    
    
    output <- demingRegression(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Region","Disordered Regions"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[4]][[i]] <- p
    
    ## Ordered Conserved vs Non-conserved
    output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Region (Conserved Residues)","Structured Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[5]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Disordered Region (Conserved Residues)","Disordered Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[6]][[i]] <- p

    ## Ordered and Disordered conserved vs non-conserved
    output <- demingRegression(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Region (Conserved Residues)","Disordered Region (Conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_conserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[7]][[i]] <- p

    output <- demingRegression(file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Region (Non-conserved Residues)","Disordered Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[8]][[i]] <- p

    # Ordered vs Disorderd Secondary Structures (Coil and Helix Only)
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Coil","Disordered Coil"),ci=T,title =bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_disordered_coil_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[9]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Structured Helix","Disordered Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_helix_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[10]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Coil","Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_ordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[11]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_ordered_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[12]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Disordered Coil","Disordered Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_coil_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[13]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Disordered Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[14]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_helix_vs_ordered_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[15]][[i]] <- p
    
    output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",include.AA=aa.groups[[i]])
    p <- plotDetaWODeming(output$df,categories = c("Disordered Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[16]][[i]] <- p

    

  }
  return(plots)
}


# aa.groups <- list(c("C","D","E","F","H","I","K","N","Q","Y","Z"),
#                   c("A","G","P","S","T","V"),
#                   c("R","L"))

#plots <- plotDemingRegressionSecondaryStructures(head.directory="../Scer/Predicted/Results/")

#plots <- plotDEtaByAA(head.directory="../Scer/Predicted/Results/",target.directory="../Images_for_pres/")

head.directory <- "../Scer/Predicted/Results/"
target.directory <- "../Images_for_pres/"
output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
p1 <- plotDetaWODeming(output$df,categories = c("Coil","Helix"),ci=T,title = bquote("Comparison of Selection  "*Delta*eta["(AC|GT)"]),file=paste0(target.directory,"/scer_coil_vs_helix_agct_codon.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
p2 <- plotDetaWODeming(output$df,categories = c("Coil","Sheet"),ci=T,title = bquote("Comparison of Selection  "*Delta*eta["(AC|GT)"]),file=paste0(target.directory,"/scer_coil_vs_sheet_agct_codon.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
p3 <- plotDetaWODeming(output$df,categories = c("Helix","Sheet"),ci=T,title = bquote("Comparison of Selection  "*Delta*eta["(AC|GT)"]),file=paste0(target.directory,"/scer_helix_vs_sheet_agct_codon.pdf"),range.xy=NULL)


#comb_top <- plot_grid(plots[[1]],plots[[2]],plots[[3]],labels = c("A","B","C"),ncol=3)

comb_top <- plot_grid(NULL,p1,labels = c("",""),ncol=2)
comb_bottom <- plot_grid(p2,p3,labels = c("",""),ncol=2)

comb <- plot_grid(comb_top,comb_bottom,nrow=2,ncol=1,align="rl")

ggsave2("../Images_for_pres/sec_struct_comp_agct.pdf",plot=comb,dpi=150,width=14,height=14)



# ##### Secondary Structures ###############################
# subplots <- plots[[1]]
# # legend <- get_legend(subplots[[1]])
# coil.helix <- plot_grid(subplots[[1]]+xlab("Coil")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain"))+ggtitle("Uncharged Polar\nAmino Acids"), 
#                       subplots[[2]] +xlab("Coil") + ylab("\n")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain"))+ ggtitle("Charged\nAmino Acids"), 
#                       subplots[[3]] +xlab("Coil")+ ylab("\n") + theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("Non-polar\nAmino Acids"),labels = c("A","B","C"),ncol=3,label_size=16)
# subplots <- plots[[2]]
# coil.sheet <- plot_grid(subplots[[1]] +xlab("Coil")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain")), 
#                         subplots[[2]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.1,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[3]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.15,0.85),axis.title.x = element_text(face="plain")),labels = c("D","E","F"),ncol=3,label_size=16)
# subplots <- plots[[3]]
# helix.sheet <- plot_grid(subplots[[1]] +theme(legend.position=c(0.2,0.75)), 
#                           subplots[[2]] + ylab("\n")+ theme(legend.position=c(0.1,0.8)), 
#                           subplots[[3]]+ ylab("\n") + theme(legend.position=c(0.15,0.8)),labels = c("G","H","I"),ncol=3,label_size=16)


# subplots <- plots[[11]]
# # legend <- get_legend(subplots[[1]])
# coil.helix <- plot_grid(subplots[[1]]+xlab("Coil")+theme(legend.position=c(0.2,0.8),axis.title.x = element_text(face="plain"))+ggtitle("2 Codon\nAmino Acids"), 
#                       subplots[[2]] +xlab("Coil") + ylab("\n")+theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("4 Codon\nAmino Acids"), 
#                       subplots[[3]] +xlab("Coil")+ ylab("\n") + theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("6 Codon\nAmino Acids"),labels = c("A","B","C"),ncol=3,label_size=16)
# subplots <- plots[[12]]
# coil.sheet <- plot_grid(subplots[[1]] +xlab("Coil")+theme(legend.position=c(0.2,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[2]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[3]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.15,0.85),axis.title.x = element_text(face="plain")),labels = c("D","E","F"),ncol=3,label_size=16)
# subplots <- plots[[15]]
# helix.sheet <- plot_grid(subplots[[1]] +theme(legend.position=c(0.2,0.8)), 
#                           subplots[[2]] + ylab("\n")+ theme(legend.position=c(0.15,0.8)), 
#                           subplots[[3]]+ ylab("\n") + theme(legend.position=c(0.15,0.8)),labels = c("G","H","I"),ncol=3,label_size=16)


# title <- ggdraw() + draw_label(
#     "Comparison of Selection Coefficients: Secondary Structures ",
#       fontface = 'bold',
#       x = 0,
#       hjust =0,size=25
#     ) +
#     theme(
#       # add margin on the left of the drawing canvas,
#       # so title is aligned with left edge of first plot
#       plot.margin = margin(0, 0, 0, 2)
#     )

# all.plots <- plot_grid(
#     title, coil.helix, coil.sheet,helix.sheet,
#     ncol = 1,
#     # rel_heights values control vertical title margins
#     rel_heights=c(0.1,1,1,1)
#   )
# ggsave2(file=paste0("../Images_for_dissertation/Examples/SS_ordered_by_aa.pdf"),plot=all.plots,width=18,height=16)



# # # ##### Ordered Disordered ###############################
# subplots <- plots[[4]]
# # legend <- get_legend(subplots[[1]])
# coil.helix <- plot_grid(subplots[[1]]+xlab("Structured Region")+theme(legend.position=c(0.25,0.8),axis.title.x = element_text(face="plain"))+ggtitle("Uncharged Polar\nAmino Acids"), 
#                       subplots[[2]] +xlab("Structured Region") + ylab("\n")+theme(legend.position=c(0.1,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("Charged\nAmino Acids"), 
#                       subplots[[3]] +xlab("Structured Region")+ ylab("\n") + theme(legend.position=c(0.20,0.82),axis.title.x = element_text(face="plain"))+ ggtitle("Non-polar\nAmino Acids"),labels = c("A","B","C"),ncol=3,label_size=16)
# subplots <- plots[[7]]
# coil.sheet <- plot_grid(subplots[[1]] +xlab("Structured Region (Conserved Residues)")+theme(legend.position=c(0.25,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[2]] +xlab("Structured Region (Conserved Residues)") + ylab("\n") +theme(legend.position=c(0.1,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[3]] +xlab("Structured Region (Conserved Residues)") + ylab("\n") +theme(legend.position=c(0.20,0.8),axis.title.x = element_text(face="plain")),labels = c("D","E","F"),ncol=3,label_size=16)
# subplots <- plots[[8]]
# helix.sheet <- plot_grid(subplots[[1]] +theme(legend.position=c(0.17,0.8),axis.title.y=element_text(size=12)), 
#                           subplots[[2]] + ylab("\n")+ theme(legend.position=c(0.1,0.8)), 
#                           subplots[[3]]+ ylab("\n") + theme(legend.position=c(0.20,0.82)),labels = c("G","H","I"),ncol=3,label_size=16)



# subplots <- plots[[4]]
# # legend <- get_legend(subplots[[1]])
# coil.helix <- plot_grid(subplots[[1]]+xlab("Structured Region")+theme(legend.position=c(0.2,0.8),axis.title.x = element_text(face="plain"))+ggtitle("2 Codon\nAmino Acids"), 
#                       subplots[[2]] +xlab("Structured Region") + ylab("\n")+theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("4 Codon\nAmino Acids"), 
#                       subplots[[3]] +xlab("Structured Region")+ ylab("\n") + theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("6 Codon\nAmino Acids"),labels = c("A","B","C"),ncol=3,label_size=16)
# subplots <- plots[[7]]
# coil.sheet <- plot_grid(subplots[[1]] +xlab("Structured Region (Conserved Residues)")+theme(legend.position=c(0.2,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[2]] +xlab("Structured Region (Conserved Residues)") + ylab("\n") +theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain")), 
#                         subplots[[3]] +xlab("Structured Region (Conserved Residues)") + ylab("\n") +theme(legend.position=c(0.15,0.85),axis.title.x = element_text(face="plain")),labels = c("D","E","F"),ncol=3,label_size=16)
# subplots <- plots[[8]]
# helix.sheet <- plot_grid(subplots[[1]] +theme(legend.position=c(0.2,0.8)), 
#                           subplots[[2]] + ylab("\n")+ theme(legend.position=c(0.15,0.8)), 
#                           subplots[[3]]+ ylab("\n") + theme(legend.position=c(0.15,0.8)),labels = c("G","H","I"),ncol=3,label_size=16)





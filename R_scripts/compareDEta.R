## Author: Alexander Cope
## A script for comparing \Delta\Etas across runs and producing plots.
## Note the input for these functions should be generated using the rescaleDEta.R script (selection_rescaled_by_mean.csv).

library(AnaCoDa)
library(deming)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggrepel)
library(ggnewscale)
library(ggpubr)


aa.groups<-list(c("N","Q","S","T","C","Z"), 
                 c("D","E","H","K","R"),
                 c("A","F","G","I","L","P","V","Y"))
groups<-c("Uncharged polar","Charged","Hydrophobic")

pc <- data.frame("AA" = unlist(aa.groups),"Property"=c(rep("Uncharged Polar",6),rep("Charged",5),rep("Hydrophobic",8)))


colors <- c("Charged" = "#F8766D","Hydrophobic" = "#00BA38", "Uncharged Polar"="#619CFF")


getSignificantCodons<-function(df.1,df.2)
{
  df.1[,"Significance"] <- rep("Not Significant",nrow(df.1)) 
  df.2[,"Significance"] <- rep("Not Significant",nrow(df.2)) 
  sig <- which((df.1[,"X0.025."] > df.2[,"X0.975."]) | (df.2[,"X0.025."] > df.1[,"X0.975."]))
  df.1[sig,"Significance"] <- "Significant"
  df.2[sig,"Significance"] <- "Significant"
  return(list(df.1,df.2))
}


## directory.1 directory containing data for x-axis
## directory.2 directory containing data for y-axis
## file.1 actual file name used for x-axis, will be combined with directory.1 (why did I separate these?)
## file.2 actual file name used for y-axis, will be combined with directory.2 (why did I separate these?)
## include.AA you can use this to only compare certain amino acids if you want
demingRegression <- function(directory.1,directory.2,file.1,file.2,include.AA=c())
{
 
  sel.1 <- read.table(paste0(directory.1,"Parameter_est/",file.1),sep=",",header=TRUE,stringsAsFactors=F)
  sel.2 <- read.table(paste0(directory.2,"Parameter_est/",file.2),sep=",",header=TRUE,stringsAsFactors=F) 

  # Code updates changed this column to Mean
  colnames(sel.1)[3] <- colnames(sel.2)[3] <- "Posterior"
 
  # Eliminate references if using non-scaled \Delta\eta
  sel.1 <- sel.1[which(sel.1$Posterior != 0),]
  sel.2 <- sel.2[which(sel.2$Posterior != 0),]
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
  b0 <- 0.0
  p.val <- pnorm(abs((b1-1)/sqrt(reg$variance[2,2])),lower.tail = F)*2
  std.err <- sqrt(reg$variance[2,2])


  # b1 <- 1
  # ci.b1 <- c(1,1)
  # b0 <- 0.0
  # p.val <- 1
  # std.err <- 0


  return(list(df=df,Slope=unname(b1),Intercept=unname(b0),Slope.CI=ci.b1,SD.1=sd.1,SD.2=sd.2,p.val=p.val,std.err = std.err))
}


## data data.frame output from demingRegression
## b1 slope from demingRegression
## b0 slope from demingRegression (usually set to 0)
## reg.ci regression confidence interval for slope. Will need to be slightly modified if want to incorporate intercept
## categories vector of length 2 giving X and Y axis labels
## ci Include posterior probability intervals for each point. Default is false
## file file name to output plot to
## title title of plot
## range.xy allows you to specify the range of the plot

plotDeta <-function(data,b1,b0,reg.ci=NULL,categories=c("X","Y"),ci = F,file="test.pdf",title="Regression",range.xy = NULL)
{
  data <- merge(data,pc,by="AA")
  data[which(data$AA == "S"),"AA"] <- "S[4]"
  data[which(data$AA == "Z"),"AA"] <- "S[2]"
  data[,"Codon"] <- unlist(lapply(data$Codon,deparse))
 

  conf.int.1 <- unlist(unname(reg.ci[1]))
  conf.int.2 <- unlist(unname(reg.ci[2]))

  l <- data.frame(s=c(b1,conf.int.2,conf.int.1, 1.0),ic=c(b0,b0,b0,0.0),Line=c("Slope","97.5% CI","2.5% CI","1:1 Line"),stringsAsFactors = F)
  l$Line <- factor(l$Line,levels=c("Slope","97.5% CI","2.5% CI","1:1 Line"))
  levels(l$Line) <- c(paste0("Slope: ",round(b1,4)),paste0("95% CI: ",round(conf.int.1,4),"-",round(conf.int.2,4)),paste0("95% CI: ",round(conf.int.1,4),"-",round(conf.int.2,4)),"1:1 Line")
  legend.colors <- c("black","black","red")
  lines.reg <- c("solid","dashed","solid")

  uniqueInitials <- paste(data$AA,data$Codon,sep=":")
  
  sig <- which(data$Significance == "Significant")
  if (length(sig) > 0)
  {
    uniqueInitials[-sig] <- NA
    data[-sig,"AA"] <- '"Not"~"Significant"'
  } else{
    uniqueInitials[1:length(uniqueInitials)] <- NA
    data[,"AA"] <- '"Not"~"Significant"'
  }
  p <- ggplot(data,aes(x=Posterior,y=Posterior.2))
  p <-(p + geom_abline(data=l,mapping=aes(slope=s,intercept=ic,linetype=Line,color=Line),size=0.15)
       + scale_colour_manual(values=legend.colors,name="",guide = guide_legend(order = 1),drop=F)
       + scale_linetype_manual(values=lines.reg,name="",guide = guide_legend(order = 1)))
  
  ## Include posterior probability intervals for each point
  if (ci)
  {
    if (is.null(range.xy))
    {
      range.xy <- range(c(data[,c("X0.025.","X0.975.")],data[,c("X0.025.2","X0.975.2")]),na.rm=T)
      #max.value <- max(range.xy)
      #range.xy <- c(-max.value,max.value)
    }
    if (range.xy[1] > -0.15)
    {
      range.xy[1] <- -0.2
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
  p <- p + geom_hline(yintercept=0,color="red",linetype="dashed",size=0.15) + geom_vline(xintercept=0,color="red",linetype="dashed",size=0.15)
 
  if (length(sig) > 0)
  {
    p <-(p
      + geom_point(data=data[-sig,],fill="grey",color="black",size=1,shape=21,alpha=0.5)
      + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),size=0.25,color="black",alpha=0.25) 
      + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),size=0.25,color="black",alpha=0.25))
    p <- (p + new_scale_color()
        + geom_point(data=data[sig,],mapping=aes(color=Property),size=3)
        + scale_colour_manual(values=colors,name="Significant")
        + geom_errorbar(data=data[sig,],mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),color="black") 
        + geom_errorbarh(data=data[sig,],mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),color="black")
        + geom_text_repel(label=uniqueInitials,size=4,max.iter=10000,force=5,box.padding=1.0,parse=T,segment.alpha=0.5,fontface="bold",max.overlaps=100))
    p <- p + guides(color = guide_legend(order=2))
 } else {
    p <-(p
        + geom_point(data=data,fill="grey",color="black",size=1,shape=21,alpha=0.5)
        + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),size=0.25,color="black",alpha=0.25) 
        + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),size=0.25,color="black",alpha=0.25))
 }
 p <- (p + labs(x=bquote(atop("   Favored"%<-%"Reference"%->%"Disfavored",.(categories[1])~Delta*eta)),y=bquote(atop("   Favored"%<-%"Reference"%->%"Disfavored",.(categories[2])~Delta*eta)),parse=T)
      + ggtitle(label=bquote(.(title)))
      )
  # p <- (p + labs(x=bquote(atop("Reference"%->%"Disfavored",.(categories[1])~Delta*eta)),y=bquote(atop("Reference"%->%"Disfavored",.(categories[2])~Delta*eta)),parse=T)
  #     + ggtitle(label=bquote(.(title)))
  #     )
  
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  p <- (p + theme_bw()
        + theme(axis.title=element_text(face="bold",size=12),axis.text=element_text(face="bold",size=8))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.32,0.82),legend.text=element_text(face="bold",size=10),legend.title=element_blank(),legend.key.width=unit(0.2,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
        + theme(plot.title = element_text(hjust = 0.5,size=12,face="bold")))
  ggsave(filename = file,plot=p,device="pdf",width=6,height=6)
  return(p)
 
}


compareDiffToMissenseError <- function(directory.1,directory.2,file.1,file.2,categories=c("X","Y"),absolute=TRUE,missense.error = "../scer_error_rates.tsv")
{
  sel.1 <- read.table(paste0(directory.1,"Parameter_est/",file.1),sep=",",header=TRUE,stringsAsFactors=F)
  sel.2 <- read.table(paste0(directory.2,"Parameter_est/",file.2),sep=",",header=TRUE,stringsAsFactors=F) 
  err <- read.table(missense.error,sep="\t",header=T,stringsAsFactors=F)
  # Code updates changed this column to Mean
  colnames(sel.1)[3] <- colnames(sel.2)[3] <- "Posterior"
  dfs <- getSignificantCodons(sel.1,sel.2)
  sel.1 <- dfs[[1]]
  sel.2 <- dfs[[2]]
  df <- plyr::join(sel.1,sel.2,by=c("AA","Codon","Significance"))
  colnames(df)[8] <- "Posterior.2"
  colnames(df)[9] <- "Std.Dev.2"
  colnames(df)[10] <- "X0.025.2"
  colnames(df)[11] <- "X0.975.2"
  rownames(sel.1) <- sel.1[,2]
  rownames(sel.2) <- sel.2[,2]

  df[,"Diff"] <- df$Posterior.2 - df$Posterior
  df[,"Disfavored"] <- "No difference"
  df[which(df$Significance == "Significant" & df$Diff > 0),"Disfavored"] <- categories[2]
  df[which(df$Significance == "Significant" & df$Diff < 0),"Disfavored"] <- categories[1]
  



  df <- merge(df,err,by=c("AA","Codon"))
  df <- merge(df,pc,by="AA")
  
  df <- df[which(df$Posterior != 0),]
  if (absolute)
  {
    df$Diff <- abs(df$Diff)
    y.label <- bquote("|"*Delta*eta[.(categories[2])]~"-"~Delta*eta[.(categories[1])]*"|")
  } else{
    y.label <- bquote(Delta*eta[.(categories[2])]~"-"~Delta*eta[.(categories[1])])
  }
  df[which(df$Significance != "Significant"),"Codon"] <- NA
  
  
  if ("e_m" %in% colnames(df) && "R_c" %in% colnames(df))
  {
    print(ppcor::pcor.test(df$Diff,df$e_m,df$R_c,method="spearman"))
    df <- df %>% pivot_longer(c("e_m","R_c"),names_to="Type",values_to="Rate")
    df$Type <- ifelse(df$Type=="e_m","Missense Error","Elongation")
    #
    p <- ggplot(df,aes(x=Rate,y=Diff,label=Codon)) + geom_point(aes(color=Disfavored)) +
         scale_colour_manual(values=c("Helix"="blue","Sheet"="green","Coil"="red","IDR"="gold","Structured"="cyan","No difference"="grey")) +
         labs(x="Expected Rate",y=y.label) +
         geom_text_repel(size=4,max.iter=10000,parse=T,segment.alpha=0.5,fontface="bold") +
         facet_wrap(~Type,scales="free_x") + 
         stat_cor(method="spearman",cor.coef.name="rho",label.sep="\n")  

    p <- (p + theme_bw() 
          + theme(axis.title=element_text(face="bold",size=12),axis.text=element_text(face="bold",size=8))
          + theme(axis.line = element_line(colour = "black"))
          + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
          + theme(legend.text=element_text(face="bold",size=8),legend.title=element_text(face="bold",size=8),legend.key.width=unit(0.2,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
          + theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"))) 
    } else if ("e_m" %in% colnames(df) && (!"R_c" %in% colnames(df))){
      p <- ggplot(df,aes(x=e_m,y=Diff,label=Codon)) + geom_point(aes(color=Disfavored)) +
         scale_colour_manual(values=c("Helix"="blue","Sheet"="green","Coil"="red","IDR"="gold","Structured"="cyan","No difference"="grey")) +
         labs(x="Expected Missense Rate",y=y.label) +
         geom_text_repel(size=4,max.iter=10000,parse=T,segment.alpha=0.5,fontface="bold") +
         stat_cor(method="spearman",cor.coef.name="rho",label.sep="\n")  

    p <- (p + theme_bw() 
          + theme(axis.title=element_text(face="bold",size=12),axis.text=element_text(face="bold",size=8))
          + theme(axis.line = element_line(colour = "black"))
          + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
          + theme(legend.text=element_text(face="bold",size=8),legend.title=element_text(face="bold",size=8),legend.key.width=unit(0.2,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
          + theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"))) 
    } else if ((!"e_m" %in% colnames(df)) && "R_c" %in% colnames(df)){
      p <- ggplot(df,aes(x=R_c,y=Diff,label=Codon)) + geom_point(aes(color=Disfavored)) +
         scale_colour_manual(values=c("Helix"="blue","Sheet"="green","Coil"="red","IDR"="gold","Structured"="cyan","No difference"="grey")) +
         labs(x="Expected Elongation Rate",y=y.label) +
         geom_text_repel(size=4,max.iter=10000,parse=T,segment.alpha=0.5,fontface="bold") +
         stat_cor(method="spearman",cor.coef.name="rho",label.sep="\n")  

    p <- (p + theme_bw() 
          + theme(axis.title=element_text(face="bold",size=12),axis.text=element_text(face="bold",size=8))
          + theme(axis.line = element_line(colour = "black"))
          + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
          + theme(legend.text=element_text(face="bold",size=8),legend.title=element_text(face="bold",size=8),legend.key.width=unit(0.2,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
          + theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"))) 
    }
    return(p)
}



plotSelectionParameters <- function(head.directory,target.directory)
{

  plots <- vector(mode="list",length=10)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Helix","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[1]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/ecoli_coil_vs_helix.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[2]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/ecoli_coil_vs_sheet.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[3]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Helix vs. Sheet")),file=paste0(target.directory,"/ecoli_helix_vs_sheet.pdf"),range.xy=NULL)
 
  output <- demingRegression(file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[4]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Structured Region","IDRs"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Structured vs. IDR")),file=paste0(target.directory,"/ecoli_ordered_vs_disorderd.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[5]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil (Structured only)","IDRs"),ci=T,title =bquote(atop("Comparison of Selection "*Delta*eta,"Coil vs. IDR")),file=paste0(target.directory,"/ecoli_coil_vs_idr.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[6]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix (Structured only)","IDRs"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Helix vs. IDR")),file=paste0(target.directory,"/ecoli_helix_vs_idr.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[7]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Sheet (Structured only)","IDRs"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Sheet vs. IDR")),file=paste0(target.directory,"/ecoli_sheet_vs_idr.pdf"),range.xy=NULL)
    
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[8]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/ecoli_coil_vs_helix_idr_removed.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[9]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/ecoli_coil_vs_sheet_idr_removed.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
  plots[[10]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/ecoli_helix_vs_sheet_idr_removed.pdf"),range.xy=NULL)
 

  return(plots)
}


plotSelectionParametersEmp <- function(head.directory,target.directory)
{

  plots <- vector(mode="list",length=10)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Helix","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[1]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/ecoli_coil_vs_helix_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[2]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/ecoli_coil_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[3]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Helix vs. Sheet")),file=paste0(target.directory,"/ecoli_helix_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","restart_5/"),file.path(head.directory,"Secondary_structures","Turn","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[4]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Turn"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Turn vs. Helix")),file=paste0(target.directory,"/ecoli_coil_vs_helix_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),file.path(head.directory,"Secondary_structures","Turn","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[5]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Sheet","Turn"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Turn vs. Sheet")),file=paste0(target.directory,"/ecoli_coil_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Turn","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[6]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Turn"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Turn vs. Coil")),file=paste0(target.directory,"/ecoli_helix_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
 
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Turn_Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Helix","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[7]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/ecoli_coil_vs_helix_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Turn_Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[8]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/ecoli_coil_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[9]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Helix vs. Sheet")),file=paste0(target.directory,"/ecoli_helix_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
  


  return(plots)
}


# head.directory <- "../Ecoli/Predicted/Results/"
# target.directory <- "../Images/"
# plots_codon_ecoli <- plotSelectionParameters(head.directory=head.directory,target.directory=target.directory)

# comb <- plot_grid(plots_codon_ecoli[[1]]+ggtitle("Coil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[2]]+ggtitle("Coil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[3]]+ggtitle("Helix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),plots_codon_ecoli[[4]]+ggtitle("Structured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/ecoli_relative_to_genome_optimal.pdf",width=12,height=14)

# comb <- plot_grid(plots_codon_ecoli[[5]]+ggtitle("Coil vs. IDR") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[6]]+ggtitle("Helix vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[7]]+ggtitle("Sheet vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C",""))
# ggsave2("../Images/ecoli_relative_to_genome_optimal_ss_vs_idr.pdf",width=12,height=14)

# comb <- plot_grid(plots_codon_ecoli[[8]]+ggtitle("Coil vs. Helix (IDRs Removed)") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[9]]+ggtitle("Coil vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[10]]+ggtitle("Helix vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C",""))
# ggsave2("../Images/ecoli_relative_to_genome_optimal_ss_idr_removed.pdf",width=12,height=14)

# helix.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Helix"),missense.error ="../ecoli_missense_rates_2021_2_26.tsv")
# sheet.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Sheet"),missense.error ="../ecoli_missense_rates_2021_2_26.tsv")
# sheet.helix <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Helix","Sheet"),missense.error ="../ecoli_missense_rates_2021_2_26.tsv")
# ord.dis <- compareDiffToMissenseError(file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("IDR","Structured"),missense.error ="../ecoli_missense_rates_2021_2_26.tsv")
  
# comb <- plot_grid(helix.coil+ggtitle("Shifts in Selection:\nCoil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.coil+ggtitle("Shifts in Selection:\nCoil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.helix+ggtitle("Shifts in Selection:\nHelix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   ord.dis+ggtitle("Shifts in Selection:\nStructured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/ecoli_missense_elongation_abs_model_estimates_2021_2_26.pdf",comb,width=12,height=14)


# helix.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Helix"),missense.error ="../ecoli_ms_missense_errors.tsv")
# sheet.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Sheet"),missense.error ="../ecoli_ms_missense_errors.tsv")
# sheet.helix <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Helix","Sheet"),missense.error ="../ecoli_ms_missense_errors.tsv")
# ord.dis <- compareDiffToMissenseError(file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("IDR","Structured"),missense.error ="../ecoli_ms_missense_errors.tsv")
  
# comb <- plot_grid(helix.coil+ggtitle("Shifts in Selection:\nCoil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.coil+ggtitle("Shifts in Selection:\nCoil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.helix+ggtitle("Shifts in Selection:\nHelix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   ord.dis+ggtitle("Shifts in Selection:\nStructured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/ecoli_missense_abs_empirical_estiamtes.pdf",comb,width=12,height=14)




# head.directory <- "../Scer/Statistical_power/Results/"
# target.directory <- "../Images/"
# output.100 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"100","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# output.50 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"50","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# output.25 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"25","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# output.10 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"10","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# output.5 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"5","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# output.1 <- demingRegression(file.path(head.directory,"RegionA","restart_5/"),file.path(head.directory,"1","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")

# df <- data.frame(slope=c(output.100$Slope,
#                          output.50$Slope,
                         
#                          output.10$Slope,
#                          output.1$Slope,
#                          1,
#                          output.100$Slope.CI[1],
#                          output.100$Slope.CI[2],
#                          output.50$Slope.CI[1],
#                          output.50$Slope.CI[2],
                         
#                          output.10$Slope.CI[1],
#                          output.10$Slope.CI[2],
#                          output.1$Slope.CI[1],
#                          output.1$Slope.CI[2]),ic=c(rep(0,13)),
#                  Group=c("100%","50%","10%","1%","0% (y=x)",rep("95% CI",8)),stringsAsFactors = F)

# df$Group <- factor(df$Group,levels=c("100%","50%","10%","1%","0% (y=x)","95% CI"))
# levels(df$Group) <- c("100%","50%","10%","1%","0% (y=x)","95% CI")
# lines.reg <- c(rep("solid",5),"dashed")
# legend.colors <- c("blue","purple","green","red","black","grey")
# p <- ggplot(df)
# p <-(p + geom_abline(data=df,mapping=aes(slope=slope,intercept=ic,colour=Group,linetype=Group))
#      + labs(x=bquote("Region A ("*Delta*eta*")"),y=bquote("Region B ("*Delta*eta*")")) 
     
#      + scale_color_manual(values=legend.colors)
#      + scale_linetype_manual(values=lines.reg)
#      + ggtitle(label="Effects of codons under different\nselective pressures"))
# p <- p + scale_x_continuous(limits = c(-1,1)) + scale_y_continuous(limits = c(-1,1))
# p <- p + guides(colour=guide_legend(title="%Region B under\nSeleciton for\nInefficiency"),
#                 linetype=guide_legend(title="%Region B under\nSeleciton for\nInefficiency"))
# p <- (p + theme_bw()
#       + theme(axis.title=element_text(size = 14,face="bold"),axis.text=element_text(size=14))
#       + theme(axis.line = element_line(colour = "black"))
#       + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#       + theme(legend.text=element_text(size=14),plot.title = element_text(hjust = 0.5,size=14)))

# ggsave2("../Images/simulated_expectation.pdf",p)




# head.directory <- "../Scer/Predicted/Results/"
# plots_codon_scer<- plotSelectionParameters(head.directory=head.directory,target.directory=target.directory)

# comb <- plot_grid(plots_codon_scer[[1]]+ggtitle("Coil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[2]]+ggtitle("Coil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[3]]+ggtitle("Helix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),plots_codon_scer[[4]]+ggtitle("Structured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/scer_relative_to_genome_optimal.pdf",width=12,height=14)

# comb <- plot_grid(plots_codon_scer[[5]]+ggtitle("Coil vs. IDR") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[6]]+ggtitle("Helix vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[7]]+ggtitle("Sheet vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C",""))
# ggsave2("../Images/scer_relative_to_genome_optimal_ss_vs_idr.pdf",width=12,height=14)

# comb <- plot_grid(plots_codon_scer[[8]]+ggtitle("Coil vs. Helix (IDRs Removed)") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[9]]+ggtitle("Coil vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[10]]+ggtitle("Helix vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C",""))
# ggsave2("../Images/scer_relative_to_genome_optimal_ss_idr_removed.pdf",width=12,height=14)


# helix.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Helix"))#,missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# sheet.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Sheet"))#,missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# sheet.helix <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Helix","Sheet"))#,missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# ord.dis <- compareDiffToMissenseError(file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("IDR","Structured"))#,missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
  

# comb <- plot_grid(helix.coil+ggtitle("Shifts in Selection:\nCoil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.coil+ggtitle("Shifts in Selection:\nCoil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.helix+ggtitle("Shifts in Selection:\nHelix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   ord.dis+ggtitle("Shifts in Selection:\nStructured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/scer_missense_elongation_abs_model_estimates.pdf",comb,width=12,height=14)

# helix.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Helix"),missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# sheet.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Sheet"),missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# sheet.helix <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Helix","Sheet"),missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
# ord.dis <- compareDiffToMissenseError(file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("IDR","Structured"),missense.error ="~/Research_projects/PANSE/Data/weingberg_ribo_dens.tsv")
  

# comb <- plot_grid(helix.coil+ggtitle("Shifts in Selection:\nCoil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.coil+ggtitle("Shifts in Selection:\nCoil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.helix+ggtitle("Shifts in Selection:\nHelix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   ord.dis+ggtitle("Shifts in Selection:\nStructured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/scer_elongation_abs_empirical_estimates.pdf",comb,width=12,height=14)

# helix.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Helix"),missense.error ="../scer_ms_missense_errors.tsv")
# sheet.coil <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Coil","Sheet"),missense.error ="../scer_ms_missense_errors.tsv")
# sheet.helix <- compareDiffToMissenseError(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("Helix","Sheet"),missense.error ="../scer_ms_missense_errors.tsv")
# ord.dis <- compareDiffToMissenseError(file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv",categories=c("IDR","Structured"),missense.error ="../scer_ms_missense_errors.tsv")
  
# comb <- plot_grid(helix.coil+ggtitle("Shifts in Selection:\nCoil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.coil+ggtitle("Shifts in Selection:\nCoil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   sheet.helix+ggtitle("Shifts in Selection:\nHelix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   ord.dis+ggtitle("Shifts in Selection:\nStructured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))
# ggsave2("../Images/scer_missense_abs_empirical_estiamtes.pdf",comb,width=12,height=14)


# ### Plots combining both species


# ## Secondary Structures
# title.1 <- ggdraw() + 
#   draw_label(
#     "   S. cerevisiae",
#     fontface = 'bold',
#     x = 0.5,
#     hjust = 0.5,  angle = 90, size=24
#   ) +
#   theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 2)
#   )

# title.2 <- ggdraw() + 
#   draw_label(
#     "   E. coli",
#     fontface = 'bold',
#     x = 0.5,
#     hjust = 0.5,  angle = 90, size=24
#   ) +
#   theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 2)
  
# )
# sc.comb <- plot_grid(plots_codon_scer[[1]]+ggtitle("Coil vs. Helix") +theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[2]]+ggtitle("Coil vs. Sheet")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[3]]+ggtitle("Helix vs. Sheet")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")), nrow=1,ncol=3,byrow=T,labels=c("A","B","C"))

# sc.comb <- plot_grid(title.1,sc.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))


# ec.comb <-plot_grid(plots_codon_ecoli[[1]]+theme(plot.title = element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[2]]+theme(plot.title =  element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[3]]+theme(plot.title =  element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=1,ncol=3,byrow=T,labels=c("D","E","F"))

# ec.comb <- plot_grid(title.2,ec.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))
# comb <- plot_grid(sc.comb,ec.comb,nrow=2,ncol=1)

# ggsave2("../Images/scer_ecoli_ss.pdf",comb,width=20,height=14)


# struct.idr <- plot_grid(plots_codon_scer[[4]]+ggtitle("S. cerevisiae\nStructured vs. IDRs") +theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#                   plots_codon_ecoli[[4]]+ggtitle("E. coli\nStructured vs. IDRs") +theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),ncol=2,nrow=1,labels=c("A","B"))

# sc.comb <- plot_grid(plots_codon_scer[[5]]+ggtitle("Coil vs. IDR") +theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[6]]+ggtitle("Helix vs. IDR")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[7]]+ggtitle("Sheet vs. IDR")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=1,ncol=3,labels=c("C","D","E"))
# sc.comb <- plot_grid(title.1,sc.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))


# ec.comb <- plot_grid(plots_codon_ecoli[[5]] +theme(plot.title = element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[6]]+theme(plot.title = element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[7]]+theme(plot.title = element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=1,ncol=3,labels=c("F","G","H"))

# ec.comb <- plot_grid(title.2,ec.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))
# comb <- plot_grid(struct.idr,sc.comb,ec.comb,nrow=3,ncol=1,align = "v")

# ggsave2("../Images/scer_ecoli_idr.pdf",comb,width=20,height=21)


# sc.comb <- plot_grid(plots_codon_scer[[8]]+ggtitle("Coil vs. Helix\nIDRs Removed") +theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[9]]+ggtitle("Coil vs. Sheet\nIDRs Removed")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_scer[[10]]+ggtitle("Helix vs. Sheet\nIDRs Removed")+theme(plot.title = element_text(size=18),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")), nrow=1,ncol=3,byrow=T,labels=c("A","B","C"))

# sc.comb <- plot_grid(title.1,sc.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))


# ec.comb <-plot_grid(plots_codon_ecoli[[8]]+theme(plot.title = element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[9]]+theme(plot.title =  element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
#   plots_codon_ecoli[[10]]+theme(plot.title =  element_blank(),plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=1,ncol=3,byrow=T,labels=c("D","E","F"))

# ec.comb <- plot_grid(title.2,ec.comb,nrow=1,ncol=2,rel_widths=c(0.05,1))
# comb <- plot_grid(sc.comb,ec.comb,nrow=2,ncol=1)

# ggsave2("../Images/scer_ecoli_ss_no_idrs.pdf",comb,width=20,height=14)


# head.directory <- "../Ecoli/Empirical/Results/"
# output <- demingRegression(file.path(head.directory,"Secondary_structures_begin_end_length_at_least_6_2_codon_for_termini","Core_helix","restart_5/"),file.path(head.directory,"Secondary_structures_begin_end_length_at_least_6_2_codon_for_termini","Start_helix_End_helix","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Core","Termini"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Helix (Core vs. Termini)")),file=paste0(target.directory,"/ecoli_helices_core_vs_termini.pdf"),range.xy=NULL)
  

# head.directory <- "../Scer/Empirical/Results/"
# output <- demingRegression(file.path(head.directory,"Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini","Core_helix","restart_5/"),file.path(head.directory,"Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini","Start_helix_End_helix","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
# p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Core","Termini"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Helix (Core vs. Termini)")),file=paste0(target.directory,"/scer_helices_core_vs_termini.pdf"),range.xy=NULL)
  

head.directory <- "../Scer/Predicted/Results/"
target.directory <- "../Images/mrna/"

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p<- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Paired","Unpaired"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil (mRNA Secondary Structure")),file=paste0(target.directory,"/scer_coil_paired_vs_unpaired.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Paired","Unpaired"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Helix (mRNA Secondary Structure)")),file=paste0(target.directory,"/scer_helix_paired_vs_unpaired.pdf"),range.xy=NULL)
 
output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Paired","Unpaired"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Sheet (mRNA Secondary Structure)")),file=paste0(target.directory,"/scer_sheet_pairved_vs_unpaired.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_paired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_coil_vs_helix_paired.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_paired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/scer_coil_vs_sheet_paired.pdf"),range.xy=NULL)
 
output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_paired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_paired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_helix_vs_sheet_paired.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_unpaired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_coil_vs_helix_unpaired.pdf"),range.xy=NULL)

output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Coil_unpaired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Sheet")),file=paste0(target.directory,"/scer_coil_vs_sheet_unpaired.pdf"),range.xy=NULL)
 
output <- demingRegression(file.path(head.directory,"Secondary_structure_paired_est_dM","Helix_unpaired","restart_5/"),file.path(head.directory,"Secondary_structure_paired_est_dM","Sheet_unpaired","restart_5/"),"selection_rescaled_to_genome_optimal.csv","selection_rescaled_to_genome_optimal.csv")
p <- plotDeta(output$df,b1=unname(output$Slope),b0 = unname(output$Intercept),reg.ci=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_helix_vs_sheet_unpaired.pdf"),range.xy=NULL)

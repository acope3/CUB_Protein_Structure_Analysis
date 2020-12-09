## Author: Alexander Cope
## A script for comparing \Delta\Etas across runs and producing plots.
## Note the input for these functions should be generated using the rescaleDEta.R script (selection_rescaled_by_mean.csv).

library(AnaCoDa)
library(deming)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggnewscale)


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
  return(list(df=df,Slope=unname(b1),Slope.CI=ci.b1,SD.1=sd.1,SD.2=sd.2,p.val=p.val,std.err = sqrt(reg$variance[2,2])))
}






plotDeta <-function(data,b1,b0,bounds=NULL,p.val,categories=c("X","Y"),ci = F,file="title",title="Regression",range.xy = NULL)
{
  data <- merge(data,pc,by="AA")
  data[which(data$AA == "S"),"AA"] <- "S[4]"
  data[which(data$AA == "Z"),"AA"] <- "S[2]"
  data[,"Codon"] <- unlist(lapply(data$Codon,deparse))
 

  
  conf.int.1 <- bounds[1]
  conf.int.2 <- bounds[2]
  l <- data.frame(s=c(b1,conf.int.2,conf.int.1, 1.0),ic=c(b0,0.0,0.0,0.0),Line=c("Slope","97.5% CI","2.5% CI","1:1 Line"),stringsAsFactors = F)
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
  
  if (ci)
  {
    if (is.null(range.xy))
    {
      range.xy <- range(c(data[,5:6],data[,10:11]),na.rm=T)
      max.value <- max(range.xy)
      range.xy <- c(-max.value,max.value)
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
  p <- p + geom_hline(yintercept=0,color="red",linetype="dashed",size=0.15) + geom_vline(xintercept=0,color="red",linetype="dashed",size=0.15)
  p<-(p
      + geom_point(data=data[-sig,],fill="grey",color="black",size=1,shape=21,alpha=0.5)
      + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),size=0.25,color="black",alpha=0.25) 
      + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),size=0.25,color="black",alpha=0.25)
      + new_scale_color()
      + geom_point(data=data[sig,],mapping=aes(color=Property),size=3)
      + scale_colour_manual(values=colors,name="Significant")
      + geom_errorbar(data=data[sig,],mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),color="black") 
      + geom_errorbarh(data=data[sig,],mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),color="black")
      + geom_text_repel(label=uniqueInitials,size=3,max.iter=10000,force=2,box.padding=0.75,parse=T,segment.alpha=0.5,fontface="bold")
      + labs(x=bquote(atop("   Favored"%<-%"Average"%->%"Disfavored",.(categories[1])~Delta*eta)),y=bquote(atop("   Favored"%<-%"Average"%->%"Disfavored",.(categories[2])~Delta*eta)),parse=T)
      + ggtitle(label=bquote(.(title)))
      )
  
  p <- p + guides(color = guide_legend(order=2))
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  p <- (p + theme_bw()
        + theme(axis.title=element_text(face="bold",size=12),axis.text=element_text(face="bold",size=8))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.2,0.85),legend.text=element_text(face="bold",size=8),legend.title=element_text(face="bold",size=8),legend.key.width=unit(0.2,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
        + theme(plot.title = element_text(hjust = 0.5,size=12,face="bold")))
  ggsave(filename = file,plot=p,device="pdf",width=6,height=6)
  return(p)
 
}





plotSelectionParamters <- function(head.directory,target.directory)
{

  plots <- vector(mode="list",length=10)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Helix","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[1]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_coil_vs_helix_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Coil","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[2]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta*,"Coil vs. Sheet")),file=paste0(target.directory,"/scer_coil_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structures","Helix","restart_5/"),file.path(head.directory,"Secondary_structures","Sheet","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[3]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_helix_vs_sheet_sig_codon_by_mean.pdf"),range.xy=NULL)
 
  output <- demingRegression(file.path(head.directory,"Ordered_disordered","Ordered","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[4]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Structured Region","Disordered Regions"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Structured vs. IDR")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_sig_codon.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[5]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Coil (Structured only)","IDRs"),ci=T,title =bquote(atop("Comparison of Selection "*Delta*eta,"Coil vs. IDR")),file=paste0(target.directory,"/scer_coil_vs_idr_sig_codon.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[6]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Helix (Structured only)","IDRs"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Helix vs. IDR")),file=paste0(target.directory,"/scer_helix_vs_idr_sig_codon.pdf"),range.xy=NULL)

  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),file.path(head.directory,"Ordered_disordered","Disordered","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[7]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Sheet (Structured only)","IDRs"),ci=T,title = bquote(atop("Comparison of Selection "*Delta*eta,"Sheet vs. IDR")),file=paste0(target.directory,"/scer_sheet_vs_idr_sig_codon.pdf"),range.xy=NULL)
    
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[8]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Coil","Helix"),ci=T,title =expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_coil_vs_helix_idr_removed_sig_codon_by_mean.pdf"),range.xy=NULL)
  
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Coil_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[9]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Coil","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta*,"Coil vs. Sheet")),file=paste0(target.directory,"/scer_coil_vs_sheet_idr_removed_sig_codon_by_mean.pdf"),range.xy=NULL)
   
  output <- demingRegression(file.path(head.directory,"Secondary_structure_order","Helix_ord","restart_5/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","restart_5/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv")
  plots[[10]] <- plotDeta(output$df,b1=unname(output$Slope),b0 = 0,bounds=output$Slope.CI,categories = c("Helix","Sheet"),ci=T,title = expression(atop("Comparison of Selection "*Delta*eta,"Coil vs. Helix")),file=paste0(target.directory,"/scer_helix_vs_sheet_idr_removed_sig_codon_by_mean.pdf"),range.xy=NULL)
 

  return(plots)
}



head.directory <- "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/"
target.directory <- "../Images/"
plots_codon <- plotSelectionParamters(head.directory=head.directory,target.directory=target.directory)

## Figure 1
comb <- plot_grid(plots_codon[[1]]+ggtitle("Coil vs. Helix") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[2]]+ggtitle("Coil vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[3]]+ggtitle("Helix vs. Sheet")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[4]]+ggtitle("Structured vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C","D"))

ggsave2("../Images/test.pdf",width=12,height=12)


## Figure 2
comb <- plot_grid(plots_codon[[5]] +ggtitle("Coil vs. IDR") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[6]]+ggtitle("Helix vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[7]]+ggtitle("Sheet vs. IDR")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C"))

ggsave2("../Images/secondary_structure_vs_idrs.pdf",width=12,height=12)

## Supplemental Figure
comb <- plot_grid(plots_codon[[8]]+ggtitle("Coil vs. Helix (IDRs Removed)") +theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[9]]+ggtitle("Coil vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),
  plots_codon[[10]]+ggtitle("Helix vs. Sheet (IDRs Removed)")+theme(plot.margin=unit(c(0.5,0.5,0.5,1),"lines")),nrow=2,ncol=2,labels=c("A","B","C"))

ggsave2("../Images/secondary_structure_comparison_no_idrs.pdf",width=12,height=12)




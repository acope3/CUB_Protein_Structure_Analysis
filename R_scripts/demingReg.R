library(AnaCoDa)
library(deming)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(ggrepel)

rescaleAGCT <- function(sel.file,parameter.file)
{
  data <- read.table(sel.file,sep=",",header=T,stringsAsFactors=F)
  parameter <- loadParameterObject(parameter.file)
  trace <- parameter$getTraceObject()
  amino <- aminoAcids()
  aa <- aminoAcids()
    for (a in aa)
    {
     if (a == "X" || a=="M" || a == "W") next
     data.tmp <- data[which(data$AA == a),]
     if (a == "R")
     {
      trace.aga <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"AGA",1,T)
      trace.agg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"AGG",1,T)
      trace.cga <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGA",1,T)
      trace.cgg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGG",1,T)
      trace.aga <- trace.aga - trace.agg
      trace.cga <- trace.cga - trace.cgg
      data.tmp[which(data.tmp$Codon == "AGA"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(mean(trace.aga),sd(trace.aga),quantile(trace.aga,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "CGA"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(mean(trace.cga),sd(trace.cga),quantile(trace.cga,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "AGG"),] <- NA
      data.tmp[which(data.tmp$Codon == "CGG"),] <- NA
      data.tmp[which(data.tmp$Codon == "CGT"),] <- NA
     }
     if (a == "L")
     {
      trace.ctc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTC",1,T)
      trace.ctt <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTT",1,T)
      trace.cta <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTA",1,T)
      trace.ctg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTG",1,T)
      trace.ctc <- trace.ctc - trace.ctt
      trace.cta <- trace.cta - trace.ctg
      data.tmp[which(data.tmp$Codon == "CTC"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(mean(trace.ctc),sd(trace.ctc),quantile(trace.ctc,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "CTA"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(mean(trace.cta),sd(trace.cta),quantile(trace.cta,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "CTT"),] <- NA
      data.tmp[which(data.tmp$Codon == "CTG"),] <- NA
      data.tmp[which(data.tmp$Codon == "TTG"),] <- NA
     }
     if (a == "I")
     {
      data.tmp[which(data.tmp$Codon == "ATA"),] <-NA
      data.tmp[which(data.tmp$Codon == "ATT"), ] <- NA
     }
     if (nrow (data.tmp) == 4)
     {
      codon.a <- data.tmp[1,"Codon"]
      codon.g <- data.tmp[3,"Codon"]
      trace.a <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.a,1,T)
      trace.g <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.g,1,T)
      trace.a <- trace.a - trace.g 
      data.tmp[1,c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(mean(trace.a),sd(trace.a),quantile(trace.a,c(0.025,0.975)))
      data.tmp[3,] <- NA
      data.tmp[4,] <- NA
     }
     if (nrow (data.tmp) == 2)
     {
      data.tmp[2,] <- NA
     }
     data[which(data$AA == a),] <- data.tmp
   }
   data <- data[which(!is.na(data$Posterior)),]
   return(data)
}

rescaleByMean<- function(sel.file,parameter.file)
{
  data <- read.table(sel.file,sep=",",header=T,stringsAsFactors=F)
  parameter <- loadParameterObject(parameter.file)
  trace <- parameter$getTraceObject()
  aa <- aminoAcids()
    for (a in aa)
    {
     if (a == "X" || a=="M" || a == "W") next
     data.tmp <- data[which(data$AA == a),]
     codons.wo.ref <- AAToCodon(a,T)
     traces <- vector(mode="list",length=(length(codons.wo.ref)+1))
     for (i in 1:length(codons.wo.ref))
     {
       traces[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codons.wo.ref[i],1,T)
     }
     traces[[(length(codons.wo.ref)+1)]] <- rep(0,length(traces[[1]]))
     df.traces <- t(as.data.frame(traces))
     mean.trace <- colMeans(df.traces)
     mean.deta <- mean(mean.trace)
     sd.trace <- sd(mean.trace)
     data.tmp[which(data.tmp$Posterior!=0),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Posterior!=0),c("Posterior","X0.025.","X0.975.")] - mean.deta
     data.tmp[which(data.tmp$Posterior==0),c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(-mean.deta,sd.trace,quantile(0-mean.trace,c(0.025,0.975)))
     data[which(data$AA == a),] <- data.tmp
   }
   return(data)
}



rescaleDEtaRun <- function(head.directory="../Scer/Predicted/Results/",targets=c("Secondary_structures","Ordered_disordered"))
{
  for (f in targets)
  {
    structure.loc <- file.path(head.directory,f)
    structures <- list.dirs(structure.loc,recursive=F)
    for (struct in structures)
    {
      sel.file.loc <- file.path(struct,"final_run","Parameter_est/")
      sel.file <- list.files(sel.file.loc,pattern="*_Selection.csv",full.names=T)
      parameter.file <- file.path(struct,"final_run","R_objects","parameter.Rda")
      print(sel.file)
      sel.rescaled <- rescaleAGCT(sel.file,parameter.file)
      write.table(sel.rescaled,file.path(sel.file.loc,"selection_rescaled_nucleotide.csv"),sep=",",col.names=T,row.names=F,quote=F)

      sel.rescaled <- rescaleByMean(sel.file,parameter.file)
      write.table(sel.rescaled,file.path(sel.file.loc,"selection_rescaled_by_mean.csv"),sep=",",col.names=T,row.names=F,quote=F)
    } 
  }
}



normalize <- function(dEta)
{
  aa <- AnaCoDa::aminoAcids()
  for (a in aa)
  {
    if (a=="W" || a=="M"||a=="X") next
    row.ind <- which(dEta[,1] == a)
    rows <- dEta[row.ind,]
    mu <- mean(rows[,"Posterior"])
    rows[,c("Posterior","X0.025.","X0.975.")] <- rows[,c("Posterior","X0.025.","X0.975.")] - mu
    dEta[row.ind,] <- rows
  }
  return(dEta)
}



## Some older runs had a different structure for the saved parameter object. This function
## fixes this, if necessary
fixParameter <- function(parameter.file.to.fix)
{
  load(parameter.file.to.fix)
  
  aminoAcids()
  paramBase$grouplist <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z")
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection", "model",  
                "mutationPrior", "mutationTrace", "selectionTrace", 
                "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
                "observedSynthesisNoiseTrace", "withPhi"),
       file=parameter.file.to.fix)
}

getSignificantCodons<-function(df.1,df.2)
{
  df.1[,"Significance"] <- rep("Not Significant",nrow(df.1)) 
  df.2[,"Significance"] <- rep("Not Significant",nrow(df.2)) 
  sig <- which((df.1[,5] > df.2[,6]) | (df.2[,5] > df.1[,6]))
  df.1[sig,"Significance"] <- "Significant"
  df.2[sig,"Significance"] <- "Significant"
  return(list(df.1,df.2))
}



demingApproach <- function(directory.1,directory.2,file.1,file.2,include.AA=c(),optimal.as.ref = F,normalize.deta=T,include.orig.ref=T)
{
  names.aa <- aminoAcids()
  sd.1 <- c()
  sd.2 <- c()
  sel.1 <- read.table(paste0(directory.1,"Parameter_est/",file.1),sep=",",header=TRUE,stringsAsFactors=F)
  sel.2 <- read.table(paste0(directory.2,"Parameter_est/",file.2),sep=",",header=TRUE,stringsAsFactors=F)
 
 
  if (optimal.as.ref)
  {
    sel.1 <- optimalAsReference(sel.1)
    sel.2 <- optimalAsReference(sel.2)
  } 
  rownames(sel.1) <- sel.1[,2]
  rownames(sel.2) <- sel.2[,2]
  if (length(include.AA) != 0)
  {
    sel.1 <- sel.1[c(which(sel.1$AA %in% include.AA)),]
    sel.2 <- sel.2[c(which(sel.2$AA %in% include.AA)),]
  }
  ## At this point, this section of the code is probably obsolete by using the output rescaleByMean function for the deming regression
  if (include.orig.ref)
  {
    for(aa in names.aa)
    {
      removed <- 0
      if(aa == "M" || aa == "W" || aa == "X") next
      sd.1.tmp <- sel.1[which(sel.1$AA == aa),"Std.Dev"]^2
      sd.2.tmp <- sel.2[which(sel.2$AA == aa),"Std.Dev"]^2
      total.sd.1 <- sum(sd.1.tmp)
      total.sd.2 <- sum(sd.2.tmp)
      sd.1.tmp <- sqrt(rep(total.sd.1/length(sd.1.tmp),length(sd.1.tmp)))
      sd.2.tmp <- sqrt(rep(total.sd.2/length(sd.2.tmp),length(sd.2.tmp)))
      sd.1 <- c(sd.1,sd.1.tmp)
      sd.2 <- c(sd.2,sd.2.tmp)
    }
    if (normalize.deta)
    {
      sel.1.for.reg <- normalize(sel.1)
      sel.2.for.reg <- normalize(sel.2)
      sel.1 <- sel.1.for.reg
      sel.2 <- sel.2.for.reg
    }
  } ## End obsolete portion
  else{
    sel.1<- sel.1[which(sel.1$Posterior != 0),]
    sel.2<- sel.2[which(sel.2$Posterior != 0),]
    sd.1 <- sel.1$Std.Dev
    sd.2 <- sel.2$Std.Dev
  }
  dfs <- getSignificantCodons(sel.1,sel.2)
  sel.1 <- dfs[[1]]
  sel.2 <- dfs[[2]]
  df <- plyr::join(sel.1,sel.2,by=c("AA","Codon","Significance"))
  colnames(df)[8] <- "Posterior.2"
  colnames(df)[9] <- "Std.Dev.2"
  colnames(df)[10] <- "X0.025.2"
  colnames(df)[11] <- "X0.975.2"
  reg <- deming(Posterior.2 ~ Posterior+0,data=df,xstd=sd.1,ystd=sd.2)
  b1 <- reg$coefficients[2]
  #b0 <- reg$coefficients[1]
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



plotResultsAll<-function(data,b1,b0,categories=c("X","Y"),ci = F,bounds=NULL,file="title",title="Regression",range.xy = NULL)
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


plotResultsNoReg<-function(data,p.val,categories=c("X","Y"),ci = F,file="title",title="Regression",range.xy = NULL)
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
  uniqueInitials[sig] <- paste0(uniqueInitials[sig],'*"*"',sep="")
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
      + geom_point(aes(color=AA),size=4)
      #+ geom_point(size=4,alpha=0.6)
      + geom_text_repel(label=uniqueInitials,size=5,max.iter=5000,box.padding = 0.5,point.padding=0.4,nudge_x=0.01,nudge_y=0.01,force=10,parse=T,segment.alpha=0.7,fontface="bold")
      #+ labs(x=bquote("Selection("~Delta*eta*"):"~.(categories[1])),y=bquote("Selection("~Delta*eta*"):"~.(categories[2])),color="3rd Position") 
      + labs(x=bquote(atop(.(categories[1]),atop("Selection Coefficients: G or T","for A/C "%<-%" Increasing Selection "%->%" for G/T"))),y=bquote(atop(.(categories[2]),atop("Selection Coefficients: G or T","for A/C"%<-%" Increasing Selection "%->%"for G/T"))),color="",parse=T)
      #+ scale_colour_manual(values=c("Pur-Pur"="blue","Pyr-Pyr"="orange","Pur-Pyr"="green","Pyr-Pur"="yellow"),guide = guide_legend(order = 2),drop=F)
      + scale_colour_discrete(breaks = unique(data$AA),labels=lapply(unique(data$AA),function(x)parse(text=x)),guide=guide_legend(nrow=ifelse(length(unique(data$AA)) > 5,5,length(unique(data$AA))),ncol=ifelse(length(unique(data$AA)) > 5,3,1)))
      #+ ggtitle(label=title)
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
        + theme(axis.title=element_text(face="bold",size=14),axis.text=element_text(face="bold",size=12))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.2,0.8),legend.text=element_text(face="bold",size=14),legend.title=element_text(face="bold",size=14),legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
        + theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")))
  ggsave(filename = file,plot=p,device="pdf",width = 7,height = 7)
  return(p)
 
}

optimalAsReference <- function(param.2)
{
  updated.param.2 <- data.frame()
  rownames(param.2) <- param.2[,"Codon"]
  aa <- unique(param.2[,"AA"])
  for (a in aa)
  {
    codons <- AAToCodon(a)
    ## Create temporary data frames for modifying values
    tmp.2 <- param.2[codons,] ## "Selection" parameter
    current.reference.row <- which(tmp.2[,"Posterior"]==0)
    optimal.parameter.value <- min(tmp.2[,"Posterior"])
    ## No reason to do anything if optimal value is 0
    if (optimal.parameter.value != 0.0)
    {
      tmp.2[,c("Posterior","X0.025.","X0.975.")] <- tmp.2[,c("Posterior","X0.025.","X0.975.")] - optimal.parameter.value
      ##Get row of the optimal codon, which should be 0
      optimal.codon.row <- which(tmp.2[,"Posterior"]==0.0)
      tmp.2[current.reference.row,"Std.Dev"] <- tmp.2[optimal.codon.row,"Std.Dev"]
      tmp.2[current.reference.row,c("X0.025.","X0.975.")] <- tmp.2[optimal.codon.row,c("X0.025.","X0.975.")] + tmp.2[current.reference.row,"Posterior"]
      ## Can now change optimal codon values to 0.0
      tmp.2[optimal.codon.row,c("Posterior","Std.Dev","X0.025.","X0.975.")] <- 0.0
      ## Find corresponding reference value for other parameter
   }

    updated.param.2 <- rbind(updated.param.2,tmp.2)
  }

  updated.param.2 <- updated.param.2[,c("AA", "Codon", "Posterior","Std.Dev","X0.025.", "X0.975.")]
  return(updated.param.2)
}


plotDemingRegressionSecondaryStructuresAndOrdered <- function(head.directory)
{
  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p5 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Coil","Disordered Coil"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil (Structured vs. Disordered)",file="../Images_for_dissertation/scer_coil_ordered_vs_disordered.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p6 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Helix","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nHelix (Structured vs. Disordered)",file="../Images_for_dissertation/scer_helix_ordered_vs_disordered.pdf")


  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p7 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil vs. Helix (Structured Regions)",file="../Images_for_dissertation/scer_ordered_coil_vs_ordered_helix.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p8 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil vs. Sheet (Structured Regions)",file="../Images_for_dissertation/scer_ordered_coil_vs_ordered_sheet.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p9 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nHelix vs. Sheet (Structured Regions)",file="../Images_for_dissertation/scer_ordered_helix_vs_ordered_sheet.pdf")


  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p10 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Disordered Helix",file="../Images_for_dissertation/scer_disordered_helix_vs_disordered_coil_disordered.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p11 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Sheet",file="../Images_for_dissertation/scer_disordered_coil_vs_sheet.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p12 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Helix vs. Sheet",file="../Images_for_dissertation/scer_disordered helix_vs_sheet.pdf")

}


plotDemingRegressionOrdered <- function(head.directory)
{
  output <- demingApproach(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p4 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nStructured vs. Disordered Regions",file="../Images_for_dissertation/scer_ordered_vs_disorderd.pdf")

  ## Ordered and Disordered conserved vs non-conserved
  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p.cons <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured vs. Disordered (Conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_vs_disorderd_conserved_sacch.pdf"))

  output <- demingApproach(file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p.non <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Structured Region","Disordered Region"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured vs Disordered (Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_vs_disorderd_nonconserved_sacch.pdf"))


  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p.ord <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Residues","Non-conserved Residues"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nStructured Regions (Conserved vs Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_ordered_conserved_vs_nonconserved_sacch.pdf"))

  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p.dis <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Conserved Residues","Non-conserved Residues"),ci=T,bounds=output$Slope.CI,title = bquote("Comparison of Selection on CUB:\nDisordered Regions (Conserved vs Non-conserved Residues)"),file=paste0("../Images_for_dissertation/scer_disorderd_conserved_vs_nonconserved_sacch.pdf"))


}


plotDemingRegressionSecondaryStructures <- function(head.directory)
{
  output <- demingApproach(file.path(head.directory,"Secondary_structures","Turn_Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p1 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Turn and Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_dissertation/Empirical_structures/scer_coil_vs_helix_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures","Turn_Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p2 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Turn and Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_dissertation/Empirical_structures/scer_coil_vs_sheet_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p3 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_dissertation/Empirical_structures/scer_helix_vs_sheet_already_normalized.pdf")


  output <- demingApproach(file.path(head.directory,"Secondary_structures_conserved","Coil","final_run/"),file.path(head.directory,"Secondary_structures_conserved","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p1 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_dissertation/scer_cons_coil_vs_helix_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures_conserved","Coil","final_run/"),file.path(head.directory,"Secondary_structures_conserved","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p2 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_dissertation/scer_cons_coil_vs_sheet_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures_conserved","Helix","final_run/"),file.path(head.directory,"Secondary_structures_conserved","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p3 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_dissertation/scer_cons_helix_vs_sheet_already_normalized.pdf")


  output <- demingApproach(file.path(head.directory,"Secondary_structures_unconserved","Coil","final_run/"),file.path(head.directory,"Secondary_structures_unconserved","Helix","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p1 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="../Images_for_dissertation/scer_var_coil_vs_helix_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures_unconserved","Coil","final_run/"),file.path(head.directory,"Secondary_structures_unconserved","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p2 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="../Images_for_dissertation/scer_var_coil_vs_sheet_already_normalized.pdf")

  output <- demingApproach(file.path(head.directory,"Secondary_structures_unconserved","Helix","final_run/"),file.path(head.directory,"Secondary_structures_unconserved","Sheet","final_run/"),"selection_rescaled_by_mean.csv","selection_rescaled_by_mean.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F)
  p3 <- plotResultsAll(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="../Images_for_dissertation/scer_var_helix_vs_sheet_already_normalized.pdf")

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
    output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Coil","Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_coil_vs_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[1]][[i]] <- p
    

    output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[2]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[3]][[i]] <- p
    
    
    output <- demingApproach(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Region","Disordered Regions"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[4]][[i]] <- p
    
    ## Ordered Conserved vs Non-conserved
    output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Region (Conserved Residues)","Structured Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[5]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Disordered Region (Conserved Residues)","Disordered Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[6]][[i]] <- p

    ## Ordered and Disordered conserved vs non-conserved
    output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Region (Conserved Residues)","Disordered Region (Conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_conserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[7]][[i]] <- p

    output <- demingApproach(file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch_no_ncast","Disordered","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Region (Non-conserved Residues)","Disordered Region (Non-conserved Residues)"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_vs_disorderd_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[8]][[i]] <- p

    # Ordered vs Disorderd Secondary Structures (Coil and Helix Only)
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Coil","Disordered Coil"),ci=T,title =bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_disordered_coil_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[9]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Structured Helix","Disordered Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_helix_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[10]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Coil","Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_ordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[11]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_coil_vs_ordered_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[12]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Disordered Coil","Disordered Helix"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_coil_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[13]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Disordered Coil","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[14]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Sheet_ord","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_ordered_helix_vs_ordered_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[15]][[i]] <- p
    
    output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"selection_rescaled_nucleotide.csv","selection_rescaled_nucleotide.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]])
    p <- plotResultsNoReg(output$df,categories = c("Disordered Helix","Sheet"),ci=T,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0(target.directory,"/scer_disordered_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
    plots[[16]][[i]] <- p

    

  }
  return(plots)
}


# aa.groups <- list(c("C","D","E","F","H","I","K","N","Q","Y","Z"),
#                   c("A","G","P","S","T","V"),
#                   c("R","L"))



plots <- plotDEtaByAA(head.directory="../Scer/Predicted/Results/",target.directory="../Images_for_dissertation/Examples/")

##### Secondary Structures ###############################
subplots <- plots[[1]]
# legend <- get_legend(subplots[[1]])
coil.helix <- plot_grid(subplots[[1]]+xlab("Coil")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain"))+ggtitle("Uncharged Polar\nAmino Acids"), 
                      subplots[[2]] +xlab("Coil") + ylab("\n")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain"))+ ggtitle("Charged\nAmino Acids"), 
                      subplots[[3]] +xlab("Coil")+ ylab("\n") + theme(legend.position=c(0.15,0.8),axis.title.x = element_text(face="plain"))+ ggtitle("Non-polar\nAmino Acids"),labels = c("A","B","C"),ncol=3,label_size=16)
subplots <- plots[[2]]
coil.sheet <- plot_grid(subplots[[1]] +xlab("Coil")+theme(legend.position=c(0.2,0.75),axis.title.x = element_text(face="plain")), 
                        subplots[[2]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.1,0.8),axis.title.x = element_text(face="plain")), 
                        subplots[[3]] +xlab("Coil") + ylab("\n") +theme(legend.position=c(0.15,0.85),axis.title.x = element_text(face="plain")),labels = c("D","E","F"),ncol=3,label_size=16)
subplots <- plots[[3]]
helix.sheet <- plot_grid(subplots[[1]] +theme(legend.position=c(0.2,0.75)), 
                          subplots[[2]] + ylab("\n")+ theme(legend.position=c(0.1,0.8)), 
                          subplots[[3]]+ ylab("\n") + theme(legend.position=c(0.15,0.8)),labels = c("G","H","I"),ncol=3,label_size=16)


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


title <- ggdraw() + draw_label(
    "Comparison of Selection Coefficients: Secondary Structures ",
      fontface = 'bold',
      x = 0,
      hjust =0,size=25
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 2)
    )

all.plots <- plot_grid(
    title, coil.helix, coil.sheet,helix.sheet,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights=c(0.1,1,1,1)
  )
ggsave2(file=paste0("../Images_for_dissertation/Examples/SS_ordered_by_aa.pdf"),plot=all.plots,width=18,height=16)



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





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


rescaleAGCT <- function(data)
{
  aa <- aminoAcids()
  for (a in aa)
  {
    if (a == "X" || a=="M" || a == "W") next
    data.tmp <- data[which(data$AA == a),]
    if (a == "R")
    {
      data.tmp[which(data.tmp$Codon == "AGA"),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Codon == "AGA"),c("Posterior","X0.025.","X0.975.")] - data.tmp[which(data.tmp$Codon == "AGG"),c("Posterior")]
      data.tmp[which(data.tmp$Codon == "CGA"),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Codon == "CGA"),c("Posterior","X0.025.","X0.975.")] - data.tmp[which(data.tmp$Codon == "CGG"),c("Posterior")]
      data.tmp[which(data.tmp$Codon == "AGG"),] <- NA
      data.tmp[which(data.tmp$Codon == "CGG"),] <- NA
    }
    if (a == "L")
    {
      data.tmp[which(data.tmp$Codon == "CTC"),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Codon == "CTC"),c("Posterior","X0.025.","X0.975.")] - data.tmp[which(data.tmp$Codon == "CTT"),c("Posterior")]
      data.tmp[which(data.tmp$Codon == "CTA"),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Codon == "CTA"),c("Posterior","X0.025.","X0.975.")] - data.tmp[which(data.tmp$Codon == "CTG"),c("Posterior")]
      data.tmp[which(data.tmp$Codon == "CTC"),] <- NA
      data.tmp[which(data.tmp$Codon == "CTG"),] <- NA
    }
    if (a == "I")
    {
      data.tmp[which(data.tmp$Codon == "ATA"),] <-NA
    }
    if (nrow (data.tmp) == 4)
    {
      data.tmp[1,c("Posterior","X0.025.","X0.975.")] <- data.tmp[1,c("Posterior","X0.025.","X0.975.")] - data.tmp[3,c("Posterior")]
      data.tmp[3,] <- NA
    }
    data[which(data$AA == a),] <- data.tmp
  }
  data[,c("Posterior","X0.025.","X0.975.")] <- -data[,c("Posterior","X0.975.","X0.025.")]
  return(data)
}


demingApproach <- function(directory.1,directory.2,file.1,file.2,include.AA=c(),exclude.codon=c(),optimal.as.ref = F,normalize.deta=T,include.orig.ref=T,rescale.agct=F)
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
  } else if (rescale.agct){
    sel.1 <- rescaleAGCT(sel.1)
    sel.1 <- sel.1[which(!is.na(sel.1$Codon)),]
    sel.2 <- rescaleAGCT(sel.2)
    sel.2 <- sel.2[which(!is.na(sel.2$Codon)),]
  }
  rownames(sel.1) <- sel.1[,2]
  rownames(sel.2) <- sel.2[,2]
  if (length(include.AA) != 0)
  {
    sel.1 <- sel.1[c(which(sel.1$AA %in% include.AA)),]
    sel.2 <- sel.2[c(which(sel.2$AA %in% include.AA)),]
  }
  if (length(exclude.codon) != 0)
  {
    sel.1 <- sel.1[-c(which(sel.1$Codon %in% exclude.codon)),]
    sel.2 <- sel.2[-c(which(sel.2$Codon %in% exclude.codon)),]
  }
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
    }
  } else{
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
  b0 <- reg$coefficients[1]
  ci.b1 <- reg$ci[2,]
  p.val <- pnorm(abs((b1-1)/sqrt(reg$variance[2,2])),lower.tail = F)*2
  return(list(df=df,Slope=b1,Slope.CI=ci.b1,SD.1=sd.1,SD.2=sd.2,p.val=p.val,std.err = sqrt(reg$variance[2,2])))
}

plotResults<-function(data,b1,b0,categories=c("X","Y"),ci = F,bounds=NULL,file="title",title="Regression",range.xy = NULL)
{
  data[which(data$AA == "S"),"AA"] <- "S[4]"
  data[which(data$AA == "Z"),"AA"] <- "S[2]"
  if (!is.null(bounds))
  {
    conf.int.1 <- bounds[1]
    conf.int.2 <- bounds[2]
    l <- data.frame(s=c(b1,conf.int.2,conf.int.1, 1.0),ic=c(b0,0.0,0.0,0.0),Line=c("Deming Slope","97.5% CI","2.5% CI","1:1 Line"),stringsAsFactors = F)
    l$Line <- factor(l$Line,levels=c("Deming Slope","97.5% CI","2.5% CI","1:1 Line"))
    levels(l$Line) <- c(paste0("MLE Slope: ",round(b1,3)),paste0("95% CI: ",round(conf.int.1,3),"-",round(conf.int.2,3)),paste0("95% CI: ",round(conf.int.1,3),"-",round(conf.int.2,3)),"1:1 Line")
    legend.colors <- c("black","black","red")
    lines.reg <- c("solid","dashed","solid")
    
  } else{
    l <- data.frame(s=c(b1,1.0),ic=c(b0,0.0),Line=c("Model II Regression","y = x"))
    legend.colors <- c("black","black")
    lines.reg <- c("solid","dashed")
  }
  codons <- data$Codon
  nuc <- strsplit(codons,split="")
  ct <- unlist(lapply(nuc,function(x){
    if (x[3] == "C" || x[3] == "T") return("C or T")
    else return("A or G")
  }))
  data["Third.Nucleotide"] <- ct
  uniqueInitials <- data$AA
  sig <- which(data$Significance == "Significant")
  uniqueInitials[sig] <- paste0(uniqueInitials[sig],"^'*'")
  initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))
  p <- ggplot(data,aes(Posterior,Posterior.2,label=uniqueInitials))
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
    p <- (p + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2,width=0),alpha=0.2,color="black") 
          + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.,height=0),alpha=0.2,color="black"))
  } else{
    range.xy <- range(c(data[,3],data[,7]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  }
  p<-(p + new_scale_color()
       + geom_text(aes(color=Third.Nucleotide,fontface="bold"),size=12,alpha=0.6,parse=T)
       #+ labs(x=bquote("Selection("~Delta*eta*"):"~.(categories[1])),y=bquote("Selection("~Delta*eta*"):"~.(categories[2])),color="3rd Position") 
       + labs(x=bquote(atop(.(categories[1]),"Selection Coefficients: G or T")),y=bquote(atop(.(categories[2]),"Selection Coefficients: G or T")),color="3rd Position")
       + scale_color_manual(values=c("blue","orange"),guide = guide_legend(order = 2))
       + scale_shape_manual(values=uniqueInitials,guide = guide_legend(order = 2))
       + ggtitle(label=title))
  p <- p + guides(color = guide_legend(override.aes = list(label="",size=1)))
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.2,label=deparse(cor.exp),parse=T,fontface="bold",size=8)
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.10,label="*p < 0.05",parse=F,size=8)
  p <- (p + theme_bw()
       + theme(axis.title=element_text(face="bold",size=14),axis.text=element_text(face="bold",size=12))
       + theme(axis.line = element_line(colour = "black"))
       + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
       + theme(legend.position = c(0.2,0.8),legend.text=element_text(face="bold",size=14),legend.title=element_text(face="bold",size=14),legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"),legend.spacing.y=unit(0.1,"cm"))
       + theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")))


  ## Modified from https://stackoverflow.com/questions/23588127/match-legend-text-color-in-geom-text-to-symbol/37805504#37805504
  pGrob <- ggplotGrob(p)
  g.b   <- pGrob[["grobs"]][[which(pGrob$layout$name=="guide-box")]]
  index <- 0
  label_text <- c()
  ## One legend has two labels, the other has 3, so find the one with 2
  while (length(label_text) != 2)
  {
    index <- index + 1
    l <- g.b[[1]][[index]][["grobs"]]
    # get grobs for legend symbols (extract colour)
    lg <- l[sapply(l, function(i) grepl("GRID.text", i))]
    clr <- mapply(FUN=function(x){x$gp$col},x=lg)
    gb  <- which(grepl("guide-box", pGrob$layout$name))
    gb2 <- which(grepl("guides", pGrob$grobs[[gb]]$layout$name))
    label_text <- which(grepl("label",pGrob$grobs[[gb]]$grobs[[gb2[index]]]$layout$name))
  }
  pGrob$grobs[[gb]]$grobs[[gb2[index]]]$grobs[label_text] <- 
    mapply(FUN = function(x, y) {x[["children"]][[1]][["children"]][[1]]$gp <- gpar(col =y); return(x)},
           x =   pGrob$grobs[[gb]]$grobs[[gb2[index]]]$grobs[label_text],
           y =  clr, SIMPLIFY = FALSE)
  ggsave(filename = file,plot=pGrob,device="pdf",width = 7,height = 7)
  return(as_ggplot(pGrob))
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



upper.panel.plot <- function(df,slope,ci,p.val,...){
  abline(0, 1, col = "blue", lty = 2)
  
  #points(x,y, ...)
  aa.counts <- table(df$AA)
  df <- df[with(df,order(AA,Posterior)),]
  x <- df$Posterior
  y <- df$Posterior.2
  new.label <- paste0(df$AA,unlist(lapply(unique(df$AA),function(i){seq(aa.counts[i])})))
  text(df$Posterior,df$Posterior.2,labels = new.label,cex = 1.05,vfont=c("sans serif","bold"))
  y.up <- df[,"X0.975.2"]
  y.low <- df[,"X0.025.2"]
  epsilon <- range(x, na.rm = T) * 0.1
  segments(x, y.low, x, y.up, ...)


  x.up <- df[,"X0.975."]
  x.low <- df[,"X0.025."]
  epsilon <- range(y, na.rm = T) * 0.1
  segments(x.low, y, x.up, y, ...)

  
  lm.line <- lm(y~x+0, na.action = "na.exclude")
  
  
  R2 <- summary(lm.line)$r.squared
  rho <- ifelse(slope > 0, sqrt(R2), -sqrt(R2)) #make sure rho has correct sign

  xlim <- range(df[,c("X0.025.","X0.975.")], na.rm = T)
  ylim <- range(df[,c("X0.025.2","X0.975.2")], na.rm = T)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  slope <- round(slope, 3)
  if(!(ci[1]<1  && ci[2]>1)){
    abline(b=slope,a=0, lwd = 2)
    eq <- paste0("y = ",sprintf("%.3f", slope), "x *")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  }else{
    abline(b=slope,a=0, lwd = 2)
    eq <- paste0("y = ", sprintf("%.3f", slope), "x")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  }
  text(xlim[1] + width * 0.0005, ylim[2] - height * 0.3, paste0("CI:",round(ci[1],3),"-",round(ci[2],3)), pos = 4, cex = 1.5)
  text(xlim[1] + width * 0.0005, ylim[2] - height * 0.4, paste0("p = ",round(p.val,3)), pos = 4, cex = 1.5)
  if(slope > 0){
    text(xlim[2] - width * 0.04, ylim[1] + height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }else{
    text(xlim[2] - width * 0.04, ylim[2] - height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }
}

library(AnaCoDa)
library(deming)
library(ggplot2)
library(viridis)
library(ggnewscale)
library(grid)
library(cowplot)
library(ggpubr)
#library(patchwork)
# runs <- list.dirs("Scer/Exp_conservative_homology_min_200/Results/Secondary_structures",recursive=F)
# num.comp <- factorial((length(runs)))/(2*factorial(length(runs)-2))
# slope.ci <- data.frame(Category.1=character(length=num.comp),
#                        Category.2=character(length=num.comp),
#                        Slope=numeric(length=num.comp),
#                        Lower.Slope.CI.JK = numeric(length=num.comp),
#                        Upper.Slope.CI.JK = numeric(length=num.comp),
#                        stringsAsFactors=F)
# index <- 1
# runs <- rev(runs)
# numMixtures <- length(runs)
# mixture.name <- c()
# mixture.name.tmp <- strsplit(runs,"/",fixed = T)
# for (k in 1:length(mixture.name.tmp))
# {
#   mixture.name <- c(mixture.name,mixture.name.tmp[[k]][5])
# }
# pdf("ss_scer_exp_cons_structure_min_200_nterm_cterm_sep.pdf",width=25,height=20)
# mat <- matrix(rep(0,numMixtures*numMixtures),
#               nrow = numMixtures, ncol = numMixtures, byrow = TRUE)
# #count <- num.comp + numMixtures
# count <- 1
# for(i in 1:numMixtures){
#   for(j in 1:numMixtures){
#     if(i<=j){
#       mat[i,j] <-count
#       count <- count + 1
#     }
#   }
# }
# nf <- layout(mat,widths=c(rep(2,numMixtures)),heights=c(rep(2,numMixtures)),respect=FALSE)
# par(mar=c(1,1,1,1))
# 
# for (i in 1:length(runs))
# {
#   for (j in 1:length(runs))
#   {
#     if(i==j)
#     {
#       plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="",xaxt='n',yaxt='n',ann=FALSE)
#       if(is.null(mixture.name)){
#         text(x = 0.5, y = 0.5, paste0("Mixture\nElement",i),
#              cex = 1.6, col = "black")
#       }else{
#         text(x = 0.5, y = 0.5, mixture.name[i],
#              cex = 1.6, col = "black")
#       }
#     } else if (i < j)
#     {
#       directory.1 <- paste(runs[i],"final_run/",sep="/")
#       directory.2 <- paste(runs[j],"final_run/",sep="/")
# 
#       run.1 <- strsplit(runs[i],"/")[[1]][5]
#       run.2 <- strsplit(runs[j],"/")[[1]][5]
#       file.1 <- paste0(tolower(run.1),"_Selection.csv")
#       file.2 <- paste0(tolower(run.2),"_Selection.csv")
#       if(run.1 != "Nterminus")
#       {
#         sec.str.1 <- strsplit(run.1,"_")[[1]]
#         cat.1 <- paste(sec.str.1[1],sec.str.1[2],sep="/")
#       } else
#       {
#         sec.str.1 <- "nterminus"
#         cat.1 <- "Nterminus"
#       }
#       if(run.2 != "Nterminus")
#       {
#         sec.str.2 <- strsplit(run.2,"_")[[1]]
#         cat.2 <- paste(sec.str.2[1],sec.str.2[2],sep="/")
#       } else
#       {
#         sec.str.2 <- "nterminus"
#         cat.2 <- "Nterminus"
#       }
#      #output <- paste0("Final_runs/Alpha/Structure_Downstream_Combo/",sec.str.1[1],"_",sec.str.1[2],"_",sec.str.2[1],"_",sec.str.2[2],"_exclude_CGG_CGA.pdf")
#     results<-demingApproach(directory.2,directory.1,file.2,file.1,optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
#     df <- results$df
#     plot(df$Posterior,df$Posterior.2,xlim=range(df[,c("X0.025.","X0.975.")]),ylim=range(df[,c("X0.025.2","X0.975.2")]),col="white")
# 
#     if (!(results$Slope.CI[1] < 1 && results$Slope.CI[2] > 1))
#     {
#       rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="red")
#     }
#     upper.panel.plot(results$df,results$Slope,results$Slope.CI,results$p.val)
#     slope.ci[index,] <- c(cat.1,cat.2,results$Slope,results$Slope.CI[1],results$Slope.CI[2])
#     index <- index + 1
#     }
#   }
# }
# 
# dev.off()



#agree<-compareRankings(c("Scer/Predicted/Results/Secondary_structures/Coil/final_run/Parameter_est/coil_Selection.csv",
#                                      "Scer/Predicted/Results/Secondary_structures/Helix/final_run/Parameter_est/helix_Selection.csv",
#                                     "Scer/Predicted/Results/Secondary_structures/Sheet/final_run/Parameter_est/sheet_Selection.csv"))                   
head.directory <- "Scer/Predicted/Results/"

# output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"coil_Selection.csv","helix_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Helix",file="Images_for_mike/scer_coil_vs_helix.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil vs. Sheet",file="Images_for_mike/scer_coil_vs_sheet.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Helix vs. Sheet",file="Images_for_mike/scer_helix_vs_sheet.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nOrdered vs. Disordered Regions",file="Images_for_mike/scer_ordered_vs_disorderd.pdf")
# 
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"coil_ord_Selection.csv","coil_dis_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil (Ordered Region)","Coil (Disordered Region)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nCoil (Ordered vs. Disordered)",file="Images_for_mike/scer_coil_ordered_vs_disordered.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"helix_ord_Selection.csv","helix_dis_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix (Ordered Region)","Helix (Disordered Region)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nHelix (Ordered vs. Disordered)",file="Images_for_mike/scer_helix_ordered_vs_disordered.pdf")
# 
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"coil_ord_Selection.csv","helix_ord_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil (Ordered Region)","Helix (Ordered Region)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nOrdered Coil vs. Ordered Helix",file="Images_for_mike/scer_ordered_coil_vs_ordered_helix.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_ord_Selection.csv","sheet_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil (Ordered Region)","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nOrdered Coil vs. Sheet",file="Images_for_mike/scer_ordered_coil_vs_sheet.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_ord_Selection.csv","sheet_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix (Ordered Region)","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nOrdered Helix vs. Sheet",file="Images_for_mike/scer_ordered_helix_vs_sheet.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"coil_dis_Selection.csv","helix_dis_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil (Disordered Region)","Helix (Disordered Region)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Disordered Helix",file="Images_for_mike/scer_disordered_helix_vs_disordered_coil_disordered.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_dis_Selection.csv","sheet_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil (Disordered Region)","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Coil vs. Sheet",file="Images_for_mike/scer_disordered_coil_vs_sheet.pdf")
# 
# output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_dis_Selection.csv","sheet_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix (Disordered Region)","Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\nDisordered Helix vs. Sheet",file="Images_for_mike/scer_disordered helix_vs_sheet.pdf")





# output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch","Disordered","final_run/"),"disordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered (Conserved Residues)","Disordered (Non-conserved Residues)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB in Disordered Regions:\nConserved vs. Non-conserved Residues",file="Images_for_mike/scer_disordered_conserved_vs_nonconserved_sacch.pdf")

# output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch","Ordered","final_run/"),"ordered_Selection.csv","ordered_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered (Conserved Residues)","Ordered (Non-conserved Residues)"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB in Ordered Regions:\nConserved vs. Non-conserved Residues",file="Images_for_mike/scer_ordered_conserved_vs_nonconserved_sacch.pdf")

# output.con <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved_sacch","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved_sacch","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output.con$df,b1=unname(output.con$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output.con$Slope.CI,title = "Comparing Selection on CUB:\nOrdered vs. Disordered (Conserved Residues)",file="Images_for_mike/scer_ordered_vs_disordered_conserved_sacch.pdf")

# output.unc <- demingApproach(file.path(head.directory,"Ordered_disordered_unconserved_sacch","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved_sacch","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = T,include.orig.ref = T)
# p <- plotResults(output.unc$df,b1=unname(output.unc$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output.unc$Slope.CI,title = "Comparing Selection on CUB:\nOrdered vs. Disordered (Non-conserved Residues)",file="Images_for_mike/scer_ordered_vs_disordered_nonconserved_sacch.pdf")

aa.groups <- list(c("C","D","E","F","H","I","K","N","Q","Y","Z"),
                  c("A","G","P","S","T","V"),
                  c("R","L"))


plots <- vector(mode="list",length=16)
for (i in 1:16)
{
  plots[[i]] <- vector(mode="list",length=3)
}
for (i in 1:length(aa.groups))
{
  if (i == 1)
  {
    codon.group <- "2"
  } else if(i == 2)
  {
    codon.group <- "4"
  } else{
    codon.group <- "6"
  }

  ## Structures
  output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Helix","final_run/"),"coil_Selection.csv","helix_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Helix"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structures/scer_coil_vs_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[1]][[i]] <- p



  output <- demingApproach(file.path(head.directory,"Secondary_structures","Coil","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structures/scer_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[2]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structures","Helix","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structures/scer_helix_vs_sheet",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[3]][[i]] <- p


  output <- demingApproach(file.path(head.directory,"Ordered_disordered","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Ordered_vs_disordered/scer_ordered_vs_disorderd_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[4]][[i]] <- p

  ## Ordered Conserved vs Non-conserved
  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved","Ordered","final_run/"),"ordered_Selection.csv","ordered_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered (Conserved Residues)","Ordered (Non-conserved Residues)"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Ordered_vs_disordered/scer_ordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[5]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved","Disordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved","Disordered","final_run/"),"disordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered (Conserved Residues)","Disordered (Non-conserved Residues)"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Ordered_vs_disordered/scer_disordered_conserved_vs_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[6]][[i]] <- p

  ## Ordered and Disordered conserved vs non-conserved
  output <- demingApproach(file.path(head.directory,"Ordered_disordered_conserved","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_conserved","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Ordered_vs_disordered/scer_ordered_vs_disorderd_conserved_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[7]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Ordered_disordered_unconserved","Ordered","final_run/"),file.path(head.directory,"Ordered_disordered_unconserved","Disordered","final_run/"),"ordered_Selection.csv","disordered_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered","Disordered"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Ordered_vs_disordered/scer_ordered_vs_disorderd_nonconserved_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[8]][[i]] <- p

  ## Ordered vs Disorderd Secondary Structures (Coil and Helix Only)
  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),"coil_ord_Selection.csv","coil_dis_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered Coil","Disordered Coil"),ci=T,bounds=output$Slope.CI,title =bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_ordered_coil_vs_disordered_coil_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[9]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"helix_ord_Selection.csv","helix_dis_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered Helix","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_ordered_helix_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[10]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),"coil_ord_Selection.csv","helix_ord_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered Coil","Ordered Helix"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_ordered_coil_vs_ordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[11]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_ord","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_ord_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_ordered_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[12]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),"coil_dis_Selection.csv","helix_dis_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Disordered Helix"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_disordered_coil_vs_disordered_helix_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[13]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Coil_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"coil_dis_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Coil","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_disordered_coil_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[14]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_ord","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_ord_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Ordered Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_ordered_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[15]][[i]] <- p

  output <- demingApproach(file.path(head.directory,"Secondary_structure_order","Helix_dis","final_run/"),file.path(head.directory,"Secondary_structures","Sheet","final_run/"),"helix_dis_Selection.csv","sheet_Selection.csv",optimal.as.ref = F,normalize.deta = F,include.orig.ref = F,include.AA=aa.groups[[i]],rescale.agct=T)
  p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Disordered Helix","Sheet"),ci=T,bounds=output$Slope.CI,title = bquote(atop("Comparison of Selection Coefficient -"*Delta*eta,.(codon.group)*" Codon AA")),file=paste0("~/CUB_Secondary_Structure/Updated_Images_Mike_Presentation/Secondary_structure_and_disorder/scer_disordered_helix_vs_sheet_",codon.group,"_codon.pdf"),range.xy=NULL)
  plots[[16]][[i]] <- p
}


# titles <- c( substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Coil vs. Helix",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Coil vs. Sheet",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Helix vs. Sheet",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered vs. Disordered",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Conserved Residues in Ordered Regions",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Conserved Residues in Disordered Regions",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered vs. Disordered (Conserved Residues)",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered vs. Disordered (Non-conserved Residues)",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered Coil vs. Disordered Coil",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered Helix vs. Disordered Helix",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered Coil vs. Ordered Helix",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered Coil vs. Sheet",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Disordered Coil vs. Disordered Helix",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Disordered Coil vs. Sheet",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Ordered Helix vs. Sheet",sep="")),
#              substitute(paste("Comparison of Selection Coefficients -",Delta,eta,": Disordered Helix vs. Sheet",sep=""))
             
# )

# for (i in 1:16)
# {
#   subplots <- plots[[i]]
#   comb.plot <- plot_grid(subplots[[1]], subplots[[2]], subplots[[3]],labels = c("A","B","C"),ncol=3,label_size=16)
#  # label <- titles[i]
#   title <- ggdraw() + draw_label(
#     titles[[i]],
#       fontface = 'bold',
#       x = 0,
#       hjust =0,size=25
#     ) +
#     theme(
#       # add margin on the left of the drawing canvas,
#       # so title is aligned with left edge of first plot
#       plot.margin = margin(0, 0, 0, 5)
#     )
#   comb<-plot_grid(
#     title, comb.plot,
#     ncol = 1,
#     # rel_heights values control vertical title margins
#     rel_heights = c(0.2, 1)
#   )
#   ggsave2(file=paste0(i,".pdf"),plot=comb,width=14,height=7)
# }

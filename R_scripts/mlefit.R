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


demingApproach <- function(directory.1,directory.2,file.1,file.2,exclude=c(),optimal.as.ref = F,normalize.deta=T,include.orig.ref=T)
{
  names.aa <- aminoAcids()
  # sd.1 <- vector(mode = "numeric",length=59 - length(exclude))
  # sd.2 <- vector(mode="numeric",length=59 - length(exclude))
  # index <- 1
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
  if (length(exclude) != 0)
  {
    sel.1 <- sel.1[-c(which(sel.1[,2] %in% exclude)),]
    sel.2 <- sel.2[-c(which(sel.2[,2] %in% exclude)),]
  }
  if (include.orig.ref)
  {
    for(aa in names.aa)
    {
      removed <- 0
      if(aa == "M" || aa == "W" || aa == "X") next
      sd.1.tmp <- sel.1[which(sel.1$AA == aa),"Std.Dev"]^2
      sd.2.tmp <- sel.2[which(sel.2$AA == aa),"Std.Dev"]^2
      ## Changing variance redistribution to not force equal weighting, but non-observed Delta.eta still assumed to have variance equal to the 
      ## average across all variances 
      # ref.1 <- which(sd.1.tmp == 0)
      # ref.2 <- which(sd.2.tmp == 0)
      total.sd.1 <- sum(sd.1.tmp)
      total.sd.2 <- sum(sd.2.tmp)
      # weights.1 <- sd.1.tmp/total.sd.1 * ((length(sd.1.tmp)-1)/length(sd.1.tmp))
      # weights.2 <- sd.2.tmp/total.sd.2 * ((length(sd.2.tmp)-1)/length(sd.2.tmp))
      # weights.1[ref.1] <- 1/length(sd.1.tmp)
      # weights.2[ref.2] <- 1/length(sd.2.tmp)
      sd.1.tmp <- sqrt(rep(total.sd.1/length(sd.1.tmp),length(sd.1.tmp)))
      sd.2.tmp <- sqrt(rep(total.sd.2/length(sd.2.tmp),length(sd.2.tmp)))
      #sd.2.tmp <- sqrt(total.sd.2 * weights.2)
      sd.1 <- c(sd.1,sd.1.tmp)
      sd.2 <- c(sd.2,sd.2.tmp)
    }
    if (normalize.deta)
    {
      sel.1 <- normalize(sel.1)
      sel.2 <- normalize(sel.2)
    }
  } else{
    sel.1 <- sel.1[which(sel.1$Posterior != 0),]
    sel.2 <- sel.2[which(sel.2$Posterior != 0),]
    sd.1 <- sel.1$Std.Dev
    sd.2 <- sel.2$Std.Dev
  }
  df <- plyr::join(sel.1,sel.2,by=c("AA","Codon"))
  colnames(df)[7] <- "Posterior.2"
  colnames(df)[8] <- "Std.Dev.2"
  colnames(df)[9] <- "X0.025.2"
  colnames(df)[10] <- "X0.975.2"
  reg <- deming(Posterior.2 ~ Posterior+0,data=df,xstd=sd.1,ystd=sd.2)
  b1 <- reg$coefficients[2]
  ci.b1 <- reg$ci[2,]
  p.val <- pnorm(abs((b1-1)/sqrt(reg$variance[2,2])),lower.tail = F)*2
  return(list(df=df,Slope=b1,Slope.CI=ci.b1,SD.1=sd.1,SD.2=sd.2,p.val=p.val))
}

plotResults<-function(data,b1,b0,categories=c("X","Y"),file="title",title="Regression",mark=NULL)
{
  
  if (!is.null(bounds))
  {
    conf.int.1 <- bounds[1]
    conf.int.2 <- bounds[2]
    print(conf.int.2)
    print(conf.int.1)
    print(b1)
    l <- data.frame(s=c(b1,conf.int.2,conf.int.1,1.0),ic=c(b0,0.0,0.0,0.0),Line=c("Deming Regression","97.5% CI","2.5% CI","y=x"),stringsAsFactors = F)
    l$Line <- factor(l$Line,levels=c("Deming Regression","97.5% CI","2.5% CI","y=x"))
    levels(l$Line) <- c("Deming Regression","95% CI","95% CI","y=x")
    legend.colors <- c("black","grey","black")
    lines.reg <- c("solid","dashed","dashed")
    
  } else{
    l <- data.frame(s=c(b1,1.0),ic=c(b0,0.0),Line=c("Model II Regression","y = x"))
    legend.colors <- c("black","black")
    lines.reg <- c("solid","dashed")
  }
  p <- ggplot(data,aes(Posterior,Posterior.2))
  p <-(p + geom_point(colour="black",size=4)
       + labs(x=bquote(.(categories[1])~"("*Delta*eta*")"),y=bquote(.(categories[2])~"("*Delta*eta*")")) 
       + geom_abline(data=l,mapping=aes(slope=s,intercept=ic,linetype=Line,color=Line))
       + scale_color_manual(values=legend.colors)
       + scale_linetype_manual(values=lines.reg)
       + ggtitle(label=title))
  if (!is.null(mark))
  {
    p<- p + geom_text(aes(label=ifelse(data$Codon %in% mark,as.character(data$Codon),'')),nudge_x =-0.2,nudge_y=0.1)
  }
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),3)
  if (ci)
  {
    p <- (p + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2),alpha=0.2,color="black") 
          + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.),alpha=0.2,color="black"))
    range.xy <- range(c(data[,5:6],data[,8:9]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  } else{
    range.xy <- range(c(data[,3],data[,7]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  }
  p <- p + scale_x_continuous(limits = range.xy) + scale_y_continuous(limits = range.xy)
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(format(rho.p,nsmall=3)))
  b1 <- round(b1,3)
  b0 <- round(b0,3)
  if (b0 != 0)
  {
    eq.exp <- bquote("y ="~ .(format(b1,nsmall=3))*"x +"~ .(format(b0,nsmall=3)))
  } else{
    eq.exp <- bquote("y ="~ .(format(b1,nsmall=3))*"x")
  }
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.10,label=deparse(cor.exp),parse=T,size=5) 
  p <- p + annotate("text",x=xlim[2]-width*0.2,y=ylim[1]+height*0.20,label=deparse(eq.exp),parse=T,size=5)
  p <- (p + theme_bw()
        + theme(axis.title=element_text(size = 14,face="bold"),axis.text=element_text(size=14))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.25,0.8),legend.title=element_blank(),legend.text=element_text(size=14),plot.title = element_text(hjust = 0.5,size=14)))
  
  
  ggsave(filename = file,device="pdf",width = 7,height = 7)
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


compareAcrossCodon <- function(runs.to.compare,legend.info,colors=NULL)
{
  if (length(runs.to.compare) != length(legend.info)) stop("runs.to.compare and legend.info should be same length")
  if (is.null(colors))
  {
    colors <- c("black","red","blue","green","purple")
  }
  
  dfs <- lapply(runs.to.compare,FUN=read.csv,header=T,stringsAsFactors=F)
  max.element <- max(unlist(lapply(dfs,FUN=function(x){max(x[,6])})))
  min.element <- min(unlist(lapply(dfs,FUN=function(x){min(x[,5])})))
  ylim <- c(min.element - 0.05,max.element+0.05)
  plot(x=NULL,y=NULL,xlim=c(1,59),ylim=ylim,xlab="Codon",ylab="Delta.eta")
  for (i in 1:length(dfs))
  {
    text(x=1:59,dfs[[i]][,3],col=colors[i],labels=dfs[[i]][,1])
    arrows(x0=1:59,y0=dfs[[i]][,5],y1=dfs[[i]][,6],length = 0,col=colors[i])
  }
  legend(x = "topleft",legend = legend.info,col = colors[1:length(runs.to.compare)],pch=c(1,1))
  abline(h=0,lty=2)
}

compareAcrossCodon.2 <- function(run.1,run.2,legend.info,aa="all",...)
{
  df.1 <- read.csv(run.1,header=T,stringsAsFactors = F)
  df.2 <- read.csv(run.2,header=T,stringsAsFactors = F)
  if (aa == "hydrophobic")
  {
    df.1 <- df.1[which(df.1$AA %in% c("A","I","L","F","Y","V")),]
    df.2 <- df.2[which(df.2$AA %in% c("A","I","L","F","Y","V")),]
  } else if (aa == "positive"){
    df.1 <- df.1[which(df.1$AA %in% c("R","H","K")),]
    df.2 <- df.2[which(df.2$AA %in% c("R","H","K")),]
  } else if (aa == "negative"){
    df.1 <- df.1[which(df.1$AA %in% c("D","E")),]
    df.2 <- df.2[which(df.2$AA %in% c("D","E")),]
  } else if (aa == "charged"){
    df.1 <- df.1[which(df.1$AA %in% c("D","E","R","H","K")),]
    df.2 <- df.2[which(df.2$AA %in% c("D","E","R","H","K")),]
  } else if (aa == "polar"){
    df.1 <- df.1[which(df.1$AA %in% c("Z","S","T","N","Q")),]
    df.2 <- df.2[which(df.2$AA %in% c("Z","S","T","N","Q")),]
  } else if (aa == "hydrophilic"){
    df.1 <- df.1[which(df.1$AA %in% c("D","E","R","H","K","Z","S","T","N","Q")),]
    df.2 <- df.2[which(df.2$AA %in% c("D","E","R","H","K","Z","S","T","N","Q")),]
  } 
  else{
    df.1 <-df.1
    df.2 <- df.2
  }
  ylim <- c(min(c(df.1[,5],df.2[,5]))-0.05,max(c(df.1[,6],df.2[,6]))+0.05)
  cur.aa <- "A"
  labels <- c()
  splits <- c()
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
  split.line <- which(labels == "")
  placeholder <- data.frame(AA="",Codon="",Posterior=NA,Std.Dev=NA,X0.025.=NA,X0.975.=NA)
  for (i in split.line)
  {
    df.1 <- rbind(df.1[1:i-1,],placeholder,df.1[i:nrow(df.1),])
    df.2 <- rbind(df.2[1:i-1,],placeholder,df.2[i:nrow(df.2),])
  }
  #par(xpd=TRUE)
  par(mar=c(5.1,4.1,5.1,9.1))
  aa.spots <- c(0,split.line)+diff(c(0,split.line))/2
  plot(x=1:nrow(df.1),y=df.1[,3],ylim=ylim,xaxt='n',xlab="Codon",ylab=expression(Delta*eta),...)
  points(x=1:nrow(df.2),df.2[,3],col="red")
  arrows(x0=1:nrow(df.1),y0=df.1[,5],y1=df.1[,6],length = 0,col="black")
  arrows(x0=1:nrow(df.2),y0=df.2[,5],y1=df.2[,6],length = 0,col="red")
  abline(h=0,lty=2)
  axis(side=1,at=1:length(labels),labels=labels,las=2,cex.axis=0.75)
  axis(side=3,at=aa.spots,labels=unique(df.1[which(df.1$AA!=""),"AA"]),cex.axis=0.75)
  legend(x = "topright",inset=c(-0.22,0),xpd=T,legend = legend.info,col = c("black","red"),pch=c(1,1))
  abline(v = split.line,lty=2,col="lightblue")
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
library(ggrepel)
head.directory <- "For_paper/Results/Merged/"
# output.100 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"100","final_run/"),"regionA_Selection.csv","100_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# output.50 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"50","final_run/"),"regionA_Selection.csv","50_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# output.25 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"25","final_run/"),"regionA_Selection.csv","25_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# output.10 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"10","final_run/"),"regionA_Selection.csv","10_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
# output.5 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"5","final_run/"),"regionA_Selection.csv","5_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)
output.1 <- demingApproach(file.path(head.directory,"RegionA","final_run/"),file.path(head.directory,"1","final_run/"),"regionA_Selection.csv","1_Selection.csv",optimal.as.ref = T,normalize.deta = T,include.orig.ref = T)

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
#      + ggtitle(label="Expected Results when Regions Under Different\nSelective Pressures"))
# p <- p + scale_x_continuous(limits = c(-1,1)) + scale_y_continuous(limits = c(-1,1))
# p <- p + guides(colour=guide_legend(title="%Region B under\nSeleciton for\nInefficiency"),
#                 linetype=guide_legend(title="%Region B under\nSeleciton for\nInefficiency"))
# p <- (p + theme_bw()
#       + theme(axis.title=element_text(size = 14,face="bold"),axis.text=element_text(size=14))
#       + theme(axis.line = element_line(colour = "black"))
#       + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#       + theme(legend.text=element_text(size=14),plot.title = element_text(hjust = 0.5,size=14)))

#ggsave(p,"../../CUB_Secondary_Structure/Images/simulated_expectation.pdf")
# pdf("../../CUB_Secondary_Structure/Images/simulated_expectation.pdf",height=7,width=7)
# p
# dev.off()
#p <- plotResults(output$df,b1=unname(output$Slope),b0 = 0,categories = c("Coil/Turn","Helix/Sheet"),ci=T,bounds=output$Slope.CI,title = "Comparing Selection on CUB:\n Coil/Turn vs. Helix/Sheet",file="Images_for_paper/scer_coil_turn_helix_sheet_exp_conservative_1022.pdf")
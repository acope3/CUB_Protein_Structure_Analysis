## Author: Alexander Cope
## Script used to rescale selection estimates $\Delta\eta$
## rescaleByMean was the primary function used in the publication for Deming Regression and subsequent plots. 

library(AnaCoDa)


changeCodon <- function(codon)
{
  nucleotides <- unlist(strsplit(codon,""))
  if (nucleotides[3] == "A")
  {
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|G",collapse="")
  } else if (nucleotides[3] == "C"){
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|T",collapse="")
  }
  return(codon)
}

changeCodon_ATCG <- function(codon)
{
  nucleotides <- unlist(strsplit(codon,""))
  if (nucleotides[3] == "A")
  {
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|T",collapse="")
  } else if (nucleotides[3] == "C"){
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|G",collapse="")
  }
  return(codon)
}



rescaleATCG <- function(sel.file,parameter.file)
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
      trace.cgc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGC",1,T)
      trace.cgg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGG",1,T)
      trace.cgc <- trace.cgc- trace.cgg
      data.tmp[which(data.tmp$Codon == "CGC"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(mean(trace.cgc),sd(trace.cgc),quantile(trace.cgc,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "AGA"),] <- NA
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
      trace.ctc <- trace.ctc - trace.ctg
      trace.cta <- trace.cta - trace.ctt
      data.tmp[which(data.tmp$Codon == "CTC"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(mean(trace.ctc),sd(trace.ctc),quantile(trace.ctc,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "CTA"),c("Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(mean(trace.cta),sd(trace.cta),quantile(trace.cta,c(0.025,0.975)))
      data.tmp[which(data.tmp$Codon == "CTT"),] <- NA
      data.tmp[which(data.tmp$Codon == "CTG"),] <- NA
      data.tmp[which(data.tmp$Codon == "TTA"),] <- NA
      data.tmp[which(data.tmp$Codon == "TTG"),] <- NA
     }
     if (a == "I")
     {
      data.tmp[which(data.tmp$Codon == "ATC"),] <-NA
      data.tmp[which(data.tmp$Codon == "ATT"), ] <- NA
     }
     if (nrow (data.tmp) == 4)
     {
      codon.c <- data.tmp[2,"Codon"]
      codon.g <- data.tmp[3,"Codon"]
      trace.c <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.c,1,T)
      trace.g <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.g,1,T)
      trace.c <- trace.c - trace.g 
      data.tmp[,c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(mean(trace.c),sd(trace.c),quantile(trace.c,c(0.025,0.975)))
      data.tmp[3,] <- NA
      data.tmp[4,] <- NA
     }
     if (nrow (data.tmp) == 2)
     {
      data.tmp[1,] <- NA
      data.tmp[2,] <- NA
     }
     data[which(data$AA == a),] <- data.tmp
   }
   data <- data[which(!is.na(data$Posterior)),]
   return(data)
}




rescaleAGCT <- function(sel.file,parameter.file)
{
	data <- read.table(sel.file,sep=",",header=T,stringsAsFactors=F,row.names = 2)
  new.data <- data.frame("AA"=c(),"Posterior"=c(),"Std.Dev"=c(),"X0.025."=c(),"X0.975."=c())
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
      trace.cgc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGC",1,T)
      trace.cgg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CGG",1,T)
      trace.aga <- trace.aga - trace.agg
      trace.cga <- trace.cga - trace.cgg
      trace.cgc <- trace.cgc- trace.cgg
      data.tmp["AGA|G",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.aga),sd(trace.aga),quantile(trace.aga,c(0.025,0.975)))
      data.tmp["CGA|G",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.cga),sd(trace.cga),quantile(trace.cga,c(0.025,0.975)))
      data.tmp["CGC|G",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.cgc),sd(trace.cgc),quantile(trace.cgc,c(0.025,0.975)))
      data.tmp["CGA|T",] <- data.tmp["CGA",]
      data.tmp["CGC|T",] <- data.tmp["CGC",]
      data.tmp["AGG",] <- NA
      data.tmp["CGG",] <- NA
      data.tmp["CGT",] <- NA
      data.tmp["CGA",] <- NA
      data.tmp["AGA",] <- NA
      data.tmp["CGC",] <- NA    
     }
     if (a == "L")
     {
      trace.ctc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTC",1,T)
      trace.ctt <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTT",1,T)
      trace.cta <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTA",1,T)
      trace.ctg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,"CTG",1,T)
      trace.ctc.t <- trace.ctc - trace.ctt
      trace.cta.g <- trace.cta - trace.ctg
      trace.ctc.g <- trace.ctc - trace.ctg
      trace.cta.t <- trace.cta - trace.ctt
      data.tmp["CTC|T",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(a,mean(trace.ctc.t),sd(trace.ctc.t),quantile(trace.ctc.t,c(0.025,0.975)))
      data.tmp["CTA|G",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(a,mean(trace.cta.g),sd(trace.cta.g),quantile(trace.cta.g,c(0.025,0.975)))
      data.tmp["CTC|G",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(a,mean(trace.ctc.g),sd(trace.ctc.g),quantile(trace.ctc.g,c(0.025,0.975)))
      data.tmp["CTA|T",c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <-  c(a,mean(trace.cta.t),sd(trace.cta.t),quantile(trace.cta.t,c(0.025,0.975)))
      data.tmp["TTA|G",] <- data.tmp["TTA",]
      data.tmp["CTA",] <- NA
      data.tmp["CTC",] <- NA
      data.tmp["CTG",] <- NA
      data.tmp["CTT",] <- NA
      data.tmp["TTA",] <- NA
      data.tmp["TTG",] <- NA
     }
     if (a == "I")
     {
      data.tmp["ATA|T",] <- data.tmp["ATA",]
      data.tmp["ATC|T",] <- data.tmp["ATC",]
      data.tmp["ATA", ] <- NA
      data.tmp["ATC", ] <- NA
      data.tmp["ATT", ] <- NA
     }
     if (nrow (data.tmp) == 4)
     {
      codon.a <- row.names(data.tmp)[1]
      codon.c <- row.names(data.tmp)[2]
      codon.g <- row.names(data.tmp)[3]
      codon.t <- row.names(data.tmp)[4]
      trace.a <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.a,1,T)
      trace.c <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.c,1,T)
      trace.g <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codon.g,1,T)
      trace.a.g <- trace.a - trace.g
      trace.c.g <- trace.c - trace.g 
      data.tmp[changeCodon(codon.a),c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.a.g),sd(trace.a.g),quantile(trace.a.g,c(0.025,0.975)))
      data.tmp[changeCodon_ATCG(codon.c),c("AA","Posterior","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.c.g),sd(trace.c.g),quantile(trace.c.g,c(0.025,0.975)))
      data.tmp[changeCodon_ATCG(codon.a),] <- data.tmp[codon.a,]
      data.tmp[changeCodon(codon.c),] <- data.tmp[codon.c,]
      data.tmp[codon.a,] <- NA
      data.tmp[codon.c,] <- NA
      data.tmp[codon.g,] <- NA
      data.tmp[codon.t,] <- NA
     }
     if (nrow (data.tmp) == 2)
     {
      codon.1 <- row.names(data.tmp)[1]
      codon.2 <- row.names(data.tmp)[2]
      data.tmp[changeCodon(codon.1),] <- data.tmp[codon.1,]
      data.tmp[codon.1,] <- NA
      data.tmp[codon.2,] <- NA
     }
     data.tmp <- data.tmp[which(!is.na(data.tmp$Posterior)),]
     #data[which(data$AA == a),] <- data.tmp
     new.data <- rbind(new.data,data.tmp)
   }
   
   new.data[,"Codon"] <- row.names(new.data)
   #print(colnames(new.data))
   new.data <- new.data[,c("AA","Codon","Posterior","Std.Dev","X0.025.","X0.975.")]
   return(new.data)
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
    # print(data.tmp)
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
     #print(data.tmp)
     data.tmp[which(data.tmp$Mean!=0),c("Mean","X2.5.","X97.5.")] <- data.tmp[which(data.tmp$Mean!=0),c("Mean","X2.5.","X97.5.")] - mean.deta
     data.tmp[which(data.tmp$Mean==0),c("Mean","Std.Dev","X2.5.","X97.5.")] <- c(-mean.deta,sd.trace,quantile(0-mean.trace,c(0.025,0.975)))
     print(data.tmp)
     #data.tmp[which(data.tmp$Posterior!=0),c("Posterior","X0.025.","X0.975.")] <- data.tmp[which(data.tmp$Posterior!=0),c("Posterior","X0.025.","X0.975.")] - mean.deta
     #data.tmp[which(data.tmp$Posterior==0),c("Posterior","Std.Dev","X0.025.","X0.975.")] <- c(-mean.deta,sd.trace,quantile(0-mean.trace,c(0.025,0.975)))
     
     data[which(data$AA == a),] <- data.tmp

   }

   colnames(data) <- c("AA","Codon","Posterior","Std.Dev","0.025%","0.975%")
   data <- data[,c("AA","Codon","Posterior","Std.Dev","0.025%","0.975%")]
   return(data)
}

rescaleByGenomeOptimal<-function(sel.file,genome.file="selection_mod_scerevisiae.csv")
{
  data <- read.table(sel.file,sep=",",header=T,stringsAsFactors=F)
  data <- data[,which(colnames(data) != "Effective.Samples")]
  colnames(data) <- c("AA","Codon","Posterior","Std.Dev","0.025%","0.975%")
    
  #data <- data[which(data$Posterior != 0),]
  ref.set <- read.table(genome.file,sep=",",header=T,stringsAsFactors=F)
  print(ref.set)
  print(data)
  aa <- aminoAcids()
  for (a in aa)
  {
   if (a == "X" || a=="M" || a == "W") next
   codons.w.ref <- AAToCodon(a,F)
   data.tmp <- data[which(data$AA == a),]
   ref.set.aa <- ref.set[which(ref.set$AA == a ),]
   current.reference.row <- which(data.tmp[,"Posterior"]==0)
   genome.opt.index <- which.min(c(ref.set.aa[,3],0))
   ##optimal codon is the reference, nothing to do
   if (genome.opt.index != length(codons.w.ref))
   {
      data.tmp[,c("Posterior","0.025%","0.975%")] <- data.tmp[,c("Posterior","0.025%","0.975%")] - data.tmp[genome.opt.index,"Posterior"]

      ##Get row of the optimal codon, which should be 0
      optimal.codon.row <- which(data.tmp[,"Posterior"]==0.0)
      data.tmp[current.reference.row,"Std.Dev"] <- data.tmp[genome.opt.index,"Std.Dev"]
      
      data.tmp[current.reference.row,c("0.025%","0.975%")] <- data.tmp[genome.opt.index,c("0.025%","0.975%")] + data.tmp[current.reference.row,"Posterior"]
      data.tmp[optimal.codon.row,c("Posterior","Std.Dev","0.025%","0.975%")] <- 0.0
   }
   
   #data.tmp <- data.tmp[which(data.tmp$Posterior != 0),]
   data[which(data$AA == a),] <- data.tmp
  }
  return(data)
}

head.directory <- "../Scer/Simulated/Results/"
genome.file <- "../selection_mod_scerevisiae.csv"
targets <- c("Secondary_structures","Ordered_disordered")

for (f in targets)
{
	structure.loc <- file.path(head.directory,f)
  print(structure.loc)
	structures <- list.dirs(structure.loc,recursive=F)
  print(structures)
	for (struct in structures)
	{
		sel.file.loc <- file.path(struct,"restart_5","Parameter_est/")
		sel.file <- list.files(sel.file.loc,pattern="*_Selection.csv",full.names=T)
		parameter.file <- file.path(struct,"restart_5","R_objects","parameter.Rda")
    print(sel.file.loc)
    sel.rescaled <- rescaleByGenomeOptimal(sel.file,genome.file=genome.file)
    write.table(sel.rescaled,file.path(sel.file.loc,"selection_rescaled_to_genome_optimal.csv"),sep=",",col.names=T,row.names=F,quote=F)
	}	
}



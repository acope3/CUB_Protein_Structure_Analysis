library(AnaCoDa)


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



head.directory <- "../Scer/Predicted/Results/"

targets <- c("Secondary_structures_conserved","Secondary_structures_unconserved")

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



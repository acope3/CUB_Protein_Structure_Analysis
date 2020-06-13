## An implementation of the model from Shah and Gilchrist Plos Genetics 2010
## Given a list of tRNA for a species, will generate list of cognates and near-cognates for each codon
## This can then be used to generate estimates of elongation rates and missense error rates


library(ggplot2)


rescaleAGCT <- function(data,metric="e_m")
{
    aa <- unique(data$AA)
    for (a in aa)
    {
     data.tmp <- data[which(data$AA == a),]
     if (a == "R")
     {
      data.tmp[which(data.tmp$Codon == "AGA"),metric] <- data.tmp[which(data.tmp$Codon == "AGA"),metric] - data.tmp[which(data.tmp$Codon == "AGG"),metric]
      data.tmp[which(data.tmp$Codon == "CGA"),metric] <- data.tmp[which(data.tmp$Codon == "CGA"),metric] - data.tmp[which(data.tmp$Codon == "CGG"),metric]
      data.tmp[which(data.tmp$Codon == "CGC"),metric] <- data.tmp[which(data.tmp$Codon == "CGC"),metric] - data.tmp[which(data.tmp$Codon == "CGT"),metric]
      data.tmp[which(data.tmp$Codon == "AGG"),] <- NA
      data.tmp[which(data.tmp$Codon == "CGG"),] <- NA
      data.tmp[which(data.tmp$Codon == "CGT"),] <- NA
     }
     if (a == "L")
     {
      data.tmp[which(data.tmp$Codon == "CTC"),metric] <-  data.tmp[which(data.tmp$Codon == "CTC"),metric] - data.tmp[which(data.tmp$Codon == "CTT"),metric]
      data.tmp[which(data.tmp$Codon == "CTA"),metric] <-  data.tmp[which(data.tmp$Codon == "CTA"),metric] - data.tmp[which(data.tmp$Codon == "CTG"),metric]
      data.tmp[which(data.tmp$Codon == "TTA"),metric] <-  data.tmp[which(data.tmp$Codon == "TTA"),metric] - data.tmp[which(data.tmp$Codon == "TTG"),metric]
      data.tmp[which(data.tmp$Codon == "CTT"),] <- NA
      data.tmp[which(data.tmp$Codon == "CTG"),] <- NA
      data.tmp[which(data.tmp$Codon == "TTG"),] <- NA
     }
     if (a == "I")
     {

      data.tmp[which(data.tmp$Codon == "ATA"),] <-NA
      data.tmp[which(data.tmp$Codon == "ATC"),metric] <-data.tmp[which(data.tmp$Codon == "ATC"),metric] - data.tmp[which(data.tmp$Codon == "ATT"),metric]
      data.tmp[which(data.tmp$Codon == "ATT"), ] <- NA
     }
     if (nrow (data.tmp) == 4)
     {
      data.tmp[1,metric] <- data.tmp[1,metric] - data.tmp[3,metric]
      data.tmp[2,metric] <- data.tmp[2,metric] - data.tmp[4,metric]
      data.tmp[3,] <- NA
      data.tmp[4,] <- NA
     }
     if (nrow (data.tmp) == 2)
     {
      data.tmp[1,metric] <- data.tmp[1,metric] - data.tmp[2,metric]
     }
     data[which(data$AA == a),] <- data.tmp
   }
   data <- data[which(!is.na(data$AA)),]
   return(data)
}



reverse.complement <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  reverse <- rev(seq.split)
  trna <- unlist(lapply(reverse,complement.rules))
  trna <- paste(trna,collapse='')
  return(trna)

}


complement.rules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return("T")
  } else if (nucleotide == "C"){
    return("G")
  } else if(nucleotide == "G"){
    return("C")
  } else if (nucleotide == "T"){
    return("A")
  }
}

wobble.rules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return(c("T","C","G","A"))
  } else if (nucleotide == "C"){
    return(c("C"))
  } else if(nucleotide == "G"){
    return(c("G","A"))
  } else if (nucleotide == "T"){
    return(c("A","G","C","T"))
  }
}

## assumes seq is the trna sequences ordered from 5'-end to 3'-end (will match codon from 3'-end to 5'-end)
wobble.neighbors <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  wobbles <- wobble.rules(seq.split[1])
  second <- rep(seq.split[2],length(wobbles))
  third <- rep(seq.split[3],length(wobbles))
  wobble.trna <- paste0(wobbles,second,third)
  return(wobble.trna)
}

get.neighbors <- function(codon)
{
  neighbors <- character(length=12)
  nuc <- unlist(strsplit(codon,split='',fixed=T))
  neighbors[1:4] <- paste0(nuc[1],c("A","G","T","C"),nuc[3])
  neighbors[5:8] <- paste0(c("A","G","T","C"),nuc[2],nuc[3])
  neighbors[9:12] <- paste0(nuc[1],nuc[2],c("A","G","T","C"))
  return(neighbors)
}

check.same.codons<- function(set.1,set.2)
{
  set.1 <- unlist(strsplit(set.1,","))
  set.2 <- unlist(strsplit(set.2,","))
  if (length(set.1) != length(set.2))
  {
    return(F)
  } else{
  disagree <- setdiff(set.1,set.2)
  if (length(disagree) >0)
  {
    return(F)
  }
  }
  return(T)
}

## This function is only to make sure the list of cognates, pseudo-cognates, and near-cognates from the
## published work matches with published literature (it does)
sanity.check <- function()
{
  mine <- read.table("ecoli_cognate_near_cognate.tsv",sep="\t",header=T,stringsAsFactors=F)
  premal <- read.table("~/CUB_structure_ecoli/shah_gilchrist_2010_ecoli.tsv.csv",sep="\t",header=T,stringsAsFactors=F)
  for (i in 1:nrow(mine))
  {
    codon <- mine[i,"Codon"]
    my.row <- mine[i,]
    premal.row <- premal[which(premal$Codon == codon),]
    check.cognates <- check.same.codons(my.row$Cognates,premal.row$Cognates)
    check.pseudo <- check.same.codons(my.row$Pseudo,premal.row$Pseudo.cognates)
    check.near <- check.same.codons(my.row$Near,premal.row$Near.cognates)
    if (!(check.cognates && check.pseudo && check.near))
    {
      print(paste("Codon",codon,"is does not match!"))
    }
  }
}



focal.vs.neighbors <- function()
{

  trna <- read.table("../scer_tRNA_not_adjusted.tsv",sep="\t",header=T,stringsAsFactors = F)
  #trna$tRNA <- ifelse(trna$tRNA%%1 == 0,trna$tRNA,0)
  codons <- trna$Codon
  trna[,"Neighbors.trna"] <- rep(0,nrow(trna))
  for (i in codons)
  {
    nuc <- unlist(strsplit(i,split='',fixed=T))
    neighbors <- get.neighbors(i)
    codon.neighbors <- trna[which((trna$Codon %in% neighbors) & (trna$Codon != i)),"tRNA"]
    trna[which(trna$Codon == i),"Neighbors.trna"] <- sum(codon.neighbors)
  }

  deg <- c(2,4,6)
  counts <- table(trna$AA)

  xlim <- range(trna$tRNA)
  ylim <- range(trna$Neighbors.trna)
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  pdf("../Missense_errors/scer_focal_abundance_vs_neighbor_abundance.pdf")
  plot(trna$tRNA,trna$Neighbors.trna,xlab="Focal tRNA Abundance",ylab="Neighbor tRNA Abundance",main="S. cerevisiae")
  corr <- cor.test(trna$tRNA,trna$Neighbors.trna,method="spearman")
  text(x=xlim[2]-0.1*width,y=ylim[1]+0.1*height,label=bquote(rho~"="~.(round(corr$estimate,3))))
  text(x=xlim[2]-0.1*width,y=ylim[1]+0.05*height,label=bquote("p ="~.(corr$p.value)))

  for (i in deg)
  {

    amino <- names(which(counts == i))
    trna.aa <- trna[which(trna$AA %in% amino),]
    xlim <- range(trna.aa$tRNA)
    ylim <- range(trna.aa$Neighbors.trna)
    width <- xlim[2] - xlim[1]
    height <- ylim[2] - ylim[1]
    plot(trna.aa$tRNA,trna.aa$Neighbors.trna,xlab="Focal tRNA Abundance",ylab="Neighbor tRNA Abundance",main=paste0("S. cerevisiae:\n",i," Codon AA"))
    corr <- cor.test(trna.aa$tRNA,trna.aa$Neighbors.trna,method="spearman")
    text(x=xlim[2]-0.1*width,y=ylim[1]+0.1*height,label=bquote(rho["s"]~"="~.(round(corr$estimate,3))))
    text(x=xlim[2]-0.1*width,y=ylim[1]+0.05*height,label=bquote("p ="~.(corr$p.value)))
  }

  dev.off()
}


## Given list of tRNA, identify cognate, pseudo-cognate, and near-cognate tRNA for all codons
getCognate <- function()
{

  trna <- read.table("../scer_tRNA_not_adjusted.tsv",sep="\t",header=T,stringsAsFactors = F)
  trna[,"AntiCodon"] <- unlist(lapply(trna$Codon,reverse.complement))
  trna.present <- trna[which(trna$tRNA > 0),]
  anticodons <- trna.present$AntiCodon
  df <- data.frame(Codon = trna$Codon,AA=trna$AA,Cognates=character(length(nrow(trna))),Pseudo=character(length(nrow(trna))),Near=character(length(nrow(trna))),stringsAsFactors=F)
  for (i in 1:length(anticodons))
  {
     aa <- trna.present[i,"AA"]
     wobbles <- unlist(lapply(wobble.neighbors(anticodons[i]),reverse.complement))
     for (j in wobbles)
     {
        if (trna[which(trna$Codon== j),"AA"] == aa)
        {
          df[which(df$Codon == j),"Cognates"] <- paste(unlist(df[which(df$Codon == j),"Cognates"]),anticodons[i],sep=",")
        } else{
          df[which(df$Codon == j),"Near"] <- paste(unlist(df[which(df$Codon == j),"Near"]),anticodons[i],sep=",")
        }
     }
     
  }

  codons <- trna$Codon
  for (i in 1:length(codons))
  {
    aa <- trna[i,"AA"]
    #wobbles <- unlist(lapply(wobble.neighbors(reverse.complement(codons[i])),reverse.complement))
    trna.neighbors <- get.neighbors(codons[i])
    #trna.neighbors <- setdiff(trna.neighbors,wobbles)

    for (j in trna.neighbors)
    {
      if(j %in% c("TAA","TAG","TGA")) next
      cognates <- unlist(strsplit(df[i,"Cognates"],","))
      print(cognates)
      reverse <- reverse.complement(j)
      if (j %in% trna.present$Codon && (!reverse %in% cognates))
      {
        
        if (trna[which(trna$Codon== j),"AA"] == aa)
        {
          tmp <- df[i,"Pseudo"]
          tmp.split <- unlist(strsplit(tmp,","))

          if (!reverse %in% tmp.split)
          {
            df[i,"Pseudo"] <- paste(df[i,"Pseudo"],reverse,sep=",")
          }
        
        } else {
          tmp <- df[i,"Near"]
          tmp.split <- unlist(strsplit(tmp,","))
          if (!reverse %in% tmp.split)
          {
            df[i,"Near"] <- paste(df[i,"Near"],reverse,sep=",")
          }
        }
      }
    }
  }
  write.table(df,"../scer_cognate_near_cognate.tsv",sep="\t",row.names=F,col.names=T,quote=F)
}



checkWobble <- function(codon,trna)
{
  reverse <- reverse.complement(codon)
  trna.split <- unlist(strsplit(trna,""))
  codon.split <- unlist(strsplit(codon,""))
  reverse.split <- unlist(strsplit(reverse,""))
  if(reverse.complement(codon) == trna)
  {
    return(1)
  } else{
    index <- which(reverse.split != trna.split)
    if(index == 1)
    {
      t.index <- 1
      c.index <- 3
    } else if (index == 2){
      t.index <- 2
      c.index <- 2
    } else if (index == 3){
      t.index <- 3
      c.index <- 1
    } else{
      print("Error")
    }
    if (trna.split[t.index] == "C" && codon.split[c.index] == "T")
    {
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "T"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "G"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "C"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "A"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "T"){
      return(canonical.wobble)
    }
  }
}

harmonic.mean <- function(x)
{
  y <- 1/x
  summation <- sum(y)
  harm.mean <- length(x)/summation
  return(harm.mean)

}


## Estimate elongation and missense error rates.
## Will need to provide a file with mappings for codons to cognate, pseudo-cognatesm, and near-cognates and a file with tRNA numbers
## Also will want to provide file for output
## TODO: make these a function parameter instead of hard-coding
runShahModel <- function()
{
  #a <- 10.992
  target.rate <- 12.5 ## aa/s 
  p.c <- 6.52*10^-1
  p.n <- 6.2*10^-4
  p.p <- p.n # pseudo-cognates, same AA recognized through non-standard wobble

  canonical.wobble <- 0.6 ## purine-purine or pyrmidine-pyrimidine
  non.canonical.wobble <- 0.64 ## purine-pyrimidine (GT/AC)


  #r_d <- 3.146 * 10^-3

  map <- read.table("../ecoli_cognate_near_cognate.tsv",sep="\t",header=T,stringsAsFactors=F)
  trna <- read.table("../ecoli_tRNA_current.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
  df <- data.frame(AA=map$AA,Codon=map$Codon,R_c=numeric(length=nrow(map)),R_n=numeric(length=nrow(map)),e_m=numeric(length=nrow(map)),e_n=numeric(length=nrow(map)))
  for (i in 1:nrow(map))
  {
    row <- map[i,]
    codon <- row$Codon
    cognates <- unlist(strsplit(row$Cognate,","))
    pseudo <- unlist(strsplit(row$Pseudo,","))
    near <- unlist(strsplit(row$Near,","))
    r_c <- 0
    r_n <- 0
    for (j in cognates)
    {
      r_c <- r_c + (a*trna[reverse.complement(j),"tRNA"] * p.c * checkWobble(codon,j))
    }
    for (j in pseudo)
    {
      r_c <- r_c + (a*trna[reverse.complement(j),"tRNA"] * p.p * checkWobble(codon,j))
    }
    for (j in near)
    {
      r_n <- r_n + (a*trna[reverse.complement(j),"tRNA"] * p.n * checkWobble(codon,j))
    }
    df[i,c("R_c","R_n")] <- c(r_c,r_n)
  }
  harm.mean<-harmonic.mean(df[,c("R_c")]+df[,c("R_n")])
  a <- target.rate/harm.mean
  df[,c("R_c","R_n")] <- a*df[,c("R_c","R_n")]
  r_d <- target.rate/4000
  df[,"e_m"] <- (df[,"R_n"])/(df[,"R_c"] + df[,"R_n"] + r_d)
  df[,"e_n"] <- r_d/(df[,"R_c"] + df[,"R_n"] + r_d)
  write.table(df,"ecoli_missense_rates.tsv",row.names=F,col.names=T,quote=F,sep="\t")

  # premal <- read.table("~/CUB_structure_ecoli/shah_gilchrist_2010_ecoli.tsv.csv",sep="\t",header=T,stringsAsFactors=F)
  # tmp<-merge(df,premal,by="Codon")
  # pdf("mine_vs_premal_ecoli.pdf")
  # plot(tmp$R_c,tmp$Rc,main="Cognate Elongation Rates",xlab="Alex",ylab="Premal")
  # #text(tmp$R_c,tmp$Rc,labels=tmp$Codon)
  # abline(a=0,b=1,lty=2,col="red")
  # plot(tmp$R_n,tmp$Rn,main="Near-cognate Elongation Rates",xlab="Alex",ylab="Premal")
  # #text(tmp$R_n,tmp$Rn,labels=tmp$Codon)
  # abline(a=0,b=1,lty=2,col="red")
  # plot(tmp$e_m,tmp$eM,main="Missense Error Rates",xlab="Alex",ylab="Premal")
  # #text(tmp$e_m,tmp$eM,labels=tmp$Codon)
  # abline(a=0,b=1,lty=2,col="red")
  # plot(tmp$e_n,tmp$eN,main="Nonsense Error Rates",xlab="Alex",ylab="Premal")
  # #text(tmp$e_n,tmp$eN,labels=tmp$Codon)
  # abline(a=0,b=1,lty=2,col="red")
  # dev.off()

}

compareEfficientAccurate <- function()
{
  df <- read.table("scer_missense_rates.tsv",sep="\t",header=T,stringsAsFactors=F)
  aa <- unique(df$AA)
  df[,"Rank_R_c"] <- 1
  df[,"Rank_e_m"] <- 1
  for (a in aa)
  {
    if(a == "M" || a == "W") next
    tmp <- df[which(df$AA==a),]
    tmp[,"Rank_R_c"] <- rank(-tmp$R_c)
    tmp[,"Rank_e_m"] <- rank(tmp$e_m)
    df[which(df$AA==a),] <- tmp
  }
  return(df)
}





distributionRatesAccuracy <- function(metric="e_m",output="output.pdf",sig=c("AGC","TCA","TCC","TTC","GCC","GGA","GTA","GGC"),abs.value=T,rescale=T,xlab="e_m",title="title",binwidth)
{
  df <- read.table("scer_missense_rates.tsv",sep="\t",header=T,stringsAsFactors=F)
  df <- df[which(!df$AA %in% c("M","W")),]
  df[,"Total_R"] <- df[,"R_c"] + df[,"R_n"]
  if(rescale)
  {
    df <- rescaleAGCT(df,metric=metric)
  }
  if(abs.value)
  {
    df[,metric] <- abs(df[,metric])
  }
  df["Significance"] <- "Not Significant"
  df[which(df$Codon %in% sig),"Significance"] <- "Significant"
  x <- df[which(df$Codon %in% sig),metric]
  y <- df[which(!df$Codon %in% sig),metric]
  print(wilcox.test(x,y,var.equal=T))
  p <- ggplot(df, aes_string(x=metric,fill="Significance",color="Significance")) + geom_histogram(alpha=0.5,binwidth=binwidth) + xlab(xlab) +ggtitle(title)
  ggsave(plot=p,file=output)
}


#focal.vs.neighbors()
#runShahModel()
df<- compareEfficientAccurate()

# distributionRatesAccuracy(metric="e_m",output="../Images_for_dissertation/scer_compare_abs_value_error_rates.pdf",abs.value=T,rescale=T,xlab="| Relative Error Rate |",title="Distribution of Relative Error Rates",,binwidth=0.001)
# distributionRatesAccuracy(metric="R_c",output="../Images_for_dissertation/scer_compare_abs_value_cognate_rates.pdf",abs.value=T,rescale=T,xlab="| Relative Elongation Rate |", title= "Distribution of Relative Elongation Rates",binwidth=2)

# distributionRatesAccuracy(metric="e_m",output="../Images_for_dissertation/scer_compare_error_rates.pdf",abs.value=F,rescale=T,xlab="Relative Error Rate",title="Distribution of Relative Error Rates",,binwidth=0.001)
# distributionRatesAccuracy(metric="R_c",output="../Images_for_dissertation/scer_compare_cognate_rates.pdf",abs.value=F,rescale=T,xlab="Relative Elongation Rate", title= "Distribution of Relative Elongation Rates",binwidth=2)


# df <- read.table("scer_missense_rates.tsv",sep="\t",header=T,stringsAsFactors=F)
# df[,"Total_R"] <- df[,"R_c"] + df[,"R_n"]
# df <- rescaleAGCT(df,metric=c("R_c","e_m"))
# sig <- c("AGC","TCA","TCC","TTC","GCC","GGA","GTA","GGC")
# pdf("rescaled_elong_err_rates.pdf")
# plot(df$R_c,df$e_m,main="Comparison of Rescaled Elongation and Error Rates",xlab="Elongation Rates",ylab="Missense Error Rates")
# text(df[which(df$Codon %in% sig),"R_c"],df[which(df$Codon %in% sig),"e_m"],labels=df[which(df$Codon %in% sig),"Codon"])
# dev.off()
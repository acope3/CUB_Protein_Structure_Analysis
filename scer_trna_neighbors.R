
p.c <- 6.52*10^-1
p.n <- 6.2*10^-4
p.p <- p.n # pseudo-cognates, same AA recognized through non-standard wobble

## Taking target rate from Kramer et. al. Nat. Struct. & Mol. Biol. 2009
target.aa.s <- 9

canonical.wobble <- 0.4 ## purine-purine or pyrmidine-pyrimidine
non.canonical.wobble <- 0.36 ## purine-pyrimidine (GT/AC)

purines <- c("G","A")
pyrimidines <- c("T","C")

get.neighbors <- function(codon)
{
  neighbors <- character(length=12)
  nuc <- unlist(strsplit(codon,split='',fixed=T))
  neighbors[1:4] <- paste0(nuc[1],c("A","G","T","C"),nuc[3])
  neighbors[5:8] <- paste0(c("A","G","T","C"),nuc[2],nuc[3])
  neighbors[9:12] <- paste0(nuc[1],nuc[2],c("A","G","T","C"))
  return(neighbors)
}


trna <- read.table("scer_tRNA.tsv",sep="\t",header=T,stringsAsFactors = F)
trna$tRNA <- ifelse(trna$tRNA%%1 == 0,trna$tRNA,0)
codons <- trna$Codon
trna[,"Neighbors.trna"] <- rep(0,nrow(trna))
for (i in codons)
{
  nuc <- unlist(strsplit(i,split='',fixed=T))
  neighbors <- get.neighbors(i)
  aa <- trna[which(trna$Codon == i),"AA"]
  codon.neighbors <- trna[which((trna$Codon %in% neighbors)),"Codon"])

  ## Initialize cognate rate with trna and p.c for the given codon
  cognate.rate <- trna[which(trna$Codon == i),"tRNA"] * p.c 
  near.cognate.rate <- 0
  for (j in codon.neighbors)
  {
    neighbor.nuc <- unlist(strsplit(j,split='',fixed=T))
    neighbor.aa <- aa <- trna[which(trna$Codon == j),"AA"]
    ## If cognate or pseudo-cognate, will only vary at the 3rd codon position
    if (aa == neighbor.aa)
    {
      if ((nuc[3] %in% purines && neighbor.nuc[3] %in% purines) || (nuc[3] %in% pyrimidines && neighbor.nuc[3] %in% pyrimidines))
      {

      }
    }
  }
}


#deg <- c(2,4)
#counts <- table(trna$AA)
#amino <- names(which(counts == deg))
# amino <- c("A","S","T","P")
# avoid <- c()
# trna.aa <- trna[which(trna$AA %in% amino & !(trna$Codon %in% avoid)),]
# plot(trna.aa$tRNA,trna.aa$Neighbors.trna)

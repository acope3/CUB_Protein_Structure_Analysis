library(AnaCoDa)

target.fasta <- "Ordered_disordered/Ordered/ordered.fasta"
phi.file <- "Ordered_disordered/Ordered/ordered_phi.csv"

genome <- initializeGenomeObject(target.fasta)
phi <- read.table(phi.file,sep=",",header=T,stringsAsFactors=F)


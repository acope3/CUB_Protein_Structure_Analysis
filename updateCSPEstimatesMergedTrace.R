library(AnaCoDa)

createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T,dir_name="./")
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref)
  }
}

input <- "Final_runs/Beta/Structure_Downstream_Combo/"
fasta.folders <- input #, "../data/cds/sampled/",  "../data/cds/sampled/", "../data/cds/filtered/")
fasta.files <- list.files(path=fasta.folders,pattern="*.fasta",full.names = F)
mixture.labels <- unlist(strsplit(fasta.files,split=".fasta"))
fasta.paths <- paste0(fasta.folders, fasta.files)
numMixtures <- length(fasta.files)

parameter <- loadParameterObject("Final_runs/Beta/Results/Structure_Downstream_Combo_with_fixed_mutation/combined_parameter.Rda")
createParameterOutput(parameter,numMixtures,15000,mixture.labels,dir_name = "Final_runs/Beta/Results/Structure_Downstream_Combo_with_fixed_mutation/")
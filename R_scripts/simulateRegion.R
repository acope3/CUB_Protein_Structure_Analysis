library(AnaCoDa)

simulate <- function(genome.file,sel.file,mut.file,phi.file,output.dir,output.fasta)
{
  dir.create(output.dir)
  genome <- initializeGenomeObject(genome.file)
  phi <- read.table(phi.file,sep=",",header=T)
  parameter <- initializeParameterObject(genome,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome)),mixture.definition = "allUnique",initial.expression.values =phi[,2],split.serine = TRUE)
  parameter$initSelectionCategories(sel.file,1,F)
  parameter$initMutationCategories(mut.file,1,F)
  model <- initializeModelObject(parameter,model="ROC",with.phi=FALSE)
  model$simulateGenome(genome)
  genome$writeFasta(paste(output.dir,output.fasta,sep="/"),simulated=T)
  file.copy(phi.file,output.dir)
  #
}
sel.file <- "../selection_mod_scerevisiae.csv"
mut.file <- "../mutation_mod_scerevisiae.csv"
genome.file <- "../Data/Fasta/complete_seq.fasta"
phi.file <- "../complete_seq_phi.csv"
output.dir <- "../Data/Fasta/"
output.fasta <- "simulated_scer_ROC.fasta"
simulate(genome.file,sel.file,mut.file,phi.file,output.dir,output.fasta)

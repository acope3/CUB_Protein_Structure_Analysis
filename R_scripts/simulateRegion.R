## Author: Alex Cope
## Code for simulating genomes. Currently setup to simulate 100 genomes for testing Pechmann and Frydman hypothesis. 

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
  ## Uncomment this line if you want to copy the phi file to the directory of the simulated genome.
  file.copy(phi.file,output.dir)
  #
}
sel.file <- "../selection_mod_scerevisiae.csv"
mut.file <- "../mutation_mod_scerevisiae.csv"
genome.file <- "../Scer/Predicted/Secondary_structures/Coil/coil.fasta"
phi.file <- "../Scer/Predicted/Secondary_structures/Coil/coil_phi.csv"
output.dir <- "../Scer/Simulated/Secondary_structures/Coil"
output.fasta <- "coil.fasta"
simulate(genome.file,sel.file,mut.file,phi.file,output.dir,output.fasta)

# genome.file <- "../Data/Fasta/complete_seq.fasta"
# phi.file <- "../Data/Fasta/complete_seq_phi.csv"
# output.dir <- "../Scer/Test_Fisher_Exact/"
# output.fasta <- "sim.fasta"
# for (i in 1:100)
# {
#   simulate(genome.file,sel.file,mut.file,phi.file,paste0(output.dir,i,"/"),output.fasta)
# }
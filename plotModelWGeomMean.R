library(AnaCoDa)

genome <- initializeGenomeObject("Data/Yeast/complete_seq.fasta")
dis.phi <- read.table("")
parameter <- initializeParameterObject(genome,1,1,rep(1,length(genome)), split.serine = TRUE, mixture.definition = "allUnique", initial.expression.values = init_phi)
model<- initializeModelObject(parameter,model="ROC")
pdf("with_genome_.pdf",height=12,width=10)
plotModel(model,genome,samples=1,mixture=1)
dev.off()
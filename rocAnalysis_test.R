library(AnaCoDa)
rm(list=ls())

args<-(commandArgs(TRUE))

if(length(args)==0)
{
  div <- 0 
  input <- "DeaneSaundersClass/"
  directory <- "Results/DS_classes/"
  thin <- 10
  adapt <- 50
  samp <- 500
  num_threads <- 8
}else{
  for(i in 1:length(args))
  {
    eval(parse(text=args[[i]]))
  }
}


createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref)
  }
}

createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)
    plot(model, genome, samples = samples*samples.percent.keep, mixture = i,main = mixture.labels[i])
  }
}

fasta.folders <- input #, "../data/cds/sampled/",  "../data/cds/sampled/", "../data/cds/filtered/")
fasta.files <- list.files(path=fasta.folders,pattern="*.fasta",full.names = F)
mixture.labels <- unlist(strsplit(fasta.files,split=".fasta"))
fasta.paths <- paste0(fasta.folders, fasta.files)
phi.files <- list.files(path=fasta.folders,pattern="*_phi.csv",full.names = F)
phi.path <- paste0(fasta.folders, phi.files)
numMixtures <- length(fasta.files)
mixture.sizes <- rep(0, numMixtures)

## Note: writing a for loop to deal with all mixtures (1 - n.mixtures) is tricky.
## Part of the issue is the appending of the object defined in the command and the assignment of the output
mixture.index <- 1;

genome <- initializeGenomeObject(file=fasta.paths[mixture.index],match.expression.by.id = FALSE,append = FALSE)

mixture.sizes[mixture.index] <- length(genome)
if(numMixtures > 1){
  for(mixture.index in 2:numMixtures)
  {
    tmp.length <- length(genome)
    genome <- initializeGenomeObject(file=fasta.paths[mixture.index],genome=genome,match.expression.by.id = FALSE,append = TRUE)
    mixture.sizes[mixture.index] <- length(genome) - tmp.length
  }
}

if(length(genome) != sum(mixture.sizes)){
  stop("length(genomeObj) != sum(mixture.sizes), but it should.")
}else{
  print("FASTA successfully files loaded:");
  print(fasta.files[1:numMixtures])
}

cat("Genome loaded\n")
#initialize parameter object

sphi_init <- rep(0.01,numMixtures)
with.phi <- F
mixDef <- "mutationShared"
percent.to.keep <- 0.5
size <- length(genome)
cat(size,"\n")
index <- c(1:size)
#geneAssignment <- c(rep(1,size.tmp),rep(2,size.tmp.2-size.tmp),rep(3,size-size.tmp.2))
geneAssignment <- rep(1:numMixtures, mixture.sizes)

init_phi <- c()
for (i in phi.path)
{
  segment_exp <- read.table(file=i,sep=",",header=TRUE)
  init_phi <- c(init_phi,segment_exp[,2])
}

if(length(genome) != length(init_phi)){
  stop("length(genomeObj) != length(init_phi), but it should.")
}else{
  print("Initial Phi values successfully files loaded:");
}

parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi)
parameter$initMutationCategories(c("mutation_mod_scerevisiae.csv"),1,TRUE)
parameter$initSelectionCategories(c("selection_mod_scerevisiae.csv"),numMixtures,FALSE)

# initialize MCMC object
samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=FALSE, est.csp=TRUE, est.hyper=FALSE,est.mix = FALSE)
mcmc$setStepsToAdapt((samples*thinning)/2)
# get model object
model <- initializeModelObject(parameter, "ROC", with.phi)

system.time(
  runMCMC(mcmc, genome, model, num_threads,divergence.iteration = div)
)

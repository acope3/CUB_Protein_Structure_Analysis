#!/usr/bin/env Rscript
library(AnaCoDa)
library(profmem)
library(argparse)
rm(list=ls())


parser <- ArgumentParser()
parser$add_argument("-i","--input",type="character",default="./")
parser$add_argument("-o","--output",type="character",default="./")
parser$add_argument("-d","--div",type="integer",default=0)
parser$add_argument("-s","--samp",type="integer",default=5000)
parser$add_argument("-a","--adapt",type="integer",default=50)
parser$add_argument("-t","--thin",type="integer",default=20)
parser$add_argument("-n","--threads",type="integer",default=1)
args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thin <- args$thin
adapt <- args$adapt
samp <- args$samp
num_threads <- args$threads

print(args)
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

genome <- initializeGenomeObject(file=fasta.paths[mixture.index],match.expression.by.id = FALSE,append = FALSE, observed.expression.file=phi.path[mixture.index])
print(fasta.files)
print(phi.files)
mixture.sizes[mixture.index] <- length(genome)
if(numMixtures > 1){
  for(mixture.index in 2:numMixtures)
    {
    tmp.length <- length(genome)
    genome <- initializeGenomeObject(file=fasta.paths[mixture.index],genome=genome,match.expression.by.id = FALSE,append = TRUE, observed.expression.file=phi.path[mixture.index])
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


with.phi <- T
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
sphi_init <- rep(sd(log(init_phi)),numMixtures)
print(sphi_init)
if(length(genome) != length(init_phi)){
  stop("length(genomeObj) != length(init_phi), but it should.")
}else{
  print("Initial Phi values successfully files loaded:");
}


parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment,init.sepsilon = 0.05,split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi)
#parameter<-initializeParameterObject(init.with.restart.file = paste("Final_runs/Beta/Results/Structure_Downstream_Combo/run_7/Restart_files/rstartFile.rst_final",sep="/"))

#parameter<-initializeParameterObject(init.with.restart.file = "Current_3/Results/Structure_Downstream_Combo_fixed/run_5/Restart_files/test.rst")
#parameter$fixDM()
#print(parameter$getMutationPriorMean())
#print(parameter$getMutationPriorStandardDeviation())
parameter$initMutationCategories(c("mutation_mod_scerevisiae.csv"),1,F)
#initialize MCMC object
samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=T, est.csp=TRUE, est.hyper=F,est.mix = FALSE)
mcmc$setStepsToAdapt((samples*thinning)/2)
# get model object
model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise=T)

run_number <- 1
dir.create(directory)
dir_name <- paste0(directory,"/run_",run_number)
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
#run mcmc on genome with parameter using model
  sys.runtime<-system.time(
  runMCMC(mcmc, genome, model, num_threads,divergence.iteration = div)
)
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)
  
createParameterOutput(parameter = parameter,numMixtures = numMixtures,mixture.labels = mixture.labels,samples = samples,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.5)
# expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*0.5)
# write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples*percent.to.keep,mixture.name=mixture.labels)
plot(parameter,what="Selection",samples=samples*percent.to.keep,mixture.name=mixture.labels)
dev.off()

trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
dev.off()


pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
#plot(trace, what = "ExpectedPhi")
aa <- aminoAcids()
done.adapt <- TRUE
for(a in aa)
{
  if (a=="M"||a=="X"||a=="W") next
  accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
  len <- length(accept.trace)
  mean.acceptance <- mean(accept.trace[(len-len*0.5):len])
  if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
  plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
}
acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
dev.off()

param.conv <- T
for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
}

writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.2)
z<-abs(diag$z)
done <- (z < 1.96) && param.conv
rm(parameter)
rm(trace)
rm(model)
while((!done) && (run_number <= 3))
{
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
  #parameter$fixDM()
  run_number <- run_number + 1
  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))

  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=T, est.csp=TRUE, est.hyper=F,est.mix=FALSE)
  if(!done.adapt)
  {
    mcmc$setStepsToAdapt((samples*thinning)/2)
    percent.to.keep <- 0.5
  } else{
    mcmc$setStepsToAdapt(0)
    percent.to.keep <- 1.0
  }
  model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise=T)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  #run mcmc on genome with parameter using model
  #p<-profmem({
  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=0)
  )
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)

  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)

  # mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*..to.keep)
  # expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep)
  # write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  #
  #
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
 # plot(trace, what = "ExpectedPhi")
  aa <- aminoAcids()
  done.adapt <- TRUE
  param.conv <- TRUE
  if (!done.adapt)
  {
    for(a in aa)
    {
      if (a=="M"||a=="X"||a=="W") next
      accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
      len <- length(accept.trace)
      mean.acceptance <- mean(accept.trace[(len-len*percent.to.keep):len])
      if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
      plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
    }
  }
  acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
  dev.off()
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.2)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
  }
  pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
  plot(parameter,what="Mutation",samples=samples*percent.to.keep,mixture.name=mixture.labels)
  plot(parameter,what="Selection",samples=samples*percent.to.keep,mixture.name=mixture.labels)
  dev.off()


  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = percent.to.keep,mixture.labels = mixture.labels)
  dev.off()
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  #})
  #p <- total(p)
  #write(p,file=paste(dir_name,"memory_allocation_total.csv",sep="/"),ncolumns=1)
  #write.table(as.data.frame(p),file =paste(dir_name,"memory_allocation.csv",sep="/"),col.names = T,row.names = F,quote = F )
  #rm(p)
  diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)
  done <- (z < 1.96) && param.conv
  rm(parameter)
  rm(trace)
  rm(model)
}

# if (!param.conv)
# {
#   samples <- samp + samp/2
# } else{
#   samples <- samp
# }
samples <- 10000
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
#parameter$fixDM()
run_number <- run_number + 1
dir_name <- paste0(directory,"/final_run")
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=T, est.csp=TRUE, est.hyper=F,est.mix=FALSE)

mcmc$setStepsToAdapt(0)

model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise=T)
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
#run mcmc on genome with parameter using model
#p<-profmem({
sys.runtime <- system.time(
  runMCMC(mcmc, genome, model, num_threads)
)
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)
  

createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = 1,relative.to.optimal.codon = F,report.original.ref = T)

# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.5)
# expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*0.5)
# write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


#plots different aspects of trace
trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
#plot(trace, what = "ExpectedPhi")
acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples)
dev.off()

pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples,mixture.name=mixture.labels)
plot(parameter,what="Selection",samples=samples,mixture.name=mixture.labels)
dev.off()


pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
dev.off()

for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples,thin = thinning,what="Selection",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
}
rm(trace)
rm(model)
writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))

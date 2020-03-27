#!/usr/bin/env Rscript


library(profmem)
library(argparse)
rm(list=ls())

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file containing protein-coding sequences",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-d","--div",help="Number of steps to diverge from starting values. Will be applied at beginning of each run, with the exception of the last.",type="integer",default=0)
parser$add_argument("-s","--samp",help="Number of samples",type="integer",default=5000)
parser$add_argument("-a","--adapt",help="Adaptive Width",type="integer",default=50)
parser$add_argument("-t","--thin",help="Thinning value",type="integer",default=20)
parser$add_argument("-n","--threads",help="Number of threads to use for MCMC",type="integer",default=1)
parser$add_argument("--phi",help="Initial Phi values. Assumes csv format with Gene IDs in first column and Phi values in second column", type="character",default="ROC_Scer_Est/yeast_phi.csv")
parser$add_argument("--a1",type="integer",default=8000,help="Initial value for initiation cost a_1")
parser$add_argument("--mut",type="character",default="ROC_Scer_Est/yeast_dM.csv",help="Initial value for dM parameters. Expects csv file. First column must be the amino acid. Second column must be the codon. Third column should be the value for dM.")
parser$add_argument("--sel",type="character",default="omega_10_5.csv",help="Initial value for dOmega parameters. Expects csv file. First column must be the amino acid. Second column must be the codon. Third column should be the value for dOmega.")
parser$add_argument("--est_csp",help="Use this flag to indicate estimation of Phi. Otherwise, Phi will not be estimated",action="store_true")
parser$add_argument("--est_phi",help="Use this flag to indicate estimation of CSP. Otherwise, CSP will not be estimated",action="store_true")
parser$add_argument("--est_hyp",help="Use this flag to indicate estimation of Hyperparameters. Otherwise, Hyperparameters will not be estimated",action="store_true")
parser$add_argument("--max_num_runs",help="Max number of runs to do. Note that each fitting will consists of at least two runs. The first can serve as a burn-in.",type="integer",default = 6)

library(AnaCoDa)

args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thin <- args$thin
adapt <- args$adapt
samp <- args$samp
a1 <- args$a1
sel <- args$sel
mut <- args$mut
num_threads <- args$threads
phi.files <- args$phi
est.phi <- args$est_phi
est.csp <- args$est_csp
est.hyp <- args$est_hyp
max_num_runs <- args$max_num_runs


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
numMixtures <- length(fasta.files)
mixture.sizes <- rep(0, numMixtures)

## Note: writing a for loop to deal with all mixtures (1 - n.mixtures) is tricky.
## Part of the issue is the appending of the object defined in the command and the assignment of the output
mixture.index <- 1;

genome <- initializeGenomeObject(file=fasta.paths[mixture.index],observed.expression.file="Data/yeast_arith_mean_gene_expression_time_course.csv",match.expression.by.id = FALSE,append = FALSE,positional = T)

mixture.sizes[mixture.index] <- length(genome)
if(numMixtures > 1){
  for(mixture.index in 2:numMixtures)
  {
    tmp.length <- length(genome)
    genome <- initializeGenomeObject(file=fasta.paths[mixture.index],genome=genome,match.expression.by.id = FALSE,append = TRUE,positional = T)
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
segment_exp <- read.table(file=phi.files,sep=",",header=TRUE)
init_phi <-segment_exp[,2]

sphi_init <- rep(1,numMixtures)
with.phi <- F
mixDef <- "allUnique"
percent.to.keep <- 0.5 ## percent of steps to adapt during first few runs
size <- length(genome)
cat(size,"\n")
index <- c(1:size)
geneAssignment <- rep(1:numMixtures, mixture.sizes)

segment_exp <- read.table(file=phi.files,sep=",",header=TRUE)
init_phi <-segment_exp[,2]

sphi_init <- rep(sd(init_phi),numMixtures)

if(length(genome) != length(init_phi)){
  stop("length(genomeObj) != length(init_phi), but it should.")
}else{
  print("Initial Phi values successfully files loaded:");
}

#parameter<-initializeParameterObject(init.with.restart.file = "Genome_fits/Yeast_Full_Genome_Fit_wphi/final_run/Restart_files/rstartFile.rst_final",model="ROC")

parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment,split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi,model="ROC",init.initiation.cost=a1,,init.csp.variance=0.0000000025)

parameter$initMutationCategories(c(mut),1,F)
parameter$initSelectionCategories(c(sel),1,F)
#initialize MCMC object

samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix = FALSE)
mcmc$setStepsToAdapt((samples*thinning)*percent.to.keep)
# get model object
model <- initializeModelObject(parameter,"ROC",with.phi=with.phi,fix.observation.noise=F)

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
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)

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
plot(trace, what = "ExpectedPhi")
plot(trace,what="Sphi")
plot(trace,what="Sepsilon")
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
acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
dev.off()


for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_omega_",i,".txt"),ncolumns = 1)
}


for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
}


writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.2)
z<-abs(diag$z)
done <- (z < 1.96) && param.conv
rm(parameter)
rm(trace)
rm(model)
while((!done) && (run_number < max_num_runs - 1))
{
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
  run_number <- run_number + 1
  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix=FALSE)

  model <- initializeModelObject(parameter, "ROC",with.phi=with.phi,fix.observation.noise=F)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=div)
  )
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)
  
  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
  
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)

  #
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  

  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  plot(trace,what="Sphi")
  plot(trace,what="Sepsilon")
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
  acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
  dev.off()
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_omega_",i,".txt"),ncolumns = 1)
  }
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
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
  
  diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)
  done <- (z < 1.96) && param.conv
  rm(parameter)
  rm(trace)
  rm(model)
}

samples <- samp
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
run_number <- run_number + 1
dir_name <- paste0(directory,"/final_run")
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix=FALSE)

mcmc$setStepsToAdapt(0)

model <- initializeModelObject(parameter, "ROC",with.phi=with.phi,fix.observation.noise=F)
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
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


#plots different aspects of trace
trace <- parameter$getTraceObject()

pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace,what="ExpectedPhi")
plot(trace,what="Sphi")
plot(trace,what="Sepsilon")
acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples)
acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples)
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

for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
}
rm(trace)
rm(model)
writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))



#!/usr/bin/env Rscript

library(profmem)
library(argparse)
rm(list=ls())

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-d","--div",help="Number of steps to diverge from starting values. Will be applied at beginning of each run, with the exception of the last.",type="integer",default=0)
parser$add_argument("-s","--samp",help="Number of samples",type="integer",default=5000)
parser$add_argument("-a","--adapt",help="Adaptive Width By Samples, i.e. will adapt every i samples",type="integer",default=50)
parser$add_argument("-t","--thin",help="Thinning value. Total number of iterations will be samples * thinning",type="integer",default=20)
parser$add_argument("-n","--threads",help="Number of threads to use for MCMC",type="integer",default=1)
parser$add_argument("--dEta",help="Initial dEta values. Assumes csv format with columns AA,Codon,DEta. First line should be a header.",type="character")
parser$add_argument("--dM",help="Initial dM values. Assumes csv format with columns AA,Codon,DM. First line should be a header.",type="character")
parser$add_argument("--est_csp",help="Use this flag to indicate estimation of Phi. Otherwise, Phi will not be estimated",action="store_true")
parser$add_argument("--est_phi",help="Use this flag to indicate estimation of CSP. Otherwise, CSP will not be estimated",action="store_true")
parser$add_argument("--est_hyp",help="Use this flag to indicate estimation of Hyperparameters. Otherwise, Hyperparameters will not be estimated",action="store_true")
parser$add_argument("--max_num_runs",help="Max number of runs to do.",type="integer",default = 6)
parser$add_argument("--fix_dEta",help="Use this flag to fix dEta at starting value.",action="store_true")
parser$add_argument("--fix_dM",help="Use this flag to fix dM at starting value",action="store_true")
parser$add_argument("--development",help="Run a developmental version of AnaCoDa",action="store_true")


args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thinning <- args$thin
adaptiveWidth <- args$adapt
samples <- args$samp
num_threads <- args$threads
dEta.file <- args$dEta
dM.file <- args$dM
est.phi <- args$est_phi
est.csp <- args$est_csp
est.hyp <- args$est_hyp
max_num_runs <- args$max_num_runs
fix_dEta <- args$fix_dEta
fix_dM <- args$fix_dM
dev <- args$development

if(dev)
{
	library(AnaCoDa,lib.loc="~/R_dev/")
} else{
	library(AnaCoDa)
}

print(args)


## Outputs CSP estimates
createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T,thin=thin)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref,thin=thin)
  }
}


## Outputs traces for CSPs and plots expected frequecies as function of log10(\phi) 
createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)
    plot(trace,what="AcceptanceRatio")
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
numMixtures <- 1


with.phi <- F ## Change this if you want to estimate phi using empirical phi values. Will need to specify file with empirical phi when initialzing Genome Object

print(fasta.paths[1])
print(phi.path[1])
if (with.phi)
{
  genome <- initializeGenomeObject(file=fasta.paths[1],match.expression.by.id=F,observed.expression.file=phi.path[1])
} else{
  #genome <- initializeGenomeObject(file=input,codon_table=12)
  genome <- initializeGenomeObject(file=fasta.paths[1])
}

mixDef <- "allUnique" ## options are allUnique, mutationShared, selectionShared, but this only matters if have multiple categories

percent.to.keep <- 0.5 ## Script is set up to stop adapting after 50% of samples and calculate posterior estimates, convergence, mean acceptance rate, etc. from last 50% of samples.
size <- length(genome)
index <- c(1:size)

geneAssignment <- rep(1,size)

init_phi <- NULL

segment_exp <- read.table(file=phi.path[1],sep=",",header=TRUE)
init_phi <- c(init_phi,segment_exp[,2])
if(length(genome) != length(init_phi))
{
  stop("length(genomeObj) != length(init_phi), but it should.")
} else{
print("Initial Phi values successfully files loaded:");
 }
##1.238 is posterior mean estimate for s_phi from full genome
sphi_init <- rep(1.238,numMixtures)
mixture.labels <- unlist(strsplit(input,"/"))
mixture.labels <- mixture.labels[length(mixture.labels)]



dir.create(directory)
done <- FALSE
done.adapt <- FALSE
run_number <- 1
# if (run_number != 1 && !est.hyp)
# {
#   parameter <- loadParameterObject(file.path(dir_name,"R_objects","parameter.Rda"))
#   trace <- parameter$getTraceObject()
#   sphi.trace <- trace$getStdDevSynthesisRateTraces()
#   sphi_init <- mean(unlist(sphi.trace))
# }
# print(sphi_init)
while((!done) && (run_number <= max_num_runs))
{

  if (run_number == 1)
  {
    dir_name <- paste0(directory,"/restart_",run_number)
    if (length(dM.file) > 0)
    {
      mutation <- read.table(dM.file,sep=",",header=T,stringsAsFactors=F)
      mutation.prior.mean <- mutation[,3]
      mutation.prior.sd <- mutation[,4]
    } else{
      mutation.prior.mean <- 0
      mutation.prior.sd <- 0.35
    }
    percent.to.keep <- 0.5
    parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment,init.sepsilon = 0.05,split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi)#,propose.by.prior=T,mutation.prior.mean=mutation.prior.mean,mutation.prior.sd=mutation.prior.sd)

    ## This assumes only one mixture category
    ## TODO: generalize to allow for more than one mixture
    if (length(dM.file) > 0)
    {
      parameter$initMutationCategories(dM.file,1,fix_dM)
    } 

    if (length(dEta.file) > 0)
    {
      parameter$initSelectionCategories(dEta.file,1,fix_dEta)
    }
    steps.to.adapt <- (samples*thinning)*percent.to.keep
  } else if (run_number == max_num_runs){
    percent.to.keep <- 1
    parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
    # parameter$setStdDevSynthesisRate(sphi_init[1],0)
    if (fix_dM)
    {
      parameter$fixDM()
    }
    if (fix_dEta)
    {
      parameter$fixDEta()
    }
    dir_name <- paste0(directory,"/final_run")
    steps.to.adapt <- 0
    div <- 0
  } else{
    percent.to.keep <- 0.5
    parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
    dir_name <- paste0(directory,"/restart_",run_number)
    if (fix_dM)
    {
      parameter$fixDM()
    }
    if (fix_dEta)
    {
      parameter$fixDEta()
    }
    steps.to.adapt <- (samples*thinning)*percent.to.keep
    div <- 0
    dir_name <- paste0(directory,"/restart_",run_number)
  }
  
  dir_name <- paste0(directory,"/restart_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))

  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix=FALSE)
 
  mcmc$setStepsToAdapt(steps.to.adapt)

  #mcmc$setStepsToAdapt(0)
  model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise=T)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  

  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=div)
  )

  ## Output runtime for mcmc
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)


  ## Creates R objects, which can be later loaded for re-analzying already completed runs
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
 

  ## Output CSP file
  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T,thin=thinning)

  ## Output phi file
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep,genome = genome)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  if (est.hyp)
  {
    plot(trace,what="Sphi")
  }

  if (est.csp)
  {

    ## Calculate auto-correlation and convergence of CSP traces
    param.conv <- TRUE
    if (!fix_dEta)
    {
      ##acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
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
    }
    if (!fix_dM)
    {
      ##acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
      for (i in 1:numMixtures)
      {
        param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.2)
        z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
        if (length(z.scores) > 0)
        {
          param.conv <- FALSE
        }
        write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_M_",i,".txt"),ncolumns = 1)
      }
    }
  }
  dev.off()
  


  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = percent.to.keep,mixture.labels = mixture.labels)
  dev.off()
  
 
  diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)

  ## Can end if overall log(posterior) and CSP parameters have converged
  done <- (z < 1.96) && param.conv
  done <- FALSE
  rm(parameter)
  rm(trace)
  rm(model)
  run_number <- run_number + 1
}

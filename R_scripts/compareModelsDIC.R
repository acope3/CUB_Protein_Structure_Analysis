## Author: Alexander Cope
## Code to calculate DIC scores for different models. 


library(AnaCoDa)


#' @details Calculates the loglikelihood for the given model fit at the provided values of dM (mutation), dEta (selection), and phi (protein production rates)
#'
#' @param genome a Genome object initialized by the 
#' @param model a model object created from a parameter object from a previous run
#' @param mutation a data.frame from getCSPEstimates for a ROC model fit representing mutation estimates
#' @param selection a data.frame from getCSPEstimates for a ROC model fit representing selection estimates
#' @param phi a vector of phi values (on the natural scale) produced from getExpressionEstimates (must be in same order as genome)
#' 
#' @return Loglikelihood at Posterior mean estimates
calculateLogLikAtPosteriorMean <- function(genome,model,mutation,selection,phi)
{
  aa <- aminoAcids()
  loglik <- 0
  for (a in aa)
  {
    if (a == "M" || a == "X" || a == "W") next
    codon.counts <- getCodonCountsForAA(a,genome)
    codons <- AAToCodon(a,F)
    ## Use numerica indexing to avoid issue of updating getCSPEstimates to use Mean instead of Posterior
    dEta <- selection[which(selection$AA == a),3]
    dM <- mutation[which(mutation$AA == a),3]
    codon.prob.per.gene <- data.frame(matrix(unlist(lapply(phi,function(x){
      log(model$CalculateProbabilitiesForCodons(dM, dEta, x))
    })),nrow=length(genome),ncol=length(codons),byrow=T))
    loglik <- loglik + sum(codon.prob.per.gene*codon.counts,na.rm = T)
  }
  return(loglik)
}


#' @details Calculates DIC, pDIC, and the loglikelihood at the posterior mean values for a model
#'
#' @param model.fits a vector of directories containing the runs. Note this assumes my directory structure for ROC fits, so will likely need to be changed
#' @param genome.path a vector of FASTA files corresponding to the genomes analyzed by the runs indicated by model.fits 
#' @param samples how many samples to keep for calculating DIC and pDIC
#'
#' 
#' @return list with DIC, pDIC, and likelihood at posterior mean estimates
calculate.for.all.fits<-function(model.fits,genome.path,samples=10000)
{
  subfits.likelihood <- numeric(length=length(model.fits))
  subfits.likelihood.trace <- matrix(rep(0,length(model.fits)*samples),nrow=length(model.fits),ncol=samples)
  fasta.files <- list.files(path=genome.path,pattern="*.fasta",full.names = T,recursive = F)
  fasta.paths <- fasta.files
  fasta.paths<-sort(fasta.paths)
  model.fits <- sort(model.fits)
  for (i in 1:length(model.fits))
  {
    parameter <- loadParameterObject(paste0(model.fits[i],"/restart_5/R_objects/parameter.Rda"))
    model <- initializeModelObject(parameter, "ROC", F)
    genome <- initializeGenomeObject(fasta.paths[i])
    size <- length(genome)
    mcmc <- loadMCMCObject(paste0(model.fits[i],"/restart_5/R_objects/mcmc.Rda"))
    csp <- getCSPEstimates(parameter,samples=samples,relative.to.optimal.codon = F,report.original.ref = F)
    dM <- csp$Mutation
    dEta <- csp$Selection
    phi <- getExpressionEstimates(parameter,c(1:size),samples)

    ## Note that this assumes the gene names were not included in the phi data.frame
    subfits.likelihood[i] <- calculateLogLikAtPosteriorMean(genome,model,mutation = dM,selection = dEta,phi=phi[,1])
    loglike.trace <- mcmc$getLogLikelihoodTrace()
    length.trace <- length(loglike.trace)
    subfits.likelihood.trace[i,] <- loglike.trace[(length.trace - samples+1):length.trace]
  }
  ## sum up loglikelihood value across runs at each sample to get a trace representing all runs
  ## will produce vector that has a length equal to samples
  combined.lik.trace <- colSums(subfits.likelihood.trace)

  ## get sum of loglikelihoods at parameter posterior means across all model fits
  lik.posterior.mean <- sum(subfits.likelihood)

  ## pDIC is the effective number of parameters. For my runs where I'm just estimating dEta, this typically is close to the actual number of parameters
  pdic <- 2 *(lik.posterior.mean - (1/samples)*sum(combined.lik.trace[1:samples]))
  DIC <- -2 * lik.posterior.mean + 2 * pdic
  return(list(DIC=DIC,pdic=pdic,loglik.at.mean.posterior=lik.posterior.mean))
}





within.structure.diff <- function(target.dir,head.dir="../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E")
{


  genome.files <- c(file.path(head.dir,target.dir,"Start_helix/"),
                  file.path(head.dir,target.dir,"Core_helix/"),
                  file.path(head.dir,target.dir,"End_helix/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_helix/"),
                  file.path(head.dir,"Results",target.dir,"Core_helix/"),
                  file.path(head.dir,"Results",target.dir,"End_helix/"))
  DIC.helix.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
  

  genome.files <- c(file.path(head.dir,target.dir,"Start_helix_End_helix/"),
                  file.path(head.dir,target.dir,"Core_helix/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_helix_End_helix/"),
                  file.path(head.dir,"Results",target.dir,"Core_helix/"))
  DIC.helix.termini.tog <- calculate.for.all.fits(model.files,genome.files,samples=10000)

  genome.files <- c(file.path(head.dir,target.dir,"Helix/"))
  model.files <- c(file.path(head.dir,"Results",target.dir,"Helix/"))
  DIC.helix <- calculate.for.all.fits(model.files,genome.files,samples=10000)


  helix.df <- data.frame(Structure=rep("Helix",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.helix$DIC,DIC.helix.termini.tog$DIC,DIC.helix.termini$DIC),
              Delta.DIC = c(DIC.helix$DIC - DIC.helix$DIC,DIC.helix$DIC - DIC.helix.termini.tog$DIC,DIC.helix$DIC - DIC.helix.termini$DIC))

  genome.files <- c(file.path(head.dir,target.dir,"Start_coil/"),
                  file.path(head.dir,target.dir,"Core_coil/"),
                  file.path(head.dir,target.dir,"End_coil/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_coil/"),
                  file.path(head.dir,"Results",target.dir,"Core_coil/"),
                  file.path(head.dir,"Results",target.dir,"End_coil/"))
  DIC.coil.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
  
  genome.files <- c(file.path(head.dir,target.dir,"Start_coil_End_coil/"),
                  file.path(head.dir,target.dir,"Core_coil/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_coil_End_coil/"),
                  file.path(head.dir,"Results",target.dir,"Core_coil/"))
  DIC.coil.termini.tog <- calculate.for.all.fits(model.files,genome.files,samples=10000)

  genome.files <- c(file.path(head.dir,target.dir,"Turn_Coil/"))
  model.files <- c(file.path(head.dir,"Results",target.dir,"Turn_Coil/"))
  DIC.coil <- calculate.for.all.fits(model.files,genome.files,samples=10000)

  coil.df <- data.frame(Structure=rep("Coil",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.coil$DIC,DIC.coil.termini.tog$DIC,DIC.coil.termini$DIC),
              Delta.DIC = c(DIC.coil$DIC - DIC.coil$DIC,DIC.coil$DIC - DIC.coil.termini.tog$DIC,DIC.coil$DIC - DIC.coil.termini$DIC))



  genome.files <- c(file.path(head.dir,target.dir,"Start_sheet/"),
                  file.path(head.dir,target.dir,"Core_sheet/"),
                  file.path(head.dir,target.dir,"End_sheet/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_sheet/"),
                  file.path(head.dir,"Results",target.dir,"Core_sheet/"),
                  file.path(head.dir,"Results",target.dir,"End_sheet/"))
  DIC.sheet.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
  
  genome.files <- c(file.path(head.dir,target.dir,"Start_sheet_End_sheet/"),
                  file.path(head.dir,target.dir,"Core_sheet/"))
 
  model.files <- c(file.path(head.dir,"Results",target.dir,"Start_sheet_End_sheet/"),
                  file.path(head.dir,"Results",target.dir,"Core_sheet/"))
  DIC.sheet.termini.tog <- calculate.for.all.fits(model.files,genome.files,samples=10000)

  genome.files <- c(file.path(head.dir,target.dir,"Sheet/"))
  model.files <- c(file.path(head.dir,"Results",target.dir,"Sheet/"))
  DIC.sheet <- calculate.for.all.fits(model.files,genome.files,samples=10000)

  sheet.df <- data.frame(Structure=rep("Sheet",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.sheet$DIC,DIC.sheet.termini.tog$DIC,DIC.sheet.termini$DIC),
              Delta.DIC = c(DIC.sheet$DIC - DIC.sheet$DIC,DIC.sheet$DIC - DIC.sheet.termini.tog$DIC,DIC.sheet$DIC - DIC.sheet.termini$DIC))

  return(list(Helix=helix.df,
              Coil=coil.df,
              Sheet=sheet.df))

}





genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
        "Scer/Predicted/Secondary_structures/Coil/",
        "Scer/Predicted/Secondary_structures/Sheet/",
        "Scer/Predicted/Secondary_structures/Nterminus/")
model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
        "Scer/Predicted/Results/Secondary_structures/Coil/",
        "Scer/Predicted/Results/Secondary_structures/Sheet/",
        "Scer/Predicted/Results/Secondary_structures/Nterminus/")
DIC.results.ss <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print("Done with Secondary Structures")

model.fits <- c("Secondary_structure_paired_any_nucleotide/")
genome.locs <- c("Secondary_structure_paired_any_nucleotide")
genome.locs <- paste0("Scer/Predicted/",genome.locs,sep="/")
results.dir <- paste("Scer/Predicted/","Results/",sep="/")
model.locs <- paste0(results.dir,model.fits)
categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
model.files <- paste0(model.locs,categories)
genome.files <- paste0(genome.locs,categories)
DIC.results.ss_all_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Paired/Unpaired",DIC.results.ss_all_dis$DIC))



# genome.files <- c("Scer/Predicted/Secondary_structures/Helix_Coil/",
#         "Scer/Predicted/Secondary_structures/Sheet/",
#         "Scer/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix_Coil/",
#         "Scer/Predicted/Results/Secondary_structures/Sheet/",
#         "Scer/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.hc <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Coil merged")

# genome.files <- c("Scer/Predicted/Secondary_structures/Helix_Sheet/",
#         "Scer/Predicted/Secondary_structures/Coil/",
#         "Scer/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix_Sheet/",
#         "Scer/Predicted/Results/Secondary_structures/Coil/",
#         "Scer/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.he <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Sheet merged")

# genome.files <- c("Scer/Predicted/Secondary_structures/Sheet_Coil/",
#         "Scer/Predicted/Secondary_structures/Helix/",
#         "Scer/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Sheet_Coil/",
#         "Scer/Predicted/Results/Secondary_structures/Helix/",
#         "Scer/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.ec <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Coil and Sheet merged")


# # model.fits <- c("Nterm/")
# # model.type <- c("Nterm")
# # genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# # results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# # model.locs <- paste0(results.dir,model.fits)
# # categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# # model.files <- paste0(model.locs,categories)
# # genome.files <- paste0(genome.locs,categories)
# # DIC.results.comp <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# # print("Done with Full Genome (i.e. one category)")

# model.fits <- c("Ordered_disordered/")
# model.type <- c("Effects of disordered regions")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)

# DIC.results.ord <- calculate.for.all.fits(model.files,genome.files,samples=10000)


# genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
# 				"Scer/Predicted/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_dis/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
# 				"Scer/Predicted/Results/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_dis/")
# DIC.results.ss_coil_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil Ordered/Disordered",DIC.results.ss_coil_dis$DIC))

# genome.files <- c("Scer/Predicted/Secondary_structure_order/Helix_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Helix_dis/",
# 				"Scer/Predicted/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structures/Coil/")
# model.files <- c("Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Helix_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Helix_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structures/Coil/")
# DIC.results.ss_helix_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Helix Ordered/Disordered",DIC.results.ss_helix_dis$DIC))

# genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_dis/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structures/Coil/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Sheet_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Sheet_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structures/Coil/")
# DIC.results.ss_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Sheet Ordered/Disordered",DIC.results.ss_sheet_dis$DIC))


# genome.files <- c("Scer/Predicted/Secondary_structure_order/Helix_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Helix_dis/",
# 				"Scer/Predicted/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_dis/")
# model.files <- c("Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Helix_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Helix_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_dis/")
# DIC.results.ss_coil_dis_helix_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil/Helix Ordered/Disordered",DIC.results.ss_coil_dis_helix_dis$DIC))

# genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_dis/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_dis/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Sheet_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Sheet_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_dis/")
# DIC.results.ss_coil_dis_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil/Sheet Ordered/Disordered",DIC.results.ss_coil_dis_sheet_dis$DIC))


# model.fits <- c("Secondary_structure_order_02_05_2021/")
# genome.locs <- c("Secondary_structure_order")
# genome.locs <- paste0("Scer/Predicted/",genome.locs,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)
# DIC.results.ss_all_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + All Ordered/Disordered",DIC.results.ss_all_dis$DIC))


# genome.files <- c("Scer/Predicted/Secondary_structure_order/Helix_ord/",
#         "Scer/Predicted/Secondary_structure_order/Sheet_ord/",
#         "Scer/Predicted/Secondary_structure_order/Coil_ord/",
#         "Scer/Predicted/Secondary_structures/Nterminus/",
#         "Scer/Predicted/Ordered_disordered/Disordered/")
# model.files <- c("Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Helix_ord/",
#         "Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Sheet_ord/",
#         "Scer/Predicted/Results/Secondary_structure_order_02_05_2021/Coil_ord/",
#         "Scer/Predicted/Results/Secondary_structures/Nterminus/",
#         "Scer/Predicted/Results/Ordered_disordered/Disordered/")
# DIC.results.ss_vs_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures or Disordered",DIC.results.ss_vs_dis$DIC))


# struct.4.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini")
# struct.5.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_5_2_codon_for_termini")
# struct.6.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_6_2_codon_for_termini")
# struct.6.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_6")
# struct.7.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_7_2_codon_for_termini")
# struct.7.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_7")
# struct.8.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_8_2_codon_for_termini")
# struct.8.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_8")
# struct.9.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_9")
# struct.10.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_10")



# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Start_helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Remainder_helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Start_helix/",
#                  "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Remainder_helix/")
# DIC.pechmann <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_pechmann_exclude_less_than_6_exclude_GI/Helix/")
# DIC.helix <- calculate.for.all.fits(model.files,genome.files,samples=10000)




# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_helix_types/Helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_helix_types/GI_Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_helix_types/Helix/",
#                  "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_helix_types/GI_Helix/")
# DIC.helix.gi <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_helix_types/Helix_GI_Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_helix_types/Helix_GI_Helix/")
# DIC.helix <- calculate.for.all.fits(model.files,genome.files,samples=10000)



# # genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Start_helix/",
# #                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Remainder_helix/")
# # model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/Start_helix/",
# #                  "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/Remainder_helix/")
# # DIC.pechmann <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/GI_Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/Helix/",
#                  "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/GI_Helix/")
# DIC.helix.gi <- calculate.for.all.fits(model.files,genome.files,samples=10000)


# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_length_4_2_codon_for_termini/Helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_5_2_codon_for_termini/Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_begin_end_length_4_2_codon_for_termini/Helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_begin_end_exclude_less_than_5_2_codon_for_termini/Helix/")
# DIC.helix.4.vs.rest <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Helix/")
# DIC.helix.min.4 <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Core_helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Start_helix_End_helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Secondary_structures_begin_end_6_min_2_termini/Helix")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/Core_helix/",
#                  "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_terminii/Start_helix_End_helix/",
#                   "../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_6_min_2_termini/Helix/")
# DIC.helix.4.or.5 <- calculate.for.all.fits(model.files,genome.files,samples=10000)

# genome.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Helix/")
# model.files <- c("../Scer/Exp_conservative_homology_remove_X_G_I_as_H_B_as_E/Test/Results/Secondary_structures_begin_end_exclude_less_than_4_2_codon_for_termini/Helix/")
# DIC.helix.min.4 <- calculate.for.all.fits(model.files,genome.files,samples=10000)
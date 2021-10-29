## Author: Alexander Cope
## Code to calculate DIC scores for different models. 


library(AnaCoDa)

calculateUnscaledLogLik <- function(genome,model,mutation,selection,phi)
{
  aa <- aminoAcids()
  loglik <- 0
  for (a in aa)
  {
    if (a == "M" || a == "X" || a == "W") next
    codon.counts <- getCodonCountsForAA(a,genome)
    codons <- AAToCodon(a,F)
    dEta <- selection[which(selection$AA == a),3]
    dM <- mutation[which(mutation$AA == a),3]
    codon.prob.per.gene <- data.frame(matrix(unlist(lapply(phi,function(x){
      log(model$CalculateProbabilitiesForCodons(dM, dEta, x))
    })),nrow=length(genome),ncol=length(codons),byrow=T))
    loglik <- loglik + sum(codon.prob.per.gene*codon.counts,na.rm = T)
  }
  return(loglik)
}



calculate.for.all.fits<-function(model.fits,genome.phi,samples=10000)
{
  subfits.likelihood <- numeric(length=length(model.fits))
  subfits.likelihood.trace <- matrix(rep(0,length(model.fits)*samples),nrow=length(model.fits),ncol=samples)
  fasta.files <- list.files(path=genome.phi,pattern="*.fasta",full.names = T,recursive = F)
  fasta.paths <- fasta.files
  fasta.paths<-sort(fasta.paths)
  model.fits <- sort(model.fits)
  for (i in 1:length(model.fits))
  {
  	if (file.exists(paste0(model.fits[i],"/restart_5/R_objects/parameter.Rda")))
  	{
  		parameter <- loadParameterObject(paste0(model.fits[i],"/restart_5/R_objects/parameter.Rda"))
    	mcmc <- loadMCMCObject(paste0(model.fits[i],"/restart_5/R_objects/mcmc.Rda"))
  	} else{
  		parameter <- loadParameterObject(paste0(model.fits[i],"/final_run/R_objects/parameter.Rda"))
    	mcmc <- loadMCMCObject(paste0(model.fits[i],"/final_run/R_objects/mcmc.Rda"))
  	}
    
    model <- initializeModelObject(parameter, "ROC", F)
    genome <- initializeGenomeObject(fasta.paths[i])
    size <- length(genome)
    
    csp <- getCSPEstimates(parameter,samples=samples,relative.to.optimal.codon = F,report.original.ref = F)
    dM <- csp$Mutation
    dEta <- csp$Selection
    phi <- getExpressionEstimates(parameter,c(1:size),samples)
    subfits.likelihood[i] <- calculateUnscaledLogLik(genome,model,mutation = dM,selection = dEta,phi=phi[,1])
    loglike.trace <- mcmc$getLogLikelihoodTrace()
    length.trace <- length(loglike.trace)
    subfits.likelihood.trace[i,] <- loglike.trace[(length.trace - samples+1):length.trace]
  }
  combined.lik.trace <- colSums(subfits.likelihood.trace)
  mean.lik.posterior.mean <- sum(subfits.likelihood)
  pdic <- 2 *(mean.lik.posterior.mean - (1/samples)*sum(combined.lik.trace[1:samples]))
  DIC <- -2 * mean.lik.posterior.mean + 2 * pdic
  return(list(DIC=DIC,pdic=pdic,loglik.at.mean.posterior=mean.lik.posterior.mean))
}





within.structure.diff <- function(target.dir,head.dir="../Ecoli/Empirical")
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

  min.scores <- min(c(DIC.helix$DIC,DIC.helix.termini.tog$DIC,DIC.helix.termini$DIC))

  helix.df <- data.frame(Structure=rep("Helix",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.helix$DIC,DIC.helix.termini.tog$DIC,DIC.helix.termini$DIC),
              Delta.DIC = c(DIC.helix$DIC - min.scores, DIC.helix.termini.tog$DIC - min.scores, DIC.helix.termini$DIC - min.scores))

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

  min.scores <- min(c(DIC.coil$DIC,DIC.coil.termini.tog$DIC,DIC.coil.termini$DIC))

  coil.df <- data.frame(Structure=rep("Coil",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.coil$DIC,DIC.coil.termini.tog$DIC,DIC.coil.termini$DIC),
              Delta.DIC = c(DIC.coil$DIC - min.scores, DIC.coil.termini.tog$DIC - min.scores, DIC.coil.termini$DIC - min.scores))



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

  min.scores <- min(c(DIC.sheet$DIC,DIC.sheet.termini.tog$DIC,DIC.sheet.termini$DIC))
  
  sheet.df <- data.frame(Structure=rep("Sheet",3),
              Models=c("Whole Structure","Core Termini","Core N-terminus C-terminus"),
              DIC=c(DIC.sheet$DIC,DIC.sheet.termini.tog$DIC,DIC.sheet.termini$DIC),
              Delta.DIC = c(DIC.sheet$DIC - min.scores, DIC.sheet.termini.tog$DIC - min.scores, DIC.sheet.termini$DIC - min.scores))

  return(list(Helix=helix.df,
              Coil=coil.df,
              Sheet=sheet.df))

}



genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix/",
        "Ecoli/Predicted/Secondary_structures/Coil/",
        "Ecoli/Predicted/Secondary_structures/Sheet/",
        "Ecoli/Predicted/Secondary_structures/Nterminus/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix/",
        "Ecoli/Predicted/Results/Secondary_structures/Coil/",
        "Ecoli/Predicted/Results/Secondary_structures/Sheet/",
        "Ecoli/Predicted/Results/Secondary_structures/Nterminus/")
DIC.results.ss <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures",DIC.results.ss$DIC))

# genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix_Coil/",
#         "Ecoli/Predicted/Secondary_structures/Sheet/",
#         "Ecoli/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix_Coil/",
#         "Ecoli/Predicted/Results/Secondary_structures/Sheet/",
#         "Ecoli/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.hc <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Coil merged")

# genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix_Sheet/",
#         "Ecoli/Predicted/Secondary_structures/Coil/",
#         "Ecoli/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix_Sheet/",
#         "Ecoli/Predicted/Results/Secondary_structures/Coil/",
#         "Ecoli/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.he <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Sheet merged")

# genome.files <- c("Ecoli/Predicted/Secondary_structures/Sheet_Coil/",
#         "Ecoli/Predicted/Secondary_structures/Helix/",
#         "Ecoli/Predicted/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Sheet_Coil/",
#         "Ecoli/Predicted/Results/Secondary_structures/Helix/",
#         "Ecoli/Predicted/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.ec <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Coil and Sheet merged")


genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix_Sheet_Coil/",
        "Ecoli/Predicted/Secondary_structures/Nterminus/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix_Sheet_Coil/",
        "Ecoli/Predicted/Results/Secondary_structures/Nterminus/")
DIC.results.all <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures, All merged",DIC.results.all$DIC))


genome.files <- c("Ecoli/Predicted/Ordered_disordered/Ordered/",
        "Ecoli/Predicted/Ordered_disordered/Disordered/",
        "Ecoli/Predicted/Ordered_disordered/Nterminus/")
model.files <- c("Ecoli/Predicted/Results/Ordered_disordered/Ordered/",
        "Ecoli/Predicted/Results/Ordered_disordered/Disordered/",
        "Ecoli/Predicted/Results/Ordered_disordered/Nterminus/")

DIC.results.ord.2 <- calculate.for.all.fits(model.files,genome.files,samples=10000)


# genome.files <- c("Ecoli/Empirical/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Secondary_structures/Coil/",
#         "Ecoli/Empirical/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Secondary_structures/Turn/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Results/Secondary_structures/Coil/",
#         "Ecoli/Empirical/Results/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Results/Secondary_structures/Turn/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures")

# genome.files <- c("Ecoli/Empirical/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Secondary_structures/Turn_Coil/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Results/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Results/Secondary_structures/Turn_Coil/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.tc <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Turn and Coil Merged")


# genome.files <- c("Ecoli/Empirical/Secondary_structures/Helix_Coil_Turn/",
#         "Ecoli/Empirical/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Helix_Coil_Turn/",
#         "Ecoli/Empirical/Results/Secondary_structures/Sheet/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.hc <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Coil merged")

# genome.files <- c("Ecoli/Empirical/Secondary_structures/Helix_Sheet/",
#         "Ecoli/Empirical/Secondary_structures/Turn_Coil/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Helix_Sheet/",
#         "Ecoli/Empirical/Results/Secondary_structures/Turn_Coil/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.he <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Helix and Sheet merged")

# genome.files <- c("Ecoli/Empirical/Secondary_structures/Sheet_Coil_Turn/",
#         "Ecoli/Empirical/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Sheet_Coil_Turn/",
#         "Ecoli/Empirical/Results/Secondary_structures/Helix/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.ec <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, Coil and Sheet merged")

# genome.files <- c("Ecoli/Empirical/Secondary_structures/Helix_Coil_Sheet_Turn/",
#         "Ecoli/Empirical/Secondary_structures/Nterminus/")
# model.files <- c("Ecoli/Empirical/Results/Secondary_structures/Helix_Coil_Sheet_Turn/",
#         "Ecoli/Empirical/Results/Secondary_structures/Nterminus/")
# DIC.results.ss.hec <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures, null model")

model.fits <- c("Ordered_disordered/")
model.type <- c("Effects of disordered regions")
genome.locs <- paste0("Ecoli/Predicted/",model.fits,sep="/")
results.dir <- paste("Ecoli/Predicted/","Results/",sep="/")
model.locs <- paste0(results.dir,model.fits)
categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
model.files <- paste0(model.locs,categories)
genome.files <- paste0(genome.locs,categories)
DIC.results.ord <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Ordered vs. Disordered",DIC.results.ord$DIC))

genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix/",
				"Ecoli/Predicted/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_dis/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix/",
				"Ecoli/Predicted/Results/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_dis/")
DIC.results.ss_coil_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Coil Ordered/Disordered",DIC.results.ss_coil_dis$DIC))

genome.files <- c("Ecoli/Predicted/Secondary_structure_order/Helix_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Helix_dis/",
				"Ecoli/Predicted/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Secondary_structures/Coil/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structure_order/Helix_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Helix_dis/",
				"Ecoli/Predicted/Results/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Results/Secondary_structures/Coil/")
DIC.results.ss_helix_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Helix Ordered/Disordered",DIC.results.ss_helix_dis$DIC))

genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix/",
				"Ecoli/Predicted/Secondary_structure_order/Sheet_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Sheet_dis/",
				"Ecoli/Predicted/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Secondary_structures/Coil/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Sheet_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Sheet_dis/",
				"Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Results/Secondary_structures/Coil/")
DIC.results.ss_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Sheet Ordered/Disordered",DIC.results.ss_sheet_dis$DIC))


genome.files <- c("Ecoli/Predicted/Secondary_structure_order/Helix_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Helix_dis/",
				"Ecoli/Predicted/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_dis/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structure_order/Helix_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Helix_dis/",
				"Ecoli/Predicted/Results/Secondary_structures/Sheet/",
				"Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_dis/")
DIC.results.ss_coil_dis_helix_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Coil/Helix Ordered/Disordered",DIC.results.ss_coil_dis_helix_dis$DIC))

genome.files <- c("Ecoli/Predicted/Secondary_structures/Helix/",
				"Ecoli/Predicted/Secondary_structure_order/Sheet_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Sheet_dis/",
				"Ecoli/Predicted/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Secondary_structure_order/Coil_dis/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structures/Helix/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Sheet_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Sheet_dis/",
				"Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_ord/",
				"Ecoli/Predicted/Results/Secondary_structure_order/Coil_dis/")
DIC.results.ss_coil_dis_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + Coil/Sheet Ordered/Disordered",DIC.results.ss_coil_dis_sheet_dis$DIC))


model.fits <- c("Secondary_structure_order/")
model.type <- c("Secondary structure_order")
genome.locs <- paste0("Ecoli/Predicted/",model.fits,sep="/")
results.dir <- paste("Ecoli/Predicted/","Results/",sep="/")
model.locs <- paste0(results.dir,model.fits)
categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
model.files <- paste0(model.locs,categories)
genome.files <- paste0(genome.locs,categories)
DIC.results.ss_all_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures + All Ordered/Disordered",DIC.results.ss_all_dis$DIC))


genome.files <- c("Ecoli/Predicted/Secondary_structure_order/Helix_ord/",
        "Ecoli/Predicted/Secondary_structure_order/Sheet_ord/",
        "Ecoli/Predicted/Secondary_structure_order/Coil_ord/",
        "Ecoli/Predicted/Secondary_structures/Nterminus/",
        "Ecoli/Predicted/Ordered_disordered/Disordered/")
model.files <- c("Ecoli/Predicted/Results/Secondary_structure_order/Helix_ord/",
        "Ecoli/Predicted/Results/Secondary_structure_order/Sheet_ord/",
        "Ecoli/Predicted/Results/Secondary_structure_order/Coil_ord/",
        "Ecoli/Predicted/Results/Secondary_structures/Nterminus/",
        "Ecoli/Predicted/Results/Ordered_disordered/Disordered/")
DIC.results.ss_vs_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures or Disordered",DIC.results.ss_vs_dis$DIC))


#struct.4.2 <- within.structure.diff("Secondary_structures_begin_end_length_at_least_4_2_codon_for_termini",head.dir="Ecoli/Empirical")
#struct.5.2 <- within.structure.diff("Secondary_structures_begin_end_length_at_least_5_2_codon_for_termini",head.dir="Ecoli/Empirical")
#struct.6.2 <- within.structure.diff("Secondary_structures_begin_end_length_at_least_6_2_codon_for_termini",head.dir="Ecoli/Empirical")
# struct.6.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_6")
#struct.7.2 <- within.structure.diff("Secondary_structures_begin_end_length_at_least_7_2_codon_for_termini",head.dir="Ecoli/Empirical")
# struct.7.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_7")
# struct.8.2 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_8_2_codon_for_termini")
# struct.8.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_8")
# struct.9.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_9")
# struct.10.3 <- within.structure.diff("Secondary_structures_begin_end_exclude_less_than_10")

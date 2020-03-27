library(AnaCoDa)
library(parallel)

calculateUnscaledLogLik <- function(genome,model,mutation,selection,phi)
{
  aa <- aminoAcids()
  loglik <- 0
  for (a in aa)
  {
    if (a == "M" || a == "X" || a == "W") next
    codon.counts <- getCodonCountsForAA(a,genome)
    codons <- AAToCodon(a,F)
    dEta <- selection[which(selection$AA == a),"Posterior"]
    dM <- mutation[which(mutation$AA == a),"Posterior"]
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
    parameter <- loadParameterObject(paste0(model.fits[i],"/final_run/R_objects/parameter.Rda"))
    model <- initializeModelObject(parameter, "ROC", F)
    genome <- initializeGenomeObject(fasta.paths[i])
    size <- length(genome)
    mcmc <- loadMCMCObject(paste0(model.fits[i],"/final_run/R_objects/mcmc.Rda"))
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


translateCategoriesFromFileNames <- function(filename)
{
  if (is.na(filename) || filename == "NA")
  {
    category <- NA
  }
  else if (filename == "Nterminus")
  {
    category <- "Nterminus"
  }
  else
  {
    categories <- unlist(strsplit(filename,split = "_",fixed = T))
    exit.site <- unique(categories[c(T,F)])
    a.site <- unique(categories[c(F,T)])
    category.exit <- paste0(exit.site,collapse=",")
    category.asite <- paste0(a.site,collapse=",")
    category.exit <- paste0("(",category.exit,")",collapse="")
    category.asite <- paste0("(",category.asite,")",collapse="")
    category <- paste(category.exit,category.asite,sep="/")
  }
  return(category)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

runAsiteExit_emp <- function(models,model.loc,genome.loc,output)
{
  models <- read.table(models,sep=",",header=F,stringsAsFactors = F,fill = T,na.strings = "")
  Loglik <- numeric(length=nrow(models))
  DIC.scores <- numeric(length=nrow(models))
  pd.scores <- numeric(length=nrow(models))
  num.param <- numeric(length=nrow(models))
  for (i in 1:nrow(models))
  {
    model.fits <- unlist(models[i,])
    model.fits <- model.fits[which(is.na(model.fits) == F)]
    model.files <- paste0(model.loc,model.fits)
    genome.files <- paste0(genome.loc,model.fits)
    results <- calculate.for.all.fits(model.files,genome.files,samples=10000)
    Loglik[i] <- results$loglik.at.mean.posterior
    DIC.scores[i] <- results$DIC
    pd.scores[i] <- results$pdic
  }
  models.DIC <- cbind(models,DIC.scores)
  models.DIC <- cbind(models.DIC,pd.scores)
  models.DIC <- cbind(models.DIC,Loglik)
  model.type <- c(rep("A-site",14),rep("Exit Tunnel",14))
  test <- models.DIC
  for (i in 1:nrow(models.DIC))
  {
    for (j in 1:(ncol(models.DIC) - 3))
    {
      test[i,j] <- translateCategoriesFromFileNames(models.DIC[i,j])
    }
    num.param[i] <- length(which(!is.na(test[i,1:(ncol(models.DIC)-3)]))) * 40
  }
  test[,"Num.Parameters"] <- num.param
  a.site.min <- which.min(test[(1:14),"DIC.scores"])
  d.dic.within.asite <- test[a.site.min,"DIC.scores"] - test[(1:14),"DIC.scores"]
  exit.site.min <- which.min(test[(15:28),"DIC.scores"])
  d.dic.within.exitsite <- test[exit.site.min+14,"DIC.scores"] - test[(15:28),"DIC.scores"]
  d.dic <- c(d.dic.within.asite,d.dic.within.exitsite)
  test[,"Delta.DIC.Within"] <- d.dic
  
  overall.best <- which.min(test[,"DIC.scores"])
  d.dic.overall <- test[overall.best,"DIC.scores"] - test[,"DIC.scores"]
  test[,"Delta.DIC.Overall"] <- d.dic.overall
  
  
  test[,"Model.type"] <- model.type
  
  test <- test[order(test$Delta.DIC.Overall,decreasing = T),]
  test <- round_df(test,digits=0)
  colnames(test) <- c(paste0("Group_",1:5,"_(Exit/A)"),"DIC","p.D","LogLik.at.posterior.mean","Num.Parameter","Within.Delta.DIC","Overall.Delta.DIC","Model.type")
  test <- test[,c(paste0("Group_",1:5,"_(Exit/A)"),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")]
  write.table(test,output,quote = F,sep="\t",row.names = F,col.names = T)
  
}

runAsite_emp <- function(models,model.loc,genome.loc,output)
{
  models <- read.table(models,sep=",",header=F,stringsAsFactors = F,fill = T,na.strings = "")
  Loglik <- numeric(length=nrow(models))
  DIC.scores <- numeric(length=nrow(models))
  pd.scores <- numeric(length=nrow(models))
  num.param <- numeric(length=nrow(models))
  for (i in 1:nrow(models))
  {
    model.fits <- unlist(models[i,])
    model.fits <- model.fits[which(is.na(model.fits) == F)]
    model.files <- paste0(model.loc,model.fits)
    genome.files <- paste0(genome.loc,model.fits)
    results <- calculate.for.all.fits(model.files,genome.files,samples=10000)
    Loglik[i] <- results$loglik.at.mean.posterior
    DIC.scores[i] <- results$DIC
    pd.scores[i] <- results$pdic
  }
  models.DIC <- cbind(models,DIC.scores)
  models.DIC <- cbind(models.DIC,pd.scores)
  models.DIC <- cbind(models.DIC,Loglik)
  model.type <- c(rep("Secondary Structure",14))
  test <- models.DIC
  for (i in 1:nrow(models.DIC))
  {
    num.param[i] <- length(which(!is.na(test[i,1:(ncol(models.DIC)-3)]))) * 40
  }
  test[,"Num.Parameters"] <- num.param
  a.site.min <- which.min(test[(1:14),"DIC.scores"])
  d.dic.within.asite <- test[a.site.min,"DIC.scores"] - test[(1:14),"DIC.scores"]
  d.dic <- c(d.dic.within.asite)
  test[,"Delta.DIC.Within"] <- d.dic
  
  overall.best <- which.min(test[,"DIC.scores"])
  d.dic.overall <- test[overall.best,"DIC.scores"] - test[,"DIC.scores"]
  test[,"Delta.DIC.Overall"] <- d.dic.overall
  
  
  test[,"Model.type"] <- model.type
  
  test <- test[order(test$Delta.DIC.Overall,decreasing = T),]
  test <- round_df(test,digits=2)
  colnames(test) <- c(paste0("Group_",1:5),"DIC","p.D","LogLik.at.posterior.mean","Num.Parameter","Within.Delta.DIC","Overall.Delta.DIC","Model.type")
  test <- test[,c(paste0("Group_",1:5),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")]
  write.table(test,output,quote = F,sep="\t",row.names = F,col.names = T)
  
}


runAsiteExit_pred <- function(models,model.loc,genome.loc,output)
{
  models <- read.table(models,sep=",",header=F,stringsAsFactors = F,fill = T,na.strings = "")
  Loglik <- numeric(length=nrow(models))
  DIC.scores <- numeric(length=nrow(models))
  pd.scores <- numeric(length=nrow(models))
  num.param <- numeric(length=nrow(models))
  for (i in 1:nrow(models))
  {
    model.fits <- unlist(models[i,])
    model.fits <- model.fits[which(is.na(model.fits) == F)]
    model.files <- paste0(model.loc,model.fits)
    genome.files <- paste0(genome.loc,model.fits)
    results <- calculate.for.all.fits(model.files,genome.files,samples=10000)
    Loglik[i] <- results$loglik.at.mean.posterior
    DIC.scores[i] <- results$DIC
    pd.scores[i] <- results$pdic
  }
  models.DIC <- cbind(models,DIC.scores)
  models.DIC <- cbind(models.DIC,pd.scores)
  models.DIC <- cbind(models.DIC,Loglik)
  model.type <- c(rep("A-site",4),rep("Exit Tunnel",4))
  test <- models.DIC
  for (i in 1:nrow(models.DIC))
  {
    for (j in 1:(ncol(models.DIC) - 3))
    {
      test[i,j] <- translateCategoriesFromFileNames(models.DIC[i,j])
    }
    num.param[i] <- length(which(!is.na(test[i,1:(ncol(models.DIC)-3)]))) * 40
  }
  test[,"Num.Parameters"] <- num.param
  a.site.min <- which.min(test[(1:4),"DIC.scores"])
  d.dic.within.asite <- test[a.site.min,"DIC.scores"] - test[(1:4),"DIC.scores"]
  exit.site.min <- which.min(test[(5:8),"DIC.scores"])
  d.dic.within.exitsite <- test[exit.site.min+4,"DIC.scores"] - test[(5:8),"DIC.scores"]
  d.dic <- c(d.dic.within.asite,d.dic.within.exitsite)
  test[,"Delta.DIC.Within"] <- d.dic
  
  overall.best <- which.min(test[,"DIC.scores"])
  d.dic.overall <- test[overall.best,"DIC.scores"] - test[,"DIC.scores"]
  test[,"Delta.DIC.Overall"] <- d.dic.overall
  
  
  test[,"Model.type"] <- model.type
  
  test <- test[order(test$Delta.DIC.Overall,decreasing = T),]
  test <- round_df(test,digits=2)
  colnames(test) <- c(paste0("Group_",1:4,"_(Exit/A)"),"DIC","p.D","LogLik.at.posterior.mean","Num.Parameter","Within.Delta.DIC","Overall.Delta.DIC","Model.type")
  test <- test[,c(paste0("Group_",1:4,"_(Exit/A)"),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")]
  write.table(test,output,quote = F,sep="\t",row.names = F,col.names = T)
  
}


appendModelFit <- function(model.fits,DIC.results,model.type="Transitions",current.models="fixed_dic_scores_cutoff_issue_resolved_combined.tsv",output="updated_models.tsv",translate=T)
{
  df <- read.table(current.models,sep="\t",header=T,stringsAsFactors = F)
  total <- nrow(df)
  num.models <- length(model.fits)
  current.cat <- length(grep("Group_",colnames(df),fixed=T))
  if (current.cat < num.models)
  {
    for (i in current.cat + 1:(length(model.fits) - current.cat))
    {
      df[,paste0("Group_",i,"_.Exit.A.")] <- character(length=total)
    }
  }
  else if (current.cat > num.models)
  {
    for (i in 1:(current.cat-num.models))
    {
      model.fits <- c(model.fits,NA)
    }
  }
  if (translate)
  {
    cleaned.model.fits <- character(length(model.fits))
    for (i in 1:length(model.fits))
    {
      cleaned.model.fits[i] <- translateCategoriesFromFileNames(model.fits[i])
    }
    model.fits <- cleaned.model.fits
  }
  test <- df[,c(paste0("Group_",1:max(current.cat,num.models),"_.Exit.A."),"Model.type","Num.Parameter","p.D","DIC","LogLik.at.posterior.mean","Within.Delta.DIC")]
  num.param <- 40 * num.models
  
  tmp <- c(model.fits,model.type,num.param,DIC.results$pdic,DIC.results$DIC,DIC.results$loglik.at.mean.posterior,NA,-2*DIC.results$loglik.at.mean.posterior+2*num.param,NA)
  test <- rbind(test,tmp)
  test[c("Num.Parameter","p.D","DIC","Within.Delta.DIC","LogLik.at.posterior.mean")] <- lapply(test[c("Num.Parameter","p.D","DIC","Within.Delta.DIC","LogLik.at.posterior.mean")],as.numeric)
  
  overall.best <- which.min(test[,"DIC"])
  d.dic.overall <- test[overall.best,"DIC"] - test[,"DIC"]
  test[,"Overall.Delta.DIC"] <- d.dic.overall
  
  test <- test[order(test$Overall.Delta.DIC,decreasing = T),]
  test <- round_df(test,digits=0)
  test <- test[,c(paste0("Group_",1:max(current.cat,num.models),"_.Exit.A."),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")]
  colnames(test) <- c(paste0("Group_",1:max(current.cat,num.models),"_(Exit/A)"),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")
  write.table(test,output,quote = F,sep="\t",row.names = F,col.names = T)
}

appendModelFit2 <- function(model.fits,DIC.results,model.type="Transitions",current.models="fixed_dic_scores_cutoff_issue_resolved_combined.tsv",output="updated_models.tsv",translate=T)
{
  df <- read.table(current.models,sep="\t",header=T,stringsAsFactors = F)
  total <- nrow(df)
  num.models <- length(model.fits)
  current.cat <- length(grep("Group_",colnames(df),fixed=T))
  if (current.cat < num.models)
  {
    for (i in current.cat + 1:(length(model.fits) - current.cat))
    {
      df[,paste0("Group_",i)] <- character(length=total)
    }
  }
  else if (current.cat > num.models)
  {
    for (i in 1:(current.cat-num.models))
    {
      model.fits <- c(model.fits,NA)
    }
  }
  if (translate)
  {
    cleaned.model.fits <- character(length(model.fits))
    for (i in 1:length(model.fits))
    {
      cleaned.model.fits[i] <- translateCategoriesFromFileNames(model.fits[i])
    }
    model.fits <- cleaned.model.fits
  }
  test <- df[,c(paste0("Group_",1:max(current.cat,num.models)),"Model.type","Num.Parameter","p.D","DIC","LogLik.at.posterior.mean","Within.Delta.DIC")]
  num.param <- 40 * num.models
  
  tmp <- c(model.fits,model.type,num.param,DIC.results$pdic,DIC.results$DIC,DIC.results$loglik.at.mean.posterior,NA,-2*DIC.results$loglik.at.mean.posterior+2*num.param,NA)
  test <- rbind(test,tmp)
  test[c("Num.Parameter","p.D","DIC","Within.Delta.DIC","LogLik.at.posterior.mean")] <- lapply(test[c("Num.Parameter","p.D","DIC","Within.Delta.DIC","LogLik.at.posterior.mean")],as.numeric)
  
  overall.best <- which.min(test[,"DIC"])
  d.dic.overall <- test[overall.best,"DIC"] - test[,"DIC"]
  test[,"Overall.Delta.DIC"] <- d.dic.overall
  
  test <- test[order(test$Overall.Delta.DIC,decreasing = T),]
  test <- round_df(test,digits=2)
  test <- test[,c(paste0("Group_",1:max(current.cat,num.models)),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")]
  colnames(test) <- c(paste0("Group_",1:max(current.cat,num.models)),"Model.type","Num.Parameter","p.D","LogLik.at.posterior.mean","DIC","Within.Delta.DIC","Overall.Delta.DIC")
  write.table(test,output,quote = F,sep="\t",row.names = F,col.names = T)
}


getDIC <- function(model.loc,genome.loc,model.type,model.fits=NA)
{
  if (is.na(model.fits))
  {
    model.fits <- list.dirs(model.loc,full.names = F,recursive = F)
  }
  model.files <- paste0(model.loc,model.fits)
  genome.files <- paste0(genome.loc,model.fits)
  
  DIC <- calculate.for.all.fits(model.files,genome.files,samples=10000)
  row <- data.frame(Model.type=c(model.type),
                      Num.parameter=c(length(model.fits)*40),
                      p.D=DIC$pdic,
                      Loglik.at.posterior.mean=DIC$loglik.at.mean.posterior,
                      DIC=DIC$DIC,
                      stringsAsFactors = F)
  return(row)
}




appendDICToTable <- function(model.loc,genome.loc,categories,model.type,file.to.append)
{
  model.files <- paste0(model.loc,categories)
  genome.files <- paste0(genome.loc,categories)
  DIC.results <- calculate.for.all.fits(model.files,genome.files,samples=10000)
  if (is.na(DIC.results$DIC))
  {
    stop()
  }
  appendModelFit2(categories,DIC.results,current.models=file.to.append,model.type =model.type,translate = F,output = file.to.append)
}

createTable <- function(directory,output)
{
  results.dir <- paste(directory,"Results/",sep="/")
  
  model.fits <- c("Nterm/",
                  "Secondary_structures_pechmann/",
                  "Secondary_structures_begin_end/",
                  "Downstream_sheet/",
                  "Downstream_helix/",
                  "Downstream_structured/")
  model.type <- c("Nterminus vs. Remainder of Gene",
                  "Pechmann Hypothesis",
                  "Start vs. Core. vs. End",
                  "Difference if Downstream from Sheet",
                  "Difference if Downstream from Helix",
                  "Difference if Downstream from Structures")
  genome.locs <- paste0(directory,model.fits,sep="/")
  model.locs <- paste0(results.dir,model.fits)
  
  for (i in 1:length(model.fits))
  {
    categories <- list.dirs(model.locs[i],full.names = F,recursive = F)
    appendDICToTable(model.locs[i],genome.locs[i],categories,model.type[i],output)
  }
}

#runAsite_emp("models.txt","Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/","Scer/Exp_conservative_homology_missing_data/Secondary_structures/","dic_scer_exp_conservative_12_03.tsv")
#createTable("Scer/Exp_conservative_homology_missing_data/","dic_scer_exp_conservative_12_03.tsv")
#runAsiteExit_pred("runs_to_compare_predicted_chaperone_scer.txt",'Scer/Predicted/Results/Structure_Downstream_Combo/',"Scer/Predicted/Structure_Downstream_Combo/","dic_scer_pred.csv")
# model.fits <- c("Secondary_structures/")
# model.type <- c("Secondary structures")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)
# DIC.results.ss <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Secondary Structures")

# model.fits <- c("Nterm/")
# model.type <- c("Nterm")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)
# DIC.results.comp <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print("Done with Complete")

# model.fits <- c("Ordered_disordered/")
# model.type <- c("Effects of disordered regions")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)

# DIC.results.ord <- calculate.for.all.fits(model.files,genome.files,samples=10000)


# model.fits <- c("Downstream_ordered_disordered/")
# model.type <- c("Downstream_ordered_disordered")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)
# DIC.results.dod <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Downstream Ordered/Disordered",DIC.results.dod$DIC))



# genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
# 				"Scer/Predicted/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_dis/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
# 				"Scer/Predicted/Results/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_dis/")
# DIC.results.ss_coil_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil Ordered/Disordered",DIC.results.ss_coil_dis$DIC))

# genome.files <- c("Scer/Predicted/Secondary_structure_order/Helix_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Helix_dis/",
# 				"Scer/Predicted/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structures/Coil/")
# model.files <- c("Scer/Predicted/Results/Secondary_structure_order/Helix_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Helix_dis/",
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
# 				"Scer/Predicted/Results/Secondary_structure_order/Sheet_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Sheet_dis/",
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
# model.files <- c("Scer/Predicted/Results/Secondary_structure_order/Helix_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Helix_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Sheet/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_dis/")
# DIC.results.ss_coil_dis_helix_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil/Helix Ordered/Disordered",DIC.results.ss_coil_dis_helix_dis$DIC))

# genome.files <- c("Scer/Predicted/Secondary_structures/Helix/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Sheet_dis/",
# 				"Scer/Predicted/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Secondary_structure_order/Coil_dis/")
# model.files <- c("Scer/Predicted/Results/Secondary_structures/Helix/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Sheet_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Sheet_dis/",
# 				"Scer/Predicted/Results/Secondary_structures/Nterminus/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_ord/",
# 				"Scer/Predicted/Results/Secondary_structure_order/Coil_dis/")
# DIC.results.ss_coil_dis_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + Coil/Sheet Ordered/Disordered",DIC.results.ss_coil_dis_sheet_dis$DIC))


# model.fits <- c("Secondary_structure_order/")
# model.type <- c("Secondary structure_order")
# genome.locs <- paste0("Scer/Predicted/",model.fits,sep="/")
# results.dir <- paste("Scer/Predicted/","Results/",sep="/")
# model.locs <- paste0(results.dir,model.fits)
# categories <- list.dirs(model.locs[1],full.names = F,recursive = F)
# model.files <- paste0(model.locs,categories)
# genome.files <- paste0(genome.locs,categories)
# DIC.results.ss <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# print(paste("Done with Secondary Structures + All Ordered/Disordered",DIC.results.ss$DIC))


genome.files <- c("Scer/Predicted/Secondary_structure_order/Helix_ord/",
        "Scer/Predicted/Secondary_structure_order/Sheet_ord/",
        "Scer/Predicted/Secondary_structure_order/Coil_ord/",
        "Scer/Predicted/Secondary_structures/Nterminus/",
        "Scer/Predicted/Ordered_disordered/Disordered/")
model.files <- c("Scer/Predicted/Results/Secondary_structure_order/Helix_ord/",
        "Scer/Predicted/Results/Secondary_structure_order/Sheet_ord/",
        "Scer/Predicted/Results/Secondary_structure_order/Coil_ord/",
        "Scer/Predicted/Results/Secondary_structures/Nterminus/",
        "Scer/Predicted/Results/Ordered_disordered/Disordered/")
DIC.results.ss_coil_dis_sheet_dis <- calculate.for.all.fits(model.files,genome.files,samples=10000)
print(paste("Done with Secondary Structures or Disordered",DIC.results.ss_coil_dis_sheet_dis$DIC))



# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Start_helix/",
#                  "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Helix/",
#                  "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/End_helix/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Start_helix/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Helix/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/End_helix/")
# DIC.helix.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures/Helix/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/Helix/")
# DIC.helix <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Start_coil/",
#                   "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Turn_Coil/",
#                   "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/End_coil/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Start_coil/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Turn_Coil/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/End_coil/")
# DIC.coil.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures/Turn_Coil/")
# model.files <-  c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/Turn_Coil/")
# DIC.coil <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Start_sheet/",
#                   "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/Sheet/",
#                   "Scer/Exp_conservative_homology_missing_data/Secondary_structures_begin_end/End_sheet/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Start_sheet/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/Sheet/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_begin_end/End_sheet/")
# DIC.sheet.termini <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures/Sheet/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures/Sheet/")
# DIC.sheet <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 
# 
# genome.files <- c("Scer/Exp_conservative_homology_missing_data/Secondary_structures_pechmann/Start_helix/",
#                   "Scer/Exp_conservative_homology_missing_data/Secondary_structures_pechmann/Helix/")
# model.files <- c("Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_pechmann/Start_helix/",
#                  "Scer/Exp_conservative_homology_missing_data/Results/Secondary_structures_pechmann/Helix/")
# DIC.pechmann <- calculate.for.all.fits(model.files,genome.files,samples=10000)
# 

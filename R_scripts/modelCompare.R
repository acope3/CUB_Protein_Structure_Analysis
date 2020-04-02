library(AnaCoDa)
library(mvtnorm)

calculate_marginal_log_likelihood <- function(parameter, mcmc, mixture, n.samples, divisor,warnings=TRUE)
{  
  if(divisor < 1) stop("Generalized Harmonic Mean Estimation of Marginal Likelihood requires importance sampling distribution variance divisor be greater than 1")
  
  ## Collect information from AnaCoDa objects
  trace <- parameter$getTraceObject()
  ## This should be the posterior instead of the log_posterior but this causes an overflow, find fix!!!
  log_posterior <- mcmc$getLogPosteriorTrace()
  log_posterior <- log_posterior[(length(log_posterior) - n.samples+1):(length(log_posterior))]
  
  
  ### HANDLE CODON SPECIFIC PARAMETERS
  log_imp_dens_sample <- rep(0, n.samples)
  for (k in 1:mixture)
  {
  for(ptype in 0:1) # for all parameter types (mutation/selection parameters)
  {
    for(aa in AnaCoDa::aminoAcids()) # for all amino acids
    {
      if(aa == "M" || aa == "W" || aa == "X") next # ignore amino acids with only one codon or stop codons 
      codons <- AnaCoDa::AAToCodon(aa, focal = T)
      ## get covariance matrix and mean of importance distribution
      sample_mat <- matrix(NA, ncol = length(codons), nrow = n.samples)
      mean_vals <- rep(NA, length(codons))
      for(i in 1:length(codons)) # for all codons
      {
        vec <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(k, codons[i], ptype, TRUE)
        vec <- vec[(length(vec) - n.samples+1):(length(vec))]
        sample_mat[,i] <- vec
        mean_vals[i] <- mean(vec)
      }
      ## scale/shrinked covariance matrix
      cov_mat <- cov(sample_mat) / divisor
      if(all(cov_mat == 0))
      {
        if(warnings) print(paste("Covariance matrix for codons in amino acid",aa,"has 0 for all values. Skipping."))
        next
      }
      for(i in 1:n.samples)
      {
        ## calculate importance density for collected samples
        ## mikeg: We should double check with Russ that it is okay to use the full covariance matrix of the sample (even though we have no prior on the cov structure) when constructing the importance density function.
        ## It seems logical to do so
        log_imp_dens_aa <- dmvnorm(x = sample_mat[i,], mean = mean_vals, sigma = cov_mat, log = TRUE)
        log_imp_dens_sample[i] = log_imp_dens_sample[i] + log_imp_dens_aa
      }
    }
  }
  }
  ## HANDLE GENE SPECIFIC PARAMETERS
  
  # phi values are stored on natural scale.
  
  synt_trace <- trace$getSynthesisRateTrace()[[mixture]]
  n_genes = length(synt_trace);
  
  sd_vals <- rep(NA, n_genes)
  mean_vals <- rep(NA, n_genes)
  for(i in 1:n_genes) ## i is indexing across genes
  {
    vec <- synt_trace[[i]]
    vec <- vec[(length(vec) - n.samples+1):(length(vec))]
    sd_vals[i] <- sd(vec)
    if (all(sd_vals[i] == 0))
    {
      if(warnings) print(paste("Variance of gene",i,"is 0. Skipping."))
      next
    }
    mean_vals[i] <- mean(vec)
    log_mean_vals = log(mean_vals) - 0.5 * log(1+(sd_vals^2/mean_vals^2))
    log_sd_vals = sqrt(log(1+(sd_vals^2/mean_vals^2)))
    
    ## Calculate vector of importance density for entire \phi trace of gene.
    log_imp_dens_phi <- dlnorm(x = vec, meanlog = log_mean_vals[i], sdlog = log_sd_vals[i]/divisor, log = TRUE)
    
    
    ## update importance density function vector of sample with current gene;
    log_imp_dens_sample = log_imp_dens_sample + log_imp_dens_phi
    
  } ## end synth_trace loop
  
  ## Scale importance density for each sample by its posterior probability (on log scale)
  log_imp_dens_over_posterior = log_imp_dens_sample - log_posterior
  ## now scale by max term to facilitate summation
  max_log_term = max(log_imp_dens_over_posterior)
  ## Y = X - max_X
  offset_log_imp_dens_over_posterior <- log_imp_dens_over_posterior - max_log_term
  ## Z = sum(exp(vec(Y)))
  offset_sum_imp_dens_over_posterior <- sum(exp(offset_log_imp_dens_over_posterior))
  log_sum_imp_dens_over_posterior <- log(offset_sum_imp_dens_over_posterior) + max_log_term
  ## ln(ML) = ln(n) -(Z + max_X)
  log_marg_lik <- log(n.samples) - log_sum_imp_dens_over_posterior
  ##marg_lik = 1.0/(log_inv_marg_lik/n.samples) # equation 9
  return(log_marg_lik)
}



calculate.for.all.fits<-function(model.fits,turn.off.warnings = F)
{
  subfits.marginal.likelihood <- numeric(length=length(model.fits))
  for (i in 1:length(model.fits))
  {
    parameter <- loadParameterObject(paste0(model.fits[i],"/run_6/R_objects/parameter.Rda"))
    numMixtures <- length(unique(parameter$getMixtureAssignment()))
    mcmc <- loadMCMCObject(paste0(model.fits[i],"/run_6/R_objects/mcmc.Rda"))
    subfit <- calculate_marginal_log_likelihood(parameter,mcmc,numMixtures,n.samples=5000,divisor=1.0,warnings = F)
    subfits.marginal.likelihood[i] <- subfit
  }
  return(sum(subfits.marginal.likelihood))
}

# args <- commandArgs(trailingOnly = T)
# if (length(args) == 2)
# {
#   model.1.loc <- args[1]
#   model.2.loc <- args[2]
# } else{
#   stop("Please provide two directory locations for the necessary model comparisons, with the more complex model listed first")
# }

model.1.loc <- "../Simulations_BF/Results/Secondary_structures/"
#model.2.loc <- "Results/Downstream_30_5/"
#model.3.loc <- "Results/10_classes/"


model.fits.1 <- list.dirs(path=model.1.loc,full.names = T,recursive = F)
# model.fits.2 <- list.dirs(path=model.2.loc,full.names = T,recursive = F)
# model.fits.3 <- list.dirs(path=model.3.loc,full.names = T,recursive = F)

marginal.likelihood.complete <- calculate.for.all.fits("Final_runs/Beta/Results/Complete_genome/")
#marginal.likelihood.ds <- calculate.for.all.fits(model.fits.2)
#marginal.likelihood.10.classes <- calculate.for.all.fits(model.fits.3)
marginal.likelihood.ss <- calculate.for.all.fits("Final_runs/Beta/Results/Secondary_structures/")
bayes.factor.complete.ss <- marginal.likelihood.ss - marginal.likelihood.complete
# parameter <- loadParameterObject("Results/Secondary_structure_predicted_all_mixture/run_6/R_objects/parameter.Rda")
# mcmc <- loadMCMCObject("Results/Secondary_structure_predicted_all_mixture/run_6/R_objects/mcmc.Rda")
# # coil <- calculate_marginal_log_likelihood(parameter,mcmc,mixture = 1,n.samples = 15000,divisor=1.5,warnings=F)
# # helix <- calculate_marginal_log_likelihood(parameter,mcmc,mixture = 2,n.samples = 15000,divisor=1.5,warnings=F)
# # sheet <- calculate_marginal_log_likelihood(parameter,mcmc,mixture = 3,n.samples = 15000,divisor=1.5,warnings=F)
# mixture.run <- calculate_marginal_log_likelihood(parameter,mcmc,mixture = 3,n.samples = 15000,divisor=1.0,warnings=F)
# bs <- (coil + helix + sheet) - marginal.likelihood.ss
# print(bs)
#bayes.factor.ss.ds <- marginal.likelihood.ds - marginal.likelihood.ss
#bayes.factor.ds.10 <- marginal.likelihood.10.classes - marginal.likelihood.ds

# print("Domains/Linkers vs Whole-Genome")
# print(bayes.factor.complete.ss)
# print("Downstream + Secondary Structure vs Just Secondary Structure")
# print(bayes.factor.ss.ds)
# print("10 classes vs Downstream + Secondary Structure")
# print(bayes.factor.ds.10)







plotSinglePanel <- function(parameter, model, genome, expressionValues, samples, mixture, aa,codon.probability)
{
  codons <- AAToCodon(aa, T)
  # 
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T, log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T, log_scale = F)
  }
  # 
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  # codonProbability <- lapply(10^phis,  
  #                            function(phi){
  #                              model$CalculateProbabilitiesForCodons(mutation, selection, phi)
  #                            })
  codonProbability <- codon.probability
  #get codon counts
  codons <- AAToCodon(aa, F)
  codonCounts <- vector("list", length(codons))
  for(i in 1:length(codons))
  {
    codonCounts[[i]] <- genome$getCodonCountsPerGene(codons[i])
  }
  codonCounts <- do.call("cbind", codonCounts)
  # codon proportions
  codonCounts <- codonCounts / rowSums(codonCounts)
  codonCounts[is.nan(codonCounts)] <- NA # necessary if AA does not appear in gene
  
  # make empty plot
  xlimit <- range(expressionValues, na.rm = T)
  plot(NULL, NULL, xlim=xlimit, ylim=c(-0.05,1.05), 
       xlab = "", ylab="", axes = FALSE)
  # bin expression values of genes
  quantiles <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)
  for(i in 1:length(quantiles))
  {
    if(i == 1){
      tmp.id <- expressionValues < quantiles[i]
    }else if(i == length(quantiles)){
      tmp.id <- expressionValues > quantiles[i]
    }else{
      tmp.id <- expressionValues > quantiles[i] & expressionValues < quantiles[i + 1]
    }

    # plot quantiles
    means <- colMeans(codonCounts[tmp.id,], na.rm = T)
    std <- apply(codonCounts[tmp.id,], 2, sd, na.rm = T)
    for(k in 1:length(codons))
    {
      points(median(expressionValues[tmp.id]), means[k], 
             col=.codonColors[[ codons[k] ]] , pch=19, cex = 0.5)
      lines(rep(median(expressionValues[tmp.id]),2), c(means[k]-std[k], means[k]+std[k]), 
            col=.codonColors[[ codons[k] ]], lwd=0.8)
    }
  }
  
  # draw model fit
  #codonProbability <- do.call("rbind", codonProbability)
  for(i in 1:length(codons))
  {
    lines(phis, codonProbability[, i], col=.codonColors[[ codons[i] ]])
  }
  colors <- unlist(.codonColors[codons])
  
  # add indicator to optimal codon
  optim.codon.index <- which(min(c(selection, 0)) == c(selection, 0))
  codons[optim.codon.index] <- paste0(codons[optim.codon.index], "*") 
  legend("topleft", legend = codons, col=colors, bty = "n", lty=1, cex=0.75)
  
  return(xlimit)
}
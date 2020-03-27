library(AnaCoDa)
library(coda)

calcEffectiveSampleSize <- function(parameter,paramType,burn.in=1000,thin=100,mixture=1)
{
  trace <- parameter$getTraceObject()
  names.aa<- aminoAcids()
  for (aa in names.aa)
  {
    if (aa == "M" || aa == "W" || aa == "X")
      next
    codons <- AAToCodon(aa,TRUE)
    for (i in 1:length(codons))
    {
      for (j in 1:mixture)
      {
        csp.trace <-trace$getCodonSpecificParameterTraceByMixtureElementForCodon(j,codons[i],paramType,TRUE)
        csp.trace <- csp.trace[burn.in:length(csp.trace)]
        y <- mcmc(csp.trace,thin)
        cat(aa,codons[i],effectiveSize(y),"Mixture",j,"\n",sep = "\t")
      }
    }
  }
}

parameter <- loadParameterObject("test.Rda")
calcEffectiveSampleSize(parameter,1,burn.in=1,thin=50,3)
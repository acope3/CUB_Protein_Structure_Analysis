library(AnaCoDa)

calculateCovBetweenTrace <- function(trace.1,trace.2)
{
  sd.1 <- sd(trace.1)
  sd.2 <- sd(trace.2)
  corr <- cor(trace.1,trace.2)
  cov <- corr*sd.1*sd.2
  return(cov)
}


calculateCovarianceBetweenDEta <- function(parameter)
{
  trace <- parameter$getTraceObject()
  df <- data.frame(AA=c(),Total.Var=c(),Total.Covariance=c())
  aa <- aminoAcids()
  for (a in aa)
  {
    index <- 1
    if (a == "M" || a=="X" || a=="W") next
    codons <- AAToCodon(a,focal = T)
    if (length(codons) > 1)
    {
      sel.trace <- list(length=length(codons))
      total.var <- 0
      for (i in 1:length(codons))
      {
        sel.trace[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codons[i],1,T)
        total.var <- total.var + var(sel.trace[[i]])
      }
      total.cov <- 0
      for (i in 1:(length(codons)-1))
      {
        for(j in (i+1):length(codons))
        {
          cov.trace <- calculateCovBetweenTrace(sel.trace[[i]],sel.trace[[j]])
          total.cov <- total.cov + cov.trace
        }
      }
      tmp <- data.frame(AA=a,Total.Var=total.var,Total.Covariance=total.cov)
      df <- rbind(df,tmp)
    }
  }
  return(df)
}

parameter <- loadParameterObject("Scer/Predicted/Results/Ordered_disordered/Disordered/final_run/R_objects/parameter.Rda")
df <- calculateCovarianceBetweenDEta(parameter)
write.table(df,"disordered_deta_covariance.csv",sep=",",row.names=F,col.names=T,quote=F)



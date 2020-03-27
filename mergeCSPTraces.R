library(AnaCoDa)


combineTwoDimensionalTrace <- function(trace1, trace2,start=2,end=NULL)
{
  if(start < 2)
  {
    print("Start must be at least 2 because the last element of first trace is first of second trace. Setting start = 2.")
  }
  if(end <= start)
  {
    print("End must be greater than start. Setting end to length of trace2.")
    end <- trace2
  }
  if(is.null(end))
  {
    end <- length(trace2)
  }
  for (size in 1:length(trace1))
  {
    trace1[[size]]<- c(trace1[[size]], trace2[[size]][start:end])
  }
  return(trace1)
}

getTraceList <- function(trace,param.type,mixture)
{
  aa <- aminoAcids()
  traces <- vector(list,length=40)
  index <- 1
  for (a in aa)
  {
    if (a == "M" | a=="W" | a == "X") next
    codons <- AAToCodon(a,T)
    for (i in codons)
    {
      traces[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,i,param.type,T)
      index<-index+1
    }
  }
  return(traces)
}


joinCSPTracesMultipleRuns <- function(directory,param.type,starting.point,max.join=3)
{
  parameter <- loadParameterObject(paste(directory,"/run_",starting.point,"/R_objects/parameter.Rda",sep=""))
  trace <- parameter$getTraceObject()
  complete.csp.trace <- vector("list",2)
  complete.csp.trace[[1]] <- trace$getCodonSpecificParameterTrace(0)
  complete.csp.trace[[2]] <- trace$getCodonSpecificParameterTrace(1)
  for (i in starting.point-1:(1+max.join))
  {
    tmp.parameter <- loadParameterObject(paste(directory,"/run_",i,"/R_objects/parameter.Rda",sep=""))
    tmp.trace <- tmp.parameter$getTraceObject()
    for (j in 1:2)
    {
      new.trace <- tmp.trace$getCodonSpecificParameterTrace(j-1)
      numMixtures <- length(new.trace)
      for (k in 1:numMixtures)
      {
        complete.csp.trace[[j]][[k]] <- combineTwoDimensionalTrace(complete.csp.trace[[j]][[k]],new.trace[[k]],start=floor(length(new.trace[[k]][[1]])/2),end=length(new.trace[[k]][[1]]))
      }
      rm(new.trace)
    }
    rm(tmp.trace)
    rm(tmp.parameter)
  }
  trace$setCodonSpecificParameterTrace(complete.csp.trace[[1]],0)
  trace$setCodonSpecificParameterTrace(complete.csp.trace[[2]],1)
  parameter$setTraceObject(trace)
  parameter$setLastIteration(length(complete.csp.trace[[1]][[1]][[1]])-1)
  return(parameter)
}

directory <- "Final_runs/Beta/Results/Structure_Downstream_Combo/"
new.parameter <- joinCSPTracesMultipleRuns(directory,1,starting.point=7,max.join=2)
writeParameterObject(new.parameter,"Final_runs/Beta/Results/Structure_Downstream_Combo_with_fixed_mutation/combined_parameter.Rda")

fixParameter <- function(parameter.file.to.fix)
{
  load(parameter.file.to.fix)
  
  paramBase$withPhi <- F
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection", "model",  
                "mutationPrior", "mutationTrace", "selectionTrace", 
                "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
                "observedSynthesisNoiseTrace", "withPhi"),
       file=parameter.file.to.fix)
}


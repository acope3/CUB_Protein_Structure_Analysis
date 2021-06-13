##Author: Alex Cope
## Generate expected codon usage bias comparisons seen in supplementary material
library(AnaCoDa,lib="~/R_dev")

getGenesWithAA <- function(genome,aa)
{
  genes <- genome$getGenes(F)
  codons <- AAToCodon(aa=aa)
  genes.for.hist <- c()
  for (i in 1:length(genes))
  {
    seq <- genes[[i]]$seq
    seq <- unlist(strsplit(seq, ""))
    seq <- paste(seq[c(T,F,F)], seq[c(F,T,F)], seq[c(F,F,T)], sep="")
    for (codon in codons)
    {
      if (codon %in% seq)
      {
        genes.for.hist <- c(genes.for.hist,i)
        break
      }
    }
  }
  return(genes.for.hist)
}



.codonColors <- list(GCA="blue", GCC="darkorange", GCG="purple", GCT="green4", #Ala
                     TGC="darkorange", TGT="green4", #Cys
                     GAC="darkorange", GAT="green4", #Asp
                     GAA="blue", GAG="purple", #Glu
                     TTC="darkorange", TTT="green4", #Phe
                     GGA="blue", GGC="darkorange", GGG="purple", GGT="green4", #Gly
                     CAC="darkorange", CAT="green4", #His
                     ATA="blue", ATC="darkorange", ATT="green4", #Ile
                     AAA="blue", AAG="purple", #Lys
                     CTA="blue", CTC="darkorange", CTG="purple", CTT="green4", TTA="darkturquoise", TTG="deeppink3", #Leu
                     ATG="purple",
                     AAC="darkorange", AAT="green4", #Asn
                     CCA="blue", CCC="darkorange", CCG="purple", CCT="green4", #Pro
                     CAA="blue", CAG="purple", #Gln
                     CGA="blue", CGC="darkorange", CGG="purple", CGT="green4", AGA="darkturquoise", AGG="deeppink3", #Arg
                     TCA="blue", TCC="darkorange", TCG="purple", TCT="green4", #Ser4
                     ACA="blue", ACC="darkorange", ACG="purple", ACT="green4", #Thr
                     GTA="blue", GTC="darkorange", GTG="purple", GTT="green4", #Val
                     TAC="darkorange", TAT="green4", #Tyr
                     AGC="darkorange", AGT="green4", #Ser2
                     TGG="blue",
                     TAA="blue", TAG="purple", TGA="darkturquoise") #Stop


genome.1 <- initializeGenomeObject("../Data/Fasta/mod_scerevisiae_nuclear_ORF.fasta")
#genome.2 <- initializeGenomeObject("Scer/Exp_liberal_homology/Secondary_structures_pechmann/Start_helix/start_helix.fasta")
main.mutation <- read.table(file = "../mutation_mod_scerevisiae.csv", header = T, sep=",", as.is = T)
main.selection <- read.table(file = "../selection_mod_scerevisiae.csv", header = T, sep=",", as.is = T)
main.expression <- read.table(file = "../complete_seq_phi.csv", header = T, sep=",", as.is = T)[,2]
#cleft.mutation <- read.table(file = "Scer/Exp_liberal_homology/Results/Secondary_structures_pechmann/Start_helix/final_run/Parameter_est/start_helix_Mutation.csv", header = T, sep=",", as.is = T)
#cleft.selection <- read.table(file = "Scer/Exp_liberal_homology/Results/Secondary_structures_pechmann/Start_helix/final_run/Parameter_est/start_helix_Selection.csv", header = T, sep=",", as.is = T)
#cleft.expression <- read.table(file = "Scer/Exp_liberal_homology/Secondary_structures_pechmann/Start_helix/start_helix_phi.csv", header = T, sep=",", as.is = T)[,2]

mutationList = selectionList = expressionList <- vector("list", 1)
mutationList[[1]] <- main.mutation
#mutationList[[2]] <- cleft.mutation
selectionList[[1]] <- main.selection
#selectionList[[2]] <- cleft.selection
expressionList[[1]] <- main.expression
#expressionList[[2]] <- cleft.expression

nte <- read.table("../nte_pechmann.tsv",sep="\t",header=T,stringsAsFactors=F)

determineXrange <- function(expressionList)
{
  min.val <- 100000
  max.val <- -100000
  for(list.idx in 1:length(expressionList))
  {
    expression <- log10(expressionList[[list.idx]])
    min.val <- ifelse(min.val < min(expression, na.rm = T), min.val, min(expression, na.rm = T))
    max.val <- ifelse(max.val > max(expression, na.rm = T), max.val, max(expression, na.rm = T))
  }
  return(c(min.val, max.val))
}

calculateCodonProbabilities <- function(mutation, selection, expression, aa)
{
  mutation <- mutation[which(mutation[,3] !=0),]
  selection <- selection[which(selection[,3] !=0),]
  mut <- c(mutation[,3],0)
  sel <- c(selection[,3],0)
  cp <- matrix(NA, ncol=length(mut), nrow=length(expression))
  for(i in 1:length(mut))
    cp[,i] <- exp(-mut[i] - sel[i] * 10^expression)
  
  cp <- cp / rowSums(cp)
  return(cp)
}

plotOverlay <- function(mutationList, selectionList, expressionList, title, ltyList,other.metric=NULL)
{
  passed <- (length(mutationList) == length(selectionList) && length(mutationList) == length(expressionList))
  if(!passed) return()
  
  
  main <- title
  
  opar <- par(no.readonly = T)
  # mat <- matrix(rep(1,4),nrow=1,ncol=4,byrow=T)
  # mat <- rbind(mat,matrix(c(2:39,rep(40,2)),nrow = 10, ncol = 4, byrow = F))
  # mat <- rbind(mat, matrix(rep(41,4),nrow=1,ncol=4,byrow=T))
  # mat <- cbind(rep(42, 12), mat, rep(43, 12))
  # nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 5,3,5,3,5,3,5,3,5,3, 3), respect = FALSE)
  # nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 4,4,4,4,4,4,4,4,4,4, 3), respect = FALSE)
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = T)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.7, main, cex = 2)

  names.aa <- unique((mutationList[[1]])$AA)
  #names.aa <- c("A","F","K","Q","V","C","G","L","R","Y","D","H","N","S","Z","E","I","P","T")
  xlimit <- determineXrange(expressionList)
  expression <- seq(from = xlimit[1], to = xlimit[2], by = 0.01)
  max.val <- -10000
  for(list.idx in 1:length(expressionList))
  {
    hist.values <- hist(log10(expressionList[[list.idx]]), plot=FALSE, nclass=30)
    max.val <- ifelse(max(hist.values$counts) > max.val, max(hist.values$counts), max.val)
    #d <- density(log10(expressionList[[list.idx]][genes.for.hist.list[[list.idx]]]))
    #max.val <- ifelse(max(d$y) > max.val, max(d$y), max.val)
  }
  for(aa in names.aa)
  {
    
    plot(NULL, NULL, xlim=xlimit, ylim=c(-0.05,1.05), 
         xlab = "", ylab="", axes = FALSE)
    # overlay this many plots
    for(list.idx in 1:length(mutationList))
    {
      mutation <- (mutationList[[list.idx]])
      mutation <- mutation[mutation$AA == aa,]
      selection <- (selectionList[[list.idx]])
      selection <- selection[selection$AA == aa,]
      cp <- calculateCodonProbabilities(mutation, selection, expression, aa)
      codons <- AnaCoDa::AAToCodon(aa, F)
      codonCounts <- vector("list", length(codons))
	  for(i in 1:length(codons))
	  {
	    codonCounts[[i]] <- genome.1$getCodonCountsPerGene(codons[i])
	  }
	  codonCounts <- do.call("cbind", codonCounts)
	  # codon proportions
	  codonCounts <- codonCounts / rowSums(codonCounts)
	  codonCounts[is.nan(codonCounts)] <- NA
      #codons.w.ref <- AnaCoDa::AAToCodon(aa, T)
      expressionValues <- log10(expressionList[[1]])
      quantiles <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)
	  for(i in 1:length(quantiles))
	  {
	    if(i == 1){
	      tmp.id <- expressionValues< quantiles[i]
	    }else if(i == length(quantiles)){
	      tmp.id <- expressionValues> quantiles[i]
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
      for(i in 1:length(codons))
      {
        lines(expression, cp[, i], col=.codonColors[[ codons[i] ]], lty = ltyList[list.idx], lwd = 1)
      }
      colors <- unlist(.codonColors[codons])
      opt.codon <- codons[which.min(c(selection[,3],0))]
      codons.labeled <- codons
      if(!is.null(other.metric))
      {
  
      	other.metric.codons <- other.metric[which(other.metric[,1]%in% codons),]
     
      	other.metric.codons.opt <- other.metric.codons[which.max(other.metric.codons[,"Scer"]),"codons"]
      	codons.labeled[which(codons==other.metric.codons.opt)] <- paste0(codons.labeled[which(codons==other.metric.codons.opt)],"+")
      }
      codons.labeled[which(codons==opt.codon)] <- paste0(codons.labeled[which(codons==opt.codon)],"*")
      if(list.idx == length(mutationList))
      {
        legend("topleft", legend = codons.labeled, col=colors, bty = "n", lty=1, cex=0.75, lwd=1)
        
      }
    }
    box()
    main.aa <- aa #TODO map to three letter code
    text(mean(xlimit), 1, main.aa, cex = 1.5)
    if(aa %in% c("A", "F", "K", "Q", "V")){
    #if(aa %in% c("A", "C", "D", "E", "F")){
      axis(2, las=1)
    }
    if(aa %in% c("V", "Y", "Z")){
    # if(aa %in% c("F", "L", "S")){
  	 axis(1)
    }
   if(aa %in% c("A", "C", "D", "E")){
    #if(aa %in% c("A", "G","N","T")){  
      axis(3)
    }
    if(aa %in% c("E", "I", "P", "T")){
    #if(aa %in% c("T", "V", "Y", "Z")){
      axis(4, las=1)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)  
    
    #genes.for.hist.list <- list(getGenesWithAA(genome.1,aa),getGenesWithAA(genome.2,aa))
    # plot(NULL, NULL, xlim=xlimit, ylim=c(0,max.val), 
    #      xlab = "", ylab="", axes = FALSE)
    # box()
    # for(list.idx in 1:length(expressionList))
    # {
    #   hist.values <- hist(log10(expressionList[[list.idx]][genes.for.hist.list[[list.idx]]]), plot=FALSE, nclass=30)
    #   lines(hist.values$breaks, c(hist.values$counts, 0), type = "S", lwd = 1, lty = ltyList[list.idx])
    #   #lines(density(log10(expressionList[[list.idx]][genes.for.hist.list[[list.idx]]])),type = "S", lwd = 1, lty = ltyList[list.idx])
    # }
    # if(aa %in% c("F")){
    #    axis(2, las=1)
    # }
    # if(aa %in% c("V", "Y", "Z")){
    # #if(aa %in% c("F", "L", "S")){
    #   axis(1)
    # }
    # # # if(aa %in% c("A", "C", "D", "E")){
    # # if(aa %in% c("A", "G","N","T")){
    # #   axis(3)
    # # }
    # # #if(aa %in% c("E", "I", "P", "T")){
    # # if(aa %in% c("T", "V", "Y", "Z")){
    # #   axis(4, las=1,at =seq(0,max.val,max.val/2) , labels=seq(0,max.val,max.val/2))
    # # }
    # axis(1, tck = 0.02, labels = FALSE)
    # axis(2, tck = 0.02, labels = FALSE)
    # axis(3, tck = 0.02, labels = FALSE)
    # axis(4, tck = 0.02, labels = FALSE) 
  } 
  
  
  max.val <- -10000
  for(list.idx in 1:length(expressionList))
  {
    hist.values <- hist(log10(expressionList[[list.idx]]), plot=FALSE, nclass=30)
    max.val <- ifelse(max(hist.values$counts) > max.val, max(hist.values$counts), max.val)
  }
  ## adding a histogram of phi values to plot
  plot(NULL, NULL, xlim=xlimit, ylim=c(0,max.val), 
       xlab = "", ylab="", axes = FALSE)
  box()
  for(list.idx in 1:length(expressionList))
  {
    hist.values <- hist(log10(expressionList[[list.idx]]), plot=FALSE, nclass=30)
    lines(hist.values$breaks, c(hist.values$counts, 0), type = "S", lwd = 1, lty = ltyList[list.idx])
  }
  #legend("topright",legend=c("Ordered","Disordered"),lty=ltyList,cex=0.75,lwd=1)
  axis(1, las = 1)
  axis(4, las = 1)
  axis(1, tck = 0.02, labels = FALSE)
  axis(2, tck = 0.02, labels = FALSE)
  axis(3, tck = 0.02, labels = FALSE)
  axis(4, tck = 0.02, labels = FALSE)  
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.3, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"), cex = 1.5)  
  #text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.3, 0.5, "Codon Frequency", srt = 90, cex = 1.5)
}

title <- "S. cerevisiae"
ltyList <- c(1)
pdf(file = "CUB_scer_ROC_vs_nte.pdf", width = 8, height = 12)
plotOverlay(mutationList, selectionList, expressionList, title, ltyList,nte)
dev.off()


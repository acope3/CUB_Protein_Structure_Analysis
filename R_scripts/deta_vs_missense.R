library(dplyr)
library(AnaCoDa)
rescale <- function(x)
{
   x - min(x)
}

emp <- read.table("translation_errors_scer.csv",sep=",",header=T,stringsAsFactors = F)
deta <- read.table("scerevisiae_Selection.csv",sep=",",header=T,stringsAsFactors = F)

errors <- aggregate(emp$Intensities.ratio,list(emp$Origin.AA,emp$Origin.codon),FUN = function(x) c(mn = median(x), stdev = sd(x),n=length(x)))
errors <- errors[-c(1,2,3,4),]

emp_deta <- merge(errors,deta,by.x="Group.2",by.y="Codon",all=T)

present <- emp_deta %>% group_by(AA) %>% filter(length(which(!is.na(x)))>=2)
present <- present[which(!is.na(present$x)),]
error_table<-present %>% group_by(AA) %>% mutate(Relative.Error=rescale(x),Posterior.scaled=rescale(Posterior))                          

aa <- aminoAcids()
pdf("scer_relative_error_rates_vs_deta.pdf")
for (a in aa)
{
  if (a %in% unlist(error_table$AA))
  {
    index <- which(error_table$AA == a)
    plot(unlist(error_table[index,"Relative.Error"]),unlist(error_table[index,"Posterior.scaled"]),col="lightblue",main=a,xlab="Relative Error Rate",ylab="Delta.Eta")
    text(unlist(error_table[index,"Relative.Error"]),unlist(error_table[index,"Posterior.scaled"]),labels=unlist(error_table[index,"Group.2"]), cex=0.9, font=2)
  }
}
dev.off()
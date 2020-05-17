library(AnaCoDa)
df <- read.table("../Data/scer_mh_test_table_sim_ss.tsv",header=T,stringsAsFactors = F,sep="\t")
# phi <- read.table("mod_scerevisiae_expression_wo_phi_allunique.csv",header=T,stringsAsFactors = F,sep=",")
# #phi[,1]<-stringr::str_extract(phi[,1],pattern="NP_[0-9]+.[0-9]+")
# top.10.cutoff <- quantile(phi[,2],prob=0.60)
# top.10.prot <- phi[which(phi$Phi > top.10.cutoff),"Gene_id"]
# df <- df[which(df$Protein %in% top.10.prot),]
struct <- unique(df$Class)
names.aa <- aminoAcids()
p.val.list <- c()
codon.list <- c()
struct.list <- c()
odds.ratio <- c()
for (a in names.aa)
{
  if (a == "M" || a == "W" || a == "X") next
  df.aa <- df[which(df$AA == a),]
  codon.family <- AAToCodon(a,F)
  for (codon in codon.family)
  {
    df.aa.tmp <- df.aa
    df.aa.tmp[which(df.aa.tmp$Codon != codon),"Codon"] <- "Other Codon"
    for (s in struct)
    {
      df.aa.tmp.2 <- df.aa.tmp
      df.aa.tmp.2[which(df.aa.tmp.2$Class != s),"Class"] <- "Other Class"
      ctable <- table(df.aa.tmp.2$Codon,df.aa.tmp.2$Class,df.aa.tmp.2$Protein)
      ctable <- ctable[,,which(apply(ctable,MARGIN =3,FUN=sum)>1)]
      #if (length(dim(ctable)) != 3) next
      
      mh <- mantelhaen.test(ctable,exact = T)
      p.val.list <- c(p.val.list,mh$p.value)
      codon.list <- c(codon.list,codon)
      struct.list <- c(struct.list,s)
      odds.ratio <- c(odds.ratio,mh$estimate)
    }
    
  }
}
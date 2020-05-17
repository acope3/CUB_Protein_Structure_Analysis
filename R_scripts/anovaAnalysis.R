library(plyr)
library(sjPlot)
library(car)
normalize <- function(dEta)
{
  aa <- AnaCoDa::aminoAcids()
 

  for (a in aa)
  {
    if (a=="W" || a=="M"||a=="X") next
    row.ind <- which(dEta[,1] == a)
    rows <- dEta[row.ind,]
    mu <- mean(rows[,"Posterior"])
    rows[,c("Posterior","X0.025.","X0.975.")] <- rows[,c("Posterior","X0.025.","X0.975.")] - mu
    dEta[row.ind,] <- rows
  }
  return(dEta)
}



coil.ord <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Coil_ord/final_run/Parameter_est/coil_ord_Selection.csv",sep=",",header=T,stringsAsFactors=F)
helix.ord <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Helix_ord/final_run/Parameter_est/helix_ord_Selection.csv",sep=",",header=T,stringsAsFactors=F)
sheet.ord <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Sheet_ord/final_run/Parameter_est/sheet_ord_Selection.csv",sep=",",header=T,stringsAsFactors=F)

coil.dis <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Coil_dis/final_run/Parameter_est/coil_dis_Selection.csv",sep=",",header=T,stringsAsFactors=F)
helix.dis <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Helix_dis/final_run/Parameter_est/helix_dis_Selection.csv",sep=",",header=T,stringsAsFactors=F)
sheet.dis <- read.table("../Scer/Predicted/Results/Secondary_structure_order/Sheet_dis/final_run/Parameter_est/sheet_dis_Selection.csv",sep=",",header=T,stringsAsFactors=F)

ord <- read.table("../Scer/Predicted/Results/Ordered_disordered/Ordered/final_run/Parameter_est/ordered_Selection.csv",sep=",",header=T,stringsAsFactors=F)
dis <- read.table("../Scer/Predicted/Results/Ordered_disordered/Disordered/final_run/Parameter_est/disordered_Selection.csv",sep=",",header=T,stringsAsFactors=F)

coil <- read.table("../Scer/Predicted/Results/Secondary_structures/Coil/final_run/Parameter_est/coil_Selection.csv",sep=",",header=T,stringsAsFactors=F)
helix <- read.table("../Scer/Predicted/Results/Secondary_structures/Helix/final_run/Parameter_est/helix_Selection.csv",sep=",",header=T,stringsAsFactors=F)
sheet <- read.table("../Scer/Predicted/Results/Secondary_structures/Sheet/final_run/Parameter_est/sheet_Selection.csv",sep=",",header=T,stringsAsFactors=F)


# coil.ord <- normalize(coil.ord)
# helix.ord <- normalize(helix.ord)
# sheet.ord <- normalize(sheet.ord)

# coil.dis <- normalize(coil.dis)
# helix.dis <- normalize(helix.dis)
# sheet.dis <- normalize(sheet.dis)



coil.ord["Order"] <- "Ordered"
helix.ord["Order"] <- "Ordered"
sheet.ord["Order"] <- "Ordered"

coil.dis["Order"] <- "Disordered"
helix.dis["Order"] <- "Disordered"
sheet.dis["Order"] <- "Disordered"



coil.ord["SS"] <- "Coil"
coil.dis["SS"] <- "Coil"

helix.ord["SS"] <- "Helix"
helix.dis["SS"] <- "Helix"

sheet.ord["SS"] <- "Sheet"
sheet.dis["SS"] <- "Sheet"

ord["Order"] <- "Ordered"
ord["SS"] <- "Helix_Sheet_Coil"
dis["Order"] <- "Disordered"
dis["SS"] <- "Helix_Sheet_Coil"

coil["SS"] <- "Coil"
coil["Order"] <- "Ordered_Disordered"
helix["SS"] <- "Helix"
helix["Order"] <- "Ordered_Disordered"
sheet["SS"] <- "Sheet"
sheet["Order"] <- "Ordered_Disordered"



df <- do.call("rbind",list(coil.ord,helix.ord,sheet.ord,coil.dis,helix.dis,sheet.dis))

df["Properties"] <- revalue(df$AA,c(
                              "A"="Hydrophobic",
                              "C"="Other",
                              "D"="Charged",
                              "E"="Charged",
                              "F"="Hydrophobic",
                              "G"="Other",
                              "H"="Charged",
                              "I"="Hydrophobic",
                              "K"="Charged",
                              "L"="Hydrophobic",
                              "N"="Polar",
                              "P"="Other",
                              "Q"="Polar",
                              "R"="Charged",
                              "S"="Polar",
                              "T"="Polar",
                              "V"="Hydrophobic",
                              "Y"="Hydrophobic",
                              "Z"="Polar"))

df <- df[which(df$Posterior != 0),]

m1 <- aov(Posterior ~ Properties * Order * SS,data=df)
m2 <- aov(Posterior ~ Properties + Order + SS,data=df)
comp <- anova(m2,m1)

type.2 <- Anova(lm(Posterior ~ Order + SS + Properties,data=df),type=2)
print(summary(type.2))

pdf("anova_interaction.pdf")
interaction.plot(df$Order,df$SS,df$Posterior)
interaction.plot(df$SS,df$Order,df$Posterior)
interaction.plot(df$Order,df$Properties,df$Posterior)
interaction.plot(df$SS,df$Properties,df$Posterior)
dev.off()
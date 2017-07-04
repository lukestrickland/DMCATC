# # # Manifests and LMER Analyses # # #
#
#     Analysis of manifest RT and Response Proportion by PM block, 
#     Time Pressure condition and Stimulus Type
#
#     Analysis of model parameters between experiments, i.e., comparison of 
#     task importance manipulations
#

rm(list=ls())
# setwd("~/russ_model_analyses")
setwd("D:/Software/DMC_ATCPMDC")
source("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")
source("LSAnova.R")
require("lme4")
require(plyr)
require(dplyr)
# require("pander")
load("~/Modelling/x1/samples/okdats.A1.RData")  # Original data
load("~/A1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.A1 <- A1.block.B.V_cond.B.V.PMV.samples
rm(A1.block.B.V_cond.B.V.PMV.samples)


# # # Get data from samples object # # #
#
get.hdata.dmc <- function(hsamples){
  list.wind<-lapply(seq_along(hsamples), function(samples, n, i) cbind(n[[i]], samples[[i]]$data), 
                    samples= hsamples, n = names(hsamples))
  out<-do.call(rbind, list.wind)
  names(out)[1] <- "s"
  out
}
#

data.A1 <- get.hdata.dmc(samples.A1)  # Get data from samples object 
save(data.A1, file="data.A1.RData")

# all.equal(data1[order(data1$s),], data2)  # Check recovered data matches original data

head(data.A1)
tail(data.A1)

# # # Add logical S-R match factor 'C' # # #
#
data.A1$C <- rep(0,length(data.A1$RT))
for(i in 1:length(data.A1$RT)){
  if(data.A1$S[i]=="cc" & data.A1$R[i]=="C"){
    data.A1$C[i] <- 1
  } else if(data.A1$S[i]=="nn" & data.A1$R[i]=="N"){
    data.A1$C[i] <- 1
  } else if(data.A1$S[i]=="pc" & data.A1$R[i]=="P"){
    data.A1$C[i] <- 1
  } else if(data.A1$S[i]=="pn" & data.A1$R[i]=="P"){
    data.A1$C[i] <- 1
  } 
}
#


# # # Prep dataframes for analysis # # #
#
cdt <- data.A1[!(data.A1$S=="pc" | data.A1$S=="pn"),]  # Conflict detection task only (no PM)
cdt$S <- factor(as.character(cdt$S)); cdt$R <- factor(as.character(cdt$R)) 
str(cdt)

pm <- data.A1[(data.A1$S=="pc" | data.A1$S=="pn"),]  # PM trials only
pm$S <- factor(as.character(pm$S)); pm$block <- factor(as.character(pm$block))
str(pm)
#

pm.glmer.A1 <- glmer(C ~ S*cond+(1|s), data=pm, family=binomial(link="probit"))
save(pm.glmer.A1, file="pm.glmer.A1.RData")
load("pm.glmer.A1.RData")
pm.glm.A1 <- Anova(pm.glmer.A1,type="II")
pm.glm.A1
PM.TABLE.A1 <- data.frame(pm.glm.A1$Chisq, pm.glm.A1$Df, pm.glm.A1$Pr)
PM.TABLE.A1$pm.glm.A1.Chisq <- round(PM.TABLE.A1$pm.glm.A1.Chisq, 2)
PM.TABLE.A1$pm.glm.A1.Pr <- format.pval(PM.TABLE.A1$pm.glm.A1.Pr, digits=2, eps= 0.001)
PM.TABLE.A1$pm.glm.A1.Pr <- gsub("0\\.", ".", PM.TABLE.A1$pm.glm.A1.Pr)
rownames(PM.TABLE.A1) <- c("Stimulus","TP","Stimulus*TP")
PM.TABLE.A1

cdt.glmer.A1 <- glmer(C ~ block*cond+(1|s), data=cdt, family=binomial(link="probit"))
save(cdt.glmer.A1, file="cdt.glmer.A1.RData")
load("cdt.glmer.A1.RData")
cdt.glm.A1 <- Anova(cdt.glmer.A1,type="II")
cdt.glm.A1
CDT.TABLE.A1 <- data.frame(cdt.glm.A1$Chisq, cdt.glm.A1$Df, cdt.glm.A1$Pr)
CDT.TABLE.A1$cdt.glm.A1.Chisq <- round(CDT.TABLE.A1$cdt.glm.A1.Chisq, 2)
CDT.TABLE.A1$cdt.glm.A1.Pr <- format.pval(CDT.TABLE.A1$cdt.glm.A1.Pr, digits=2, eps= 0.001)
CDT.TABLE.A1$cdt.glm.A1.Pr <- gsub("0\\.", ".", CDT.TABLE.A1$cdt.glm.A1.Pr)
rownames(CDT.TABLE.A1) <- c("PM Block","TP","PM Block*TP")
CDT.TABLE.A1

# # # Prep RT dataframes for analysis (keep correct RTs only) # # #
#
dPM.RT <- arr2df(tapply(pm$RT[pm$C==1], 
                        list(s=pm$s[pm$C==1], block=pm$block[pm$C==1], cond=pm$cond[pm$C==1], S=pm$S[pm$C==1]),
                        mean))
dPM.RT
str(dPM.RT)
#

PMC.RT.lmer.A1 <- lmer(y ~ S*cond+(1|s), data=dPM.RT)
save(PMC.RT.lmer.A1, file="PMC.RT.lmer.A1.RData")
load("PMC.RT.lmer.A1.RData")
PM.RT.glm.A1 <- Anova(PMC.RT.lmer.A1,type="II")
PM.RT.glm.A1
PM.RT.TABLE.A1 <- data.frame(PM.RT.glm.A1$Chisq, PM.RT.glm.A1$Df, PM.RT.glm.A1$Pr)
PM.RT.TABLE.A1$PM.RT.glm.A1.Chisq <- round(PM.RT.TABLE.A1$PM.RT.glm.A1.Chisq, 2)
PM.RT.TABLE.A1$PM.RT.glm.A1.Pr <- format.pval(PM.RT.TABLE.A1$PM.RT.glm.A1.Pr, digits=2, eps= 0.001)
PM.RT.TABLE.A1$PM.RT.glm.A1.Pr <- gsub("0\\.", ".", PM.RT.TABLE.A1$PM.RT.glm.A1.Pr)
rownames(PM.RT.TABLE.A1) <- c("Stimulus","TP","Stimulus*TP")
PM.RT.TABLE.A1

# 
dCDT.RT <- arr2df(tapply(cdt$RT[cdt$C==1],
                         list(s=cdt$s[cdt$C==1], block=cdt$block[cdt$C==1], cond=cdt$cond[cdt$C==1], S=cdt$S[cdt$C==1]),
                         mean))
str(dCDT.RT)
#

CDTC.RT.lmer.A1 <- lmer(y ~ S*block*cond+(1|s), data=dCDT.RT)
save(CDTC.RT.lmer.A1, file="CDTC.RT.lmer.A1.RData")
load("CDTC.RT.lmer.A1.RData")
CDT.RT.glm.A1 <- Anova(CDTC.RT.lmer.A1,type="II")
CDT.RT.glm.A1
CDT.RT.TABLE.A1 <- data.frame(CDT.RT.glm.A1$Chisq, CDT.RT.glm.A1$Df, CDT.RT.glm.A1$Pr)
CDT.RT.TABLE.A1$CDT.RT.glm.A1.Chisq <- round(CDT.RT.TABLE.A1$CDT.RT.glm.A1.Chisq, 2)
CDT.RT.TABLE.A1$CDT.RT.glm.A1.Pr <- format.pval(CDT.RT.TABLE.A1$CDT.RT.glm.A1.Pr, digits=2, eps= 0.001)
CDT.RT.TABLE.A1$CDT.RT.glm.A1.Pr <- gsub("0\\.", ".", CDT.RT.TABLE.A1$CDT.RT.glm.A1.Pr)
rownames(CDT.RT.TABLE.A1) <- c("Stimulus","PM Block","TP","Stimulus*Block","Stimulus*TP","Block*TP","Stimulus*Block*TP")
CDT.RT.TABLE.A1 


# # # Manifests - Mean Reponse Proportion # # #
#

RP.A1 <- ddply(data.A1, ~s*cond*block*S, summarise, RP=mean(C,na.rm=TRUE))  # Get Response Proportions
head(RP.A1)
str(RP.A1)
save(RP.A1, file="RP.A1.RData")

mRP.A1 <- ddply(RP.A1, ~cond*block*S, summarise, M=mean(RP,na.rm=TRUE), SD=sd(RP,na.rm=TRUE))  # Means & SDs for Response Proportions
head(mRP.A1)
str(mRP.A1)
save(mRP.A1, file="mRP.A1.RData")

# # # Manifests - Mean RT # # #
#

mRT.A1 <- ddply(data.A1, ~s*cond*block*S, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))  # Mean RT
head(mRT.A1)
save(mRT.A1, file="mRT.A1.RData")





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
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")
source("LSAnova.R")
require("lme4")
require(plyr)
require(dplyr)
require("pander")
# load("~/Modelling/x1/samples/okdats.E1.RData")  # Original data
load("C:/Users/Russell Boag/Documents/GitHub/DMCATC/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.E1 <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)


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

data.E1 <- get.hdata.dmc(samples.E1)  # Get data from samples object
save(data.E1, file="data.E1.RData")

# all.equal(datE1[order(datE1$s),], data2)  # Check recovered data matches original data

head(data.E1)
tail(data.E1)

# # # Add logical S-R match factor 'C' # # #
#
data.E1$C <- rep(0,length(data.E1$RT))
for(i in 1:length(data.E1$RT)){
  if(data.E1$S[i]=="cc" & data.E1$R[i]=="C"){
    data.E1$C[i] <- 1
  } else if(data.E1$S[i]=="nn" & data.E1$R[i]=="N"){
    data.E1$C[i] <- 1
  } else if(data.E1$S[i]=="pc" & data.E1$R[i]=="P"){
    data.E1$C[i] <- 1
  } else if(data.E1$S[i]=="pn" & data.E1$R[i]=="P"){
    data.E1$C[i] <- 1
  }
}
#


# # # Prep dataframes for analysis # # #
#
cdt <- data.E1[!(data.E1$S=="pc" | data.E1$S=="pn"),]  # Conflict detection task only (no PM)
cdt$S <- factor(as.character(cdt$S)); cdt$R <- factor(as.character(cdt$R))
str(cdt)

pm <- data.E1[(data.E1$S=="pc" | data.E1$S=="pn"),]  # PM trials only
pm$S <- factor(as.character(pm$S)); pm$block <- factor(as.character(pm$block))
str(pm)
#

pm.glmer.E1 <- glmer(C ~ S*cond+(1|s), data=pm, family=binomial(link="probit"))
save(pm.glmer.E1, file="pm.glmer.E1.RData")
load("pm.glmer.E1.RData")
pm.glm.E1 <- Anova(pm.glmer.E1,type="II")
pm.glm.E1
PM.TABLE.E1 <- data.frame(pm.glm.E1$Chisq, pm.glm.E1$Df, pm.glm.E1$Pr)
PM.TABLE.E1$pm.glm.E1.Chisq <- round(PM.TABLE.E1$pm.glm.E1.Chisq, 2)
PM.TABLE.E1$pm.glm.E1.Pr <- format.pval(PM.TABLE.E1$pm.glm.E1.Pr, digits=2, eps= 0.001)
PM.TABLE.E1$pm.glm.E1.Pr <- gsub("0\\.", ".", PM.TABLE.E1$pm.glm.E1.Pr)
rownames(PM.TABLE.E1) <- c("Stimulus","TP","Stimulus*TP")
PM.TABLE.E1

cdt.glmer.E1 <- glmer(C ~ block*cond+(1|s), data=cdt, family=binomial(link="probit"))
save(cdt.glmer.E1, file="cdt.glmer.E1.RData")
load("cdt.glmer.E1.RData")
cdt.glm.E1 <- Anova(cdt.glmer.E1,type="II")
cdt.glm.E1
CDT.TABLE.E1 <- data.frame(cdt.glm.E1$Chisq, cdt.glm.E1$Df, cdt.glm.E1$Pr)
CDT.TABLE.E1$cdt.glm.E1.Chisq <- round(CDT.TABLE.E1$cdt.glm.E1.Chisq, 2)
CDT.TABLE.E1$cdt.glm.E1.Pr <- format.pval(CDT.TABLE.E1$cdt.glm.E1.Pr, digits=2, eps= 0.001)
CDT.TABLE.E1$cdt.glm.E1.Pr <- gsub("0\\.", ".", CDT.TABLE.E1$cdt.glm.E1.Pr)
rownames(CDT.TABLE.E1) <- c("PM Block","TP","PM Block*TP")
CDT.TABLE.E1

# # # Prep RT dataframes for analysis (keep correct RTs only) # # #
#
dPM.RT <- arr2df(tapply(pm$RT[pm$C==1],
                        list(s=pm$s[pm$C==1], block=pm$block[pm$C==1], cond=pm$cond[pm$C==1], S=pm$S[pm$C==1]),
                        mean))
dPM.RT
str(dPM.RT)
#

PMC.RT.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=dPM.RT)
save(PMC.RT.lmer.E1, file="PMC.RT.lmer.E1.RData")
load("PMC.RT.lmer.E1.RData")
PM.RT.glm.E1 <- Anova(PMC.RT.lmer.E1,type="II")
PM.RT.glm.E1
PM.RT.TABLE.E1 <- data.frame(PM.RT.glm.E1$Chisq, PM.RT.glm.E1$Df, PM.RT.glm.E1$Pr)
PM.RT.TABLE.E1$PM.RT.glm.E1.Chisq <- round(PM.RT.TABLE.E1$PM.RT.glm.E1.Chisq, 2)
PM.RT.TABLE.E1$PM.RT.glm.E1.Pr <- format.pval(PM.RT.TABLE.E1$PM.RT.glm.E1.Pr, digits=2, eps= 0.001)
PM.RT.TABLE.E1$PM.RT.glm.E1.Pr <- gsub("0\\.", ".", PM.RT.TABLE.E1$PM.RT.glm.E1.Pr)
rownames(PM.RT.TABLE.E1) <- c("Stimulus","TP","Stimulus*TP")
PM.RT.TABLE.E1

#
dCDT.RT <- arr2df(tapply(cdt$RT[cdt$C==1],
                         list(s=cdt$s[cdt$C==1], block=cdt$block[cdt$C==1], cond=cdt$cond[cdt$C==1], S=cdt$S[cdt$C==1]),
                         mean))
str(dCDT.RT)
#

CDTC.RT.lmer.E1 <- lmer(y ~ S*block*cond+(1|s), data=dCDT.RT)
save(CDTC.RT.lmer.E1, file="CDTC.RT.lmer.E1.RData")
load("CDTC.RT.lmer.E1.RData")
CDT.RT.glm.E1 <- Anova(CDTC.RT.lmer.E1,type="II")
CDT.RT.glm.E1
CDT.RT.TABLE.E1 <- data.frame(CDT.RT.glm.E1$Chisq, CDT.RT.glm.E1$Df, CDT.RT.glm.E1$Pr)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq <- round(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq, 2)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- format.pval(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr, digits=2, eps= 0.001)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- gsub("0\\.", ".", CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr)
rownames(CDT.RT.TABLE.E1) <- c("Stimulus","PM Block","TP","Stimulus*Block","Stimulus*TP","Block*TP","Stimulus*Block*TP")
CDT.RT.TABLE.E1


# # # Manifests - Mean Reponse Proportion # # #
#

RP.E1 <- ddply(data.E1, ~s*cond*block*S, summarise, RP=mean(C,na.rm=TRUE))  # Get Response Proportions
head(RP.E1)
str(RP.E1)
save(RP.E1, file="RP.E1.RData")

mRP.E1 <- ddply(RP.E1, ~cond*block*S, summarise, M=mean(RP,na.rm=TRUE), SD=sd(RP,na.rm=TRUE))  # Means & SDs for Response Proportions
head(mRP.E1)
str(mRP.E1)
save(mRP.E1, file="mRP.E1.RData")

# # # Manifests - Mean RT # # #
#

mRT.E1 <- ddply(data.E1, ~s*cond*block*S, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))  # Mean RT
head(mRT.E1)
save(mRT.E1, file="mRT.E1.RData")





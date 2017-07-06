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
require(gridExtra)
require("lme4")
require(car)
require(plyr)
require(dplyr)
require("pander")
# load("~/Modelling/x1/samples/okdats.E1.RData")  # Original data
load("C:/Users/Russell Boag/Documents/GitHub/DMCATC/data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
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
head(data.E1)
# save(data.E1, file="data.E1.RData")

# all.equal(datE1[order(datE1$s),], data2)  # Check recovered data matches original data


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
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC/analysis")

# Conflict Detection Task trials only (no PM)
CDT <- data.E1[!(data.E1$S=="pc" | data.E1$S=="pn"),]
CDT$S <- factor(as.character(CDT$S)); CDT$R <- factor(as.character(CDT$R))
str(CDT)
head(CDT)

# PM Task trials only
PMT <- data.E1[(data.E1$S=="pc" | data.E1$S=="pn"),]
PMT$S <- factor(as.character(PMT$S)); PMT$block <- factor(as.character(PMT$block))
str(PMT)
head(PMT)
#

CDT.glmer.E1 <- glmer(C ~ S*block*cond+(1|s), data=CDT, family=binomial(link="probit"))
save(CDT.glmer.E1, file="CDT.glmer.E1.RData")
load("CDT.glmer.E1.RData")
CDT.glm.E1 <- Anova(CDT.glmer.E1,type="II")
CDT.glm.E1
# CDT.TABLE.E1 <- data.frame(cdt.glm.E1$Chisq, cdt.glm.E1$Df, cdt.glm.E1$Pr)
# CDT.TABLE.E1$cdt.glm.E1.Chisq <- round(CDT.TABLE.E1$cdt.glm.E1.Chisq, 2)
# CDT.TABLE.E1$cdt.glm.E1.Pr <- format.pval(CDT.TABLE.E1$cdt.glm.E1.Pr, digits=2, eps= 0.001)
# CDT.TABLE.E1$cdt.glm.E1.Pr <- gsub("0\\.", ".", CDT.TABLE.E1$cdt.glm.E1.Pr)
# # rownames(CDT.TABLE.E1) <- c("PM Block","TP","PM Block*TP")
# CDT.TABLE.E1

PMT.glmer.E1 <- glmer(C ~ S*cond+(1|s), data=PMT, family=binomial(link="probit"))
save(PMT.glmer.E1, file="PMT.glmer.E1.RData")
load("PMT.glmer.E1.RData")
PMT.glm.E1 <- Anova(PMT.glmer.E1,type="II")
PMT.glm.E1
# PM.TABLE.E1 <- data.frame(pm.glm.E1$Chisq, pm.glm.E1$Df, pm.glm.E1$Pr)
# PM.TABLE.E1$pm.glm.E1.Chisq <- round(PM.TABLE.E1$pm.glm.E1.Chisq, 2)
# PM.TABLE.E1$pm.glm.E1.Pr <- format.pval(PM.TABLE.E1$pm.glm.E1.Pr, digits=2, eps= 0.001)
# PM.TABLE.E1$pm.glm.E1.Pr <- gsub("0\\.", ".", PM.TABLE.E1$pm.glm.E1.Pr)
# rownames(PM.TABLE.E1) <- c("Stimulus","TP","Stimulus*TP")
# PM.TABLE.E1


# # # Prep RT dataframes for analysis (keep correct RTs only) # # #
#
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC/analysis")

# Conflict Detection Task correct response trials only (no PM)
dCDT.RT <- arr2df(tapply(CDT$RT[CDT$C==1],
                         list(s=CDT$s[CDT$C==1], block=CDT$block[CDT$C==1], cond=CDT$cond[CDT$C==1], S=CDT$S[CDT$C==1]),
                         mean))
head(dCDT.RT)
str(dCDT.RT)

# PM Task correct response trials only
dPMT.RT <- arr2df(tapply(PMT$RT[PMT$C==1],
                        list(s=PMT$s[PMT$C==1], block=PMT$block[PMT$C==1], cond=PMT$cond[PMT$C==1], S=PMT$S[PMT$C==1]),
                        mean))
head(dPMT.RT)
str(dPMT.RT)



#
CDTC.RT.lmer.E1 <- lmer(y ~ S*block*cond+(1|s), data=dCDT.RT)
save(CDTC.RT.lmer.E1, file="CDTC.RT.lmer.E1.RData")
load("CDTC.RT.lmer.E1.RData")
CDT.RT.glm.E1 <- Anova(CDTC.RT.lmer.E1,type="II")
CDT.RT.glm.E1
# CDT.RT.TABLE.E1 <- data.frame(CDT.RT.glm.E1$Chisq, CDT.RT.glm.E1$Df, CDT.RT.glm.E1$Pr)
# CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq <- round(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq, 2)
# CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- format.pval(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr, digits=2, eps= 0.001)
# CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- gsub("0\\.", ".", CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr)
# rownames(CDT.RT.TABLE.E1) <- c("Stimulus","PM Block","TP","Stimulus*Block","Stimulus*TP","Block*TP","Stimulus*Block*TP")
# CDT.RT.TABLE.E1

#
PMC.RT.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=dPM.RT)
save(PMC.RT.lmer.E1, file="PMC.RT.lmer.E1.RData")
load("PMC.RT.lmer.E1.RData")
PM.RT.glm.E1 <- Anova(PMC.RT.lmer.E1,type="II")
PM.RT.glm.E1
# PM.RT.TABLE.E1 <- data.frame(PM.RT.glm.E1$Chisq, PM.RT.glm.E1$Df, PM.RT.glm.E1$Pr)
# PM.RT.TABLE.E1$PM.RT.glm.E1.Chisq <- round(PM.RT.TABLE.E1$PM.RT.glm.E1.Chisq, 2)
# PM.RT.TABLE.E1$PM.RT.glm.E1.Pr <- format.pval(PM.RT.TABLE.E1$PM.RT.glm.E1.Pr, digits=2, eps= 0.001)
# PM.RT.TABLE.E1$PM.RT.glm.E1.Pr <- gsub("0\\.", ".", PM.RT.TABLE.E1$PM.RT.glm.E1.Pr)
# rownames(PM.RT.TABLE.E1) <- c("Stimulus","TP","Stimulus*TP")
# PM.RT.TABLE.E1


# # # Manifests - Mean Reponse Proportion # # #
#
length(CDT$RT[CDT$S=="cc" & CDT$R=="C"])/
    length(CDT$RT[CDT$S=="cc"])  # CORRECT ACCURACY

RP.CDT <- ddply(CDT, ~s*cond*block*S, summarise, RP=mean(C))
RP.CDT <- ddply(CDT[CDT$S=="cc",], ~s*cond*block, summarise, RP=mean(C))

mean(CDT$C[CDT$S=="cc"])  # CORRECT ACCURACY
mean(RP.CDT$RP[RP.CDT$S=="cc"])

RP.CDT.S <- ddply(RP.CDT, ~S, summarise, M=mean(RP), SD=sd(RP))
RP.CDT.Block <- ddply(RP.CDT, ~S*block, summarise, M=mean(RP), SD=sd(RP))
RP.CDT.TP <- ddply(RP.CDT, ~S*cond, summarise, M=mean(RP), SD=sd(RP))
RP.CDT.Block.TP <- ddply(RP.CDT, ~S*block*cond, summarise, M=mean(RP), SD=sd(RP))

RP.PMT <- ddply(PMT, ~s*cond*block*S, summarise, RP=mean(C))
RP.PMT.S <- ddply(RP.PMT, ~S, summarise, M=mean(RP), SD=sd(RP))
RP.PMT.TP <- ddply(RP.PMT, ~cond, summarise, M=mean(RP), SD=sd(RP))

RP.CDT.S
RP.CDT.Block
RP.CDT.TP
RP.CDT.Block.TP

RP.PMT.S
RP.PMT.TP





# Means & SDs for Response Proportions
mRP.E1.Block <- ddply(RP.E1, ~block, summarise, M=mean(RP,na.rm=TRUE), SD=sd(RP,na.rm=TRUE))
mRP.E1.Block
mRP.E1.TP <- ddply(RP.E1, ~cond, summarise, M=mean(RP,na.rm=TRUE), SD=sd(RP,na.rm=TRUE))
mRP.E1.TP

head(mRP.E1)
str(mRP.E1)



# save(mRP.E1, file="mRP.E1.RData")

# # # Manifests - Mean RT # # #
#
# All responses
mRT.E1 <- ddply(data.E1[ (data.E1$S=="cc" & data.E1$C=="1") | (data.E1$S=="nn" & data.E1$C=="1"), ], ~cond*S, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))  # Mean RT
mRT.E1


mRT.E1 <- ddply(data.E1[data.E1$C=="1", ], ~S*cond, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))  # Mean RT for Correct Responses
mRT.E1
mRT.E1 <- ddply(data.E1[data.E1$C=="1", ], ~S*block, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))  # Mean RT for Correct Responses
mRT.E1

# save(mRT.E1, file="mRT.E1.RData")
dCDT.RT
# CDT Correct Responses Only
ddply(dCDT.RT, ~S*block, summarise, M=mean(y, na.rm=TRUE),SD=sd(y, na.rm = TRUE))

# PM Task Correct Responses Only
ddply(dPMT.RT, ~S*cond, summarise, M=mean(y, na.rm=TRUE))



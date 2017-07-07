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
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "car", "lme4", "xtable")
# install.packages(pkgs) #install
sapply(pkgs, require, character.only = T) #load
# load("~/Modelling/x1/samples/okdats.A4.RData")  # Original data
load("C:/Users/Russell Boag/Documents/GitHub/DMCATC/data/samples/A4.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.A4 <- A4.block.B.V_cond.B.V.PMV.samples
rm(A4.block.B.V_cond.B.V.PMV.samples)


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

data.A4 <- get.hdata.dmc(samples.A4)  # Get data from samples object
head(data.A4)
# save(data.A4, file="data.A4.RData")

# all.equal(datA4[order(datA4$s),], data2)  # Check recovered data matches original data


# # # Add logical S-R match factor 'C' # # #
#
data.A4$C <- rep(0,length(data.A4$RT))
for(i in 1:length(data.A4$RT)){
  if(data.A4$S[i]=="cc" & data.A4$R[i]=="C"){
    data.A4$C[i] <- 1
  } else if(data.A4$S[i]=="nn" & data.A4$R[i]=="N"){
    data.A4$C[i] <- 1
  } else if(data.A4$S[i]=="pc" & data.A4$R[i]=="P"){
    data.A4$C[i] <- 1
  } else if(data.A4$S[i]=="pn" & data.A4$R[i]=="P"){
    data.A4$C[i] <- 1
  }
}
#


# # # Prep dataframes for analysis # # #
#
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC/analysis")

# Conflict Detection Task trials only (no PM)
CDT <- data.A4[!(data.A4$S=="pc" | data.A4$S=="pn"),]
CDT$S <- factor(as.character(CDT$S)); CDT$R <- factor(as.character(CDT$R))
str(CDT)
head(CDT)

# PM Task trials only
PMT <- data.A4[(data.A4$S=="pc" | data.A4$S=="pn"),]
PMT$S <- factor(as.character(PMT$S)); PMT$block <- factor(as.character(PMT$block))
str(PMT)
head(PMT)
#

CDT.Acc.glmer.A4 <- glmer(C ~ S*block*cond+(1|s), data=CDT, family=binomial(link="probit"))
save(CDT.Acc.glmer.A4, file="CDT.Acc.glmer.A4.RData")
load("CDT.Acc.glmer.A4.RData")
CDT.Acc.glm.A4 <- Anova(CDT.Acc.glmer.A4,type="II")
CDT.Acc.glm.A4
# CDT.TABLE.A4 <- data.frame(cdt.glm.A4$Chisq, cdt.glm.A4$Df, cdt.glm.A4$Pr)
# CDT.TABLE.A4$cdt.glm.A4.Chisq <- round(CDT.TABLE.A4$cdt.glm.A4.Chisq, 2)
# CDT.TABLE.A4$cdt.glm.A4.Pr <- format.pval(CDT.TABLE.A4$cdt.glm.A4.Pr, digits=2, eps= 0.001)
# CDT.TABLE.A4$cdt.glm.A4.Pr <- gsub("0\\.", ".", CDT.TABLE.A4$cdt.glm.A4.Pr)
# # rownames(CDT.TABLE.A4) <- c("PM Block","TP","PM Block*TP")
# CDT.TABLE.A4

PMT.Acc.glmer.A4 <- glmer(C ~ S*cond+(1|s), data=PMT, family=binomial(link="probit"))
save(PMT.Acc.glmer.A4, file="PMT.Acc.glmer.A4.RData")
load("PMT.Acc.glmer.A4.RData")
PMT.Acc.glm.A4 <- Anova(PMT.Acc.glmer.A4,type="II")
PMT.Acc.glm.A4
# PM.TABLE.A4 <- data.frame(pm.glm.A4$Chisq, pm.glm.A4$Df, pm.glm.A4$Pr)
# PM.TABLE.A4$pm.glm.A4.Chisq <- round(PM.TABLE.A4$pm.glm.A4.Chisq, 2)
# PM.TABLE.A4$pm.glm.A4.Pr <- format.pval(PM.TABLE.A4$pm.glm.A4.Pr, digits=2, eps= 0.001)
# PM.TABLE.A4$pm.glm.A4.Pr <- gsub("0\\.", ".", PM.TABLE.A4$pm.glm.A4.Pr)
# rownames(PM.TABLE.A4) <- c("Stimulus","TP","Stimulus*TP")
# PM.TABLE.A4


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
CDTC.RT.lmer.A4 <- lmer(y ~ S*block*cond+(1|s), data=dCDT.RT)
save(CDTC.RT.lmer.A4, file="CDTC.RT.lmer.A4.RData")
load("CDTC.RT.lmer.A4.RData")
CDT.RT.glm.A4 <- Anova(CDTC.RT.lmer.A4,type="II")
CDT.RT.glm.A4
# CDT.RT.TABLE.A4 <- data.frame(CDT.RT.glm.A4$Chisq, CDT.RT.glm.A4$Df, CDT.RT.glm.A4$Pr)
# CDT.RT.TABLE.A4$CDT.RT.glm.A4.Chisq <- round(CDT.RT.TABLE.A4$CDT.RT.glm.A4.Chisq, 2)
# CDT.RT.TABLE.A4$CDT.RT.glm.A4.Pr <- format.pval(CDT.RT.TABLE.A4$CDT.RT.glm.A4.Pr, digits=2, eps= 0.001)
# CDT.RT.TABLE.A4$CDT.RT.glm.A4.Pr <- gsub("0\\.", ".", CDT.RT.TABLE.A4$CDT.RT.glm.A4.Pr)
# rownames(CDT.RT.TABLE.A4) <- c("Stimulus","PM Block","TP","Stimulus*Block","Stimulus*TP","Block*TP","Stimulus*Block*TP")
# CDT.RT.TABLE.A4

#
PMC.RT.lmer.A4 <- lmer(y ~ S*cond+(1|s), data=dPMT.RT)
save(PMC.RT.lmer.A4, file="PMC.RT.lmer.A4.RData")
load("PMC.RT.lmer.A4.RData")
PM.RT.glm.A4 <- Anova(PMC.RT.lmer.A4,type="II")
PM.RT.glm.A4
# PM.RT.TABLE.A4 <- data.frame(PM.RT.glm.A4$Chisq, PM.RT.glm.A4$Df, PM.RT.glm.A4$Pr)
# PM.RT.TABLE.A4$PM.RT.glm.A4.Chisq <- round(PM.RT.TABLE.A4$PM.RT.glm.A4.Chisq, 2)
# PM.RT.TABLE.A4$PM.RT.glm.A4.Pr <- format.pval(PM.RT.TABLE.A4$PM.RT.glm.A4.Pr, digits=2, eps= 0.001)
# PM.RT.TABLE.A4$PM.RT.glm.A4.Pr <- gsub("0\\.", ".", PM.RT.TABLE.A4$PM.RT.glm.A4.Pr)
# rownames(PM.RT.TABLE.A4) <- c("Stimulus","TP","Stimulus*TP")
# PM.RT.TABLE.A4


# # # Manifests - Mean Reponse Proportion # # #
#

# length(CDT$RT[CDT$S=="cc" & CDT$R=="C"])/
#     length(CDT$RT[CDT$S=="cc"])  # CORRECT ACCURACY
#
# mean(CDT$C[CDT$S=="cc"])  # CORRECT ACCURACY
#

levels(CDT$S) <- c("Conflict","Nonconflict")
levels(CDT$block) <- c("Control","PM")
levels(CDT$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
str(CDT)

levels(PMT$S) <- c("PM (Conflict)","PM (Nonconflict)")
levels(PMT$block) <- c("PM")
levels(PMT$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
str(PMT)


# Ongoing Task Accuracy
CDT.Acc <- ddply(CDT, ~s*cond*block*S, summarise, RP=mean(C))
CDT.Acc.S <- ddply(CDT.Acc, ~S, summarise, M=mean(RP), SD=sd(RP))
CDT.Acc.BLOCK <- ddply(CDT.Acc, ~S*block, summarise, M=mean(RP), SD=sd(RP))
CDT.Acc.COND <- ddply(CDT.Acc, ~S*cond, summarise, M=mean(RP), SD=sd(RP))

# PM Task Accuracy
PMT.Acc <- ddply(PMT, ~s*cond*block*S, summarise, RP=mean(C))
PMT.Acc.S <- ddply(PMT.Acc, ~block, summarise, M=mean(RP), SD=sd(RP))
PMT.Acc.COND <- ddply(PMT.Acc, ~cond, summarise, M=mean(RP), SD=sd(RP))

names(CDT.Acc.S) <- c("Stimulus","Mean","SD")
names(PMT.Acc.S) <- c("PM Block","Mean","SD")
names(CDT.Acc.BLOCK) <- c("Stimulus","PM Block","Mean","SD")
names(CDT.Acc.COND) <- c("Stimulus","Time Pressure","Mean","SD")
names(PMT.Acc.COND) <- c("Time Pressure","Mean","SD")

CDT.Acc.S
PMT.Acc.S
CDT.Acc.BLOCK
CDT.Acc.COND
PMT.Acc.COND

# # # Manifests - Mean RT # # #
#

# Ongoing Task Correct RTs by Block
CDT.RT.corr.BLOCK <- ddply(data.A4[ (data.A4$S=="cc" & data.A4$C=="1") |
                                        (data.A4$S=="nn" & data.A4$C=="1"), ],
                           ~S*block, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))

# PM Task Correct RTs by Block
PMT.RT.corr.BLOCK <- ddply(data.A4[ (data.A4$S=="pc" & data.A4$C==1) |
                                        (data.A4$S=="pn" & data.A4$C==1), ],
                           ~block, summarise, M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))

# Ongoing Task Correct RTs by Condition
CDT.RT.corr.COND <- ddply(data.A4[ (data.A4$S=="cc" & data.A4$C==1) | (data.A4$S=="nn" & data.A4$C==1), ], ~S*cond, summarise,
      M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))

# PM Task Correct RTs by Condition
PMT.RT.corr.COND <- ddply(data.A4[ (data.A4$S=="pc" & data.A4$C==1) | (data.A4$S=="pn" & data.A4$C==1), ], ~cond, summarise,
                     M=mean(RT, na.rm = TRUE), SD=sd(RT, na.rm = TRUE))

levels(CDT.RT.corr.BLOCK$S) <- c("Conflict","Nonconflict","PM (Conflict)","PM (Nonconflict)")
levels(CDT.RT.corr.BLOCK$block) <- c("Control","PM")
names(CDT.RT.corr.BLOCK) <- c("Stimulus","PM Block","Mean","SD")
CDT.RT.corr.BLOCK.TABLE <- xtable(CDT.RT.corr.BLOCK, auto = TRUE, digits = 2, caption = "Ongoing Task RT (s) by PM Block")
CDT.RT.corr.BLOCK.TABLE


levels(PMT.RT.corr.BLOCK$block) <- c("Control","PM")
names(PMT.RT.corr.BLOCK) <- c("PM Block","Mean","SD")
PMT.RT.corr.BLOCK.TABLE <- xtable(PMT.RT.corr.BLOCK, auto = TRUE, digits = 2, caption = "Overall PM RT (s)")
PMT.RT.corr.BLOCK.TABLE


levels(CDT.RT.corr.COND$S) <- c("Conflict","Nonconflict","PM (Conflict)","PM (Nonconflict)")
levels(CDT.RT.corr.COND$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
names(CDT.RT.corr.COND) <- c("Stimulus","Time Pressure","Mean","SD")
CDT.RT.corr.COND.TABLE <- xtable(CDT.RT.corr.COND, auto = TRUE, digits = 2, caption = "Ongoing Task RT (s) by Time Pressure")
CDT.RT.corr.COND.TABLE

levels(PMT.RT.corr.COND$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
names(PMT.RT.corr.COND) <- c("Time Pressure","Mean","SD")
PMT.RT.corr.COND.TABLE <- xtable(PMT.RT.corr.COND, auto = TRUE, digits = 2, caption = "PM RT (s) by Time Pressure")
PMT.RT.corr.COND.TABLE


CDT.Acc.S.TABLE <- xtable(CDT.Acc.S, auto = TRUE, digits = 2, caption = "Ongoing Task Accuracy")
CDT.Acc.S.TABLE

PMT.Acc.S.TABLE <- xtable(PMT.Acc.S, auto = TRUE, digits = 2, caption = "Overall PM Accuracy")
PMT.Acc.S.TABLE

CDT.Acc.BLOCK.TABLE <- xtable(CDT.Acc.BLOCK, auto = TRUE, digits = 2, caption = "Ongoing Task Accuracy by PM Block")
CDT.Acc.BLOCK.TABLE

CDT.Acc.COND.TABLE <- xtable(CDT.Acc.COND, auto = TRUE, digits = 2, caption = "Ongoing Task Accuracy by Time Pressure")
CDT.Acc.COND.TABLE

PMT.Acc.COND.TABLE <- xtable(PMT.Acc.COND, auto = TRUE, digits = 2, caption = "PM Accuracy by Time Pressure")
PMT.Acc.COND.TABLE







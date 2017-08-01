
rm(list=ls())
# setwd("~/russ_model_analyses")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")
load_model ("LBA","lbaN_B.R")
source("LSAnova.R")
require("lsr")
require("lme4")
require("car")
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
# install.packages(pkgs) #install
sapply(pkgs, require, character.only = T)

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.E1 <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)


# names(samples.E1) != "p17"
samples.E1 <- samples.E1[names(samples.E1) != "p17"]  # Exclude p17 E1 due to no PM responses

# # # Get data from samples object # # #
#
get.hdata.dmc <- function(hsamples){
  list.wind<-lapply(seq_along(hsamples), function(samples, n, i) cbind(n[[i]], samples[[i]]$data),
                    samples= hsamples, n = names(hsamples))
  out<-do.call(rbind, list.wind)
  names(out)[1] <- "s"
  out
}

# Get data from samples object
data.E1 <- get.hdata.dmc(samples.E1)
rm(samples.E1)

# head(data.E1)
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






# # # Prep dataframes for analysis # # #
#
# Conflict Detection Task trials only (no PM)
CDT <- data.E1[!(data.E1$S=="pc" | data.E1$S=="pn"),]
CDT$S <- factor(as.character(CDT$S)); CDT$R <- factor(as.character(CDT$R))
# str(CDT)
# head(CDT)

# PM Task trials only
PMT <- data.E1[(data.E1$S=="pc" | data.E1$S=="pn"),]
PMT$S <- factor(as.character(PMT$S)); PMT$block <- factor(as.character(PMT$block))
# str(PMT)
# head(PMT)
#

# CDT.Acc.glmer.E1 <- glmer(C ~ S*block*cond+(1|s), data=CDT, family=binomial(link="probit"))
# save(CDT.Acc.glmer.E1, file="analysis/CDT.Acc.glmer.E1.RData")
load("analysis/CDT.Acc.glmer.E1.RData")
CDT.Acc.glm.E1 <- Anova(CDT.Acc.glmer.E1,type="II")
# CDT.Acc.glm.E1
CDT.Acc.TABLE.E1 <- data.frame(CDT.Acc.glm.E1$Chisq, CDT.Acc.glm.E1$Df, CDT.Acc.glm.E1$Pr)
# CDT.Acc.TABLE.E1
CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Chisq <- round(CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Chisq, 2)
CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Pr <- format.pval(CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Pr, digits=2, eps= 0.001)
CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Pr <- gsub("0\\.", ".", CDT.Acc.TABLE.E1$CDT.Acc.glm.E1.Pr)
rownames(CDT.Acc.TABLE.E1) <- c("Stimulus","PM Block","Time Pressure",
                                "Stimulus by PM Block","Stimulus by Time Pressure","PM Block by Time Pressure",
                                "Stimulus by PM Block by Time Pressure")
colnames(CDT.Acc.TABLE.E1) <- c("Chi-square","df","p")
# CDT.Acc.TABLE.E1
pander(CDT.Acc.TABLE.E1)
write.csv(CDT.Acc.TABLE.E1, file="analysis/LMER.CDT.Acc.E1.csv")

# PMT.Acc.glmer.E1 <- glmer(C ~ S*cond+(1|s), data=PMT, family=binomial(link="probit"))
# save(PMT.Acc.glmer.E1, file="analysis/PMT.Acc.glmer.E1.RData")
load("analysis/PMT.Acc.glmer.E1.RData")
PMT.Acc.glm.E1 <- Anova(PMT.Acc.glmer.E1,type="II")
# PMT.Acc.glm.E1
PMT.Acc.TABLE.E1 <- data.frame(PMT.Acc.glm.E1$Chisq, PMT.Acc.glm.E1$Df, PMT.Acc.glm.E1$Pr)
# PMT.Acc.TABLE.E1
PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Chisq <- round(PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Chisq, 2)
PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Pr <- format.pval(PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Pr, digits=2, eps= 0.001)
PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Pr <- gsub("0\\.", ".", PMT.Acc.TABLE.E1$PMT.Acc.glm.E1.Pr)
rownames(PMT.Acc.TABLE.E1) <- c("Stimulus","Time Pressure",
                                "Stimulus by Time Pressure")
colnames(PMT.Acc.TABLE.E1) <- c("Chi-square","df","p")
# PMT.Acc.TABLE.E1
pander(PMT.Acc.TABLE.E1)
write.csv(PMT.Acc.TABLE.E1, file="analysis/LMER.PMT.Acc.E1.csv")




# # # Prep RT dataframes for analysis (keep correct RTs only) # # #
#
# Conflict Detection Task correct response trials only (no PM)
CDT.RT <- arr2df(tapply(CDT$RT[CDT$C==1],
                         list(s=CDT$s[CDT$C==1], block=CDT$block[CDT$C==1],
                              cond=CDT$cond[CDT$C==1], S=CDT$S[CDT$C==1]), mean))
# head(CDT.RT)
# str(CDT.RT)

# PM Task correct response trials only
PMT.RT <- arr2df(tapply(PMT$RT[PMT$C==1],
                        list(s=PMT$s[PMT$C==1], cond=PMT$cond[PMT$C==1], S=PMT$S[PMT$C==1]), mean))
# head(PMT.RT)
# str(PMT.RT)


# CDTC.RT.lmer.E1 <- lmer(y ~ S*block*cond+(1|s), data=CDT.RT)
# save(CDTC.RT.lmer.E1, file="analysis/CDTC.RT.lmer.E1.RData")
load("analysis/CDTC.RT.lmer.E1.RData")
CDT.RT.glm.E1 <- Anova(CDTC.RT.lmer.E1,type="II")
# CDT.RT.glm.E1
CDT.RT.TABLE.E1 <- data.frame(CDT.RT.glm.E1$Chisq, CDT.RT.glm.E1$Df, CDT.RT.glm.E1$Pr)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq <- round(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Chisq, 2)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- format.pval(CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr, digits=2, eps= 0.001)
CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr <- gsub("0\\.", ".", CDT.RT.TABLE.E1$CDT.RT.glm.E1.Pr)
rownames(CDT.RT.TABLE.E1) <- c("Stimulus","PM Block","Time Pressure",
                               "Stimulus by PM Block","Stimulus by Time Pressure","PM Block by Time Pressure",
                               "Stimulus by PM Block by Time Pressure")
colnames(CDT.RT.TABLE.E1) <- c("Chi-square","df","p")
# CDT.RT.TABLE.E1
pander(CDT.RT.TABLE.E1)
write.csv(CDT.RT.TABLE.E1, file="analysis/LMER.CDT.RT.E1.csv")

# PMTC.RT.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=PMT.RT)
# save(PMTC.RT.lmer.E1, file="analysis/PMTC.RT.lmer.E1.RData")
load("analysis/PMTC.RT.lmer.E1.RData")
PMT.RT.glm.E1 <- Anova(PMTC.RT.lmer.E1,type="II")
# PMT.RT.glm.E1
PMT.RT.TABLE.E1 <- data.frame(PMT.RT.glm.E1$Chisq, PMT.RT.glm.E1$Df, PMT.RT.glm.E1$Pr)
PMT.RT.TABLE.E1$PMT.RT.glm.E1.Chisq <- round(PMT.RT.TABLE.E1$PMT.RT.glm.E1.Chisq, 2)
PMT.RT.TABLE.E1$PMT.RT.glm.E1.Pr <- format.pval(PMT.RT.TABLE.E1$PMT.RT.glm.E1.Pr, digits=2, eps= 0.001)
PMT.RT.TABLE.E1$PMT.RT.glm.E1.Pr <- gsub("0\\.", ".", PMT.RT.TABLE.E1$PMT.RT.glm.E1.Pr)
rownames(PMT.RT.TABLE.E1) <- c("Stimulus","Time Pressure",
                                "Stimulus by Time Pressure")
colnames(PMT.RT.TABLE.E1) <- c("Chi-square","df","p")
# PMT.RT.TABLE.E1
pander(PMT.RT.TABLE.E1)
write.csv(PMT.RT.TABLE.E1, file="analysis/LMER.PMT.RT.E1.csv")


# # # ANOVA analysis: Ongoing Task Accuracy # # #
#
#
# Ongoing Task Accuracy Object
CDT.Acc <- arr2df(tapply(CDT$C,list(s=CDT$s, block=CDT$block, cond=CDT$cond, S=CDT$S),mean))

# str(CDT.Acc)
levels(CDT.Acc$block) <- c("Control","PM")
# levels(CDT.Acc$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.Acc$S) <- c("Conflict","Nonconflict")
# wsAnova(CDT.Acc)
# Manifests: Ongoing Task Accuracy
CDT.Acc.Means <- mneffects(CDT.Acc,list("S","block","cond",c("S","block"),c("S","cond"),c("block","cond"), digits=3))
CDT.Acc.S <- data.frame(rbind(round(CDT.Acc.Means$S,3), round(se(CDT.Acc, facs=c("S")),3)))
rownames(CDT.Acc.S)<- c("M","SE")
# CDT.Acc.S

CDT.Acc.Block <- data.frame(rbind(round(CDT.Acc.Means$block,3), round(se(CDT.Acc, facs=c("block")),3)))
rownames(CDT.Acc.Block)<- c("M","SE")
# CDT.Acc.Block

CDT.Acc.Cond <- data.frame(rbind(round(CDT.Acc.Means$cond,3), round(se(CDT.Acc, facs=c("cond")),3)))
rownames(CDT.Acc.Cond)<- c("M","SE")
# CDT.Acc.Cond

# round(se(CDT.Acc, facs=c("S","block")),3)
# round(se(CDT.Acc, facs=c("S","cond")),3)
# round(se(CDT.Acc, facs=c("block","cond")),3)

# Split up data object for cond comparisons
CDT.Acc$S <- as.numeric(factor(CDT.Acc$S))
bonf.CDT.Acc <- data.frame()
bonf.CDT.Acc <- tapply(CDT.Acc$y,
                                list(s=CDT.Acc$s,
                                     cond=factor(CDT.Acc$cond)), mean)
# bonf.CDT.Acc
tAB <- t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,2], paired=T)
dAB <- cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,2], method="paired")
AB <- data.frame(cbind(t=tAB$statistic, df=tAB$parameter, p=tAB$p.value, Difference=tAB$estimate, d=dAB))

tBC <- t.test (bonf.CDT.Acc[,2], bonf.CDT.Acc[,3], paired=T)
dBC <- cohensD(x=bonf.CDT.Acc[,2], y=bonf.CDT.Acc[,3], method="paired")
BC <- data.frame(cbind(t=tBC$statistic, df=tBC$parameter, p=tBC$p.value, Difference=tBC$estimate, d=dBC))

tCD <- t.test (bonf.CDT.Acc[,3], bonf.CDT.Acc[,4], paired=T)
dCD <- cohensD(x=bonf.CDT.Acc[,3], y=bonf.CDT.Acc[,4], method="paired")
CD <- data.frame(cbind(t=tCD$statistic, df=tCD$parameter, p=tCD$p.value, Difference=tCD$estimate, d=dCD))

# tAC <- t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,3], paired=T)
# dAC <- cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,3], method="paired")
# tBD <- t.test (bonf.CDT.Acc[,2], bonf.CDT.Acc[,4], paired=T)
# dBD <- cohensD(x=bonf.CDT.Acc[,2], y=bonf.CDT.Acc[,4], method="paired")
# tAD <- t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,4], paired=T)
# dAD <- cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,4], method="paired")

T.Contrasts.CDT.Acc <- data.frame(rbind(AB, BC, CD))
T.Contrasts.CDT.Acc$t <- round(T.Contrasts.CDT.Acc$t, digits=2)
T.Contrasts.CDT.Acc$Difference <- round(T.Contrasts.CDT.Acc$Difference*100, digits=2)
T.Contrasts.CDT.Acc$d <- round(T.Contrasts.CDT.Acc$d, digits=2)
T.Contrasts.CDT.Acc$p <- round(T.Contrasts.CDT.Acc$p, digits=3)
T.Contrasts.CDT.Acc$p <- format.pval(T.Contrasts.CDT.Acc$p, digits=2, eps= 0.001)
row.names(T.Contrasts.CDT.Acc) <- c("A-B","B-C","C-D")
# T.Contrasts.CDT.Acc
pander(T.Contrasts.CDT.Acc)
write.csv(T.Contrasts.CDT.Acc, file="analysis/T.Contrasts.CDT.Acc.E1.csv")


# # # ANOVA analysis: Ongoing Task RT # # #
#
#
# Ongoing Task RT Object
CDT.RT <- arr2df(tapply(CDT$RT[CDT$C==1],
                        list(s=CDT$s[CDT$C==1], block=CDT$block[CDT$C==1],
                             cond=CDT$cond[CDT$C==1], S=CDT$S[CDT$C==1]), mean))
# str(CDT.RT)
levels(CDT.RT$block) <- c("Control","PM")
# levels(CDT.RT$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.RT$S) <- c("Conflict","Nonconflict")
# wsAnova(CDT.RT)
# Manifests: Ongoing Task RT
CDT.RT.Means <- mneffects(CDT.RT,list("S","block","cond",c("S","block"),c("S","cond"), digits=3))

CDT.RT.S <- data.frame(rbind(round(CDT.RT.Means$S,3), round(se(CDT.RT, facs=c("S")),3)))
rownames(CDT.RT.S)<- c("M","SE")
# CDT.RT.S

CDT.RT.Block <- data.frame(rbind(round(CDT.RT.Means$block,3), round(se(CDT.RT, facs=c("block")),3)))
rownames(CDT.RT.Block)<- c("M","SE")
# CDT.RT.Block

CDT.RT.Cond <- data.frame(rbind(round(CDT.RT.Means$cond,3), round(se(CDT.RT, facs=c("cond")),3)))
rownames(CDT.RT.Cond)<- c("M","SE")
# CDT.RT.Cond

round(se(CDT.RT, facs=c("S","cond")),3)
# Split up data object for cond comparisons
CDT.RT$S <- as.numeric(factor(CDT.RT$S))
bonf.CDT.RT <- data.frame()
bonf.CDT.RT <- tapply(CDT.RT$y,
                       list(s=CDT.RT$s,
                            cond=factor(CDT.RT$cond)), mean)
# bonf.CDT.RT
tAB <- t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,2], paired=T)
dAB <- cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,2], method="paired")
AB <- data.frame(cbind(t=tAB$statistic, df=tAB$parameter, p=tAB$p.value, Difference=tAB$estimate, d=dAB))

tBC <- t.test (bonf.CDT.RT[,2], bonf.CDT.RT[,3], paired=T)
dBC <- cohensD(x=bonf.CDT.RT[,2], y=bonf.CDT.RT[,3], method="paired")
BC <- data.frame(cbind(t=tBC$statistic, df=tBC$parameter, p=tBC$p.value, Difference=tBC$estimate, d=dBC))

tCD <- t.test (bonf.CDT.RT[,3], bonf.CDT.RT[,4], paired=T)
dCD <- cohensD(x=bonf.CDT.RT[,3], y=bonf.CDT.RT[,4], method="paired")
CD <- data.frame(cbind(t=tCD$statistic, df=tCD$parameter, p=tCD$p.value, Difference=tCD$estimate, d=dCD))

# tAC <- t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,3], paired=T)
# dAC <- cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,3], method="paired")
# tBD <- t.test (bonf.CDT.RT[,2], bonf.CDT.RT[,4], paired=T)
# dBD <- cohensD(x=bonf.CDT.RT[,2], y=bonf.CDT.RT[,4], method="paired")
# tAD <- t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,4], paired=T)
# dAD <- cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,4], method="paired")

T.Contrasts.CDT.RT <- data.frame(rbind(AB, BC, CD))
T.Contrasts.CDT.RT$t <- round(T.Contrasts.CDT.RT$t, digits=2)
T.Contrasts.CDT.RT$Difference <- round(T.Contrasts.CDT.RT$Difference, digits=2)
T.Contrasts.CDT.RT$d <- round(T.Contrasts.CDT.RT$d, digits=2)
T.Contrasts.CDT.RT$p <- round(T.Contrasts.CDT.RT$p, digits=3)
T.Contrasts.CDT.RT$p <- format.pval(T.Contrasts.CDT.RT$p, digits=2, eps= 0.001)
row.names(T.Contrasts.CDT.RT) <- c("A-B","B-C","C-D")
# T.Contrasts.CDT.RT
pander(T.Contrasts.CDT.RT)
write.csv(T.Contrasts.CDT.RT, file="analysis/T.Contrasts.CDT.RT.E1.csv")


# # # ANOVA analysis: PM Accuracy # # #
#
# PM Accuracy Object
PMT.Acc <- arr2df(tapply(PMT$C,
                         list(s=PMT$s, cond=PMT$cond, S=PMT$S),mean))
# str(PMT.Acc)
# levels(PMT.Acc$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(PMT.Acc$S) <- c("PM_Conflict","PM_Nonconflict")
# wsAnova(PMT.Acc)
PMT.Acc.Means <- mneffects(PMT.Acc,list("S","cond",c("S","cond"), digits=3))

PMT.Acc.S <- data.frame(rbind(round(PMT.Acc.Means$S,3), round(se(PMT.Acc, facs=c("S")),3)))
rownames(PMT.Acc.S)<- c("M","SE")
# PMT.Acc.S

PMT.Acc.Cond <- data.frame(rbind(round(PMT.Acc.Means$cond,3), round(se(PMT.Acc, facs=c("cond")),3)))
rownames(PMT.Acc.Cond)<- c("M","SE")
# PMT.Acc.Cond

# Split up data object for cond comparisons
PMT.Acc$S <- as.numeric(factor(PMT.Acc$S))
bonf.PMT.Acc <- data.frame()
bonf.PMT.Acc <- tapply(PMT.Acc$y,
                       list(s=PMT.Acc$s,
                            cond=factor(PMT.Acc$cond)), mean)
# bonf.PMT.Acc
tAB <- t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,2], paired=T)
dAB <- cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,2], method="paired")
AB <- data.frame(cbind(t=tAB$statistic, df=tAB$parameter, p=tAB$p.value, Difference=tAB$estimate, d=dAB))

tBC <- t.test (bonf.PMT.Acc[,2], bonf.PMT.Acc[,3], paired=T)
dBC <- cohensD(x=bonf.PMT.Acc[,2], y=bonf.PMT.Acc[,3], method="paired")
BC <- data.frame(cbind(t=tBC$statistic, df=tBC$parameter, p=tBC$p.value, Difference=tBC$estimate, d=dBC))

tCD <- t.test (bonf.PMT.Acc[,3], bonf.PMT.Acc[,4], paired=T)
dCD <- cohensD(x=bonf.PMT.Acc[,3], y=bonf.PMT.Acc[,4], method="paired")
CD <- data.frame(cbind(t=tCD$statistic, df=tCD$parameter, p=tCD$p.value, Difference=tCD$estimate, d=dCD))

# tAC <- t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,3], paired=T)
# dAC <- cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,3], method="paired")
# tBD <- t.test (bonf.PMT.Acc[,2], bonf.PMT.Acc[,4], paired=T)
# dBD <- cohensD(x=bonf.PMT.Acc[,2], y=bonf.PMT.Acc[,4], method="paired")
# tAD <- t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,4], paired=T)
# dAD <- cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,4], method="paired")

T.Contrasts.PMT.Acc <- data.frame(rbind(AB, BC, CD))
T.Contrasts.PMT.Acc$t <- round(T.Contrasts.PMT.Acc$t, digits=2)
T.Contrasts.PMT.Acc$Difference <- round(T.Contrasts.PMT.Acc$Difference*100, digits=2)
T.Contrasts.PMT.Acc$d <- round(T.Contrasts.PMT.Acc$d, digits=2)
T.Contrasts.PMT.Acc$p <- round(T.Contrasts.PMT.Acc$p, digits=3)
T.Contrasts.PMT.Acc$p <- format.pval(T.Contrasts.PMT.Acc$p, digits=2, eps= 0.001)
row.names(T.Contrasts.PMT.Acc) <- c("A-B","B-C","C-D")
# T.Contrasts.PMT.Acc
pander(T.Contrasts.PMT.Acc)
write.csv(T.Contrasts.PMT.Acc, file="analysis/T.Contrasts.PMT.Acc.E1.csv")


# # # ANOVA analysis: PM RT # # #
#
# PM RT Object
PMT.RT <- arr2df(tapply(PMT$RT[PMT$C==1],
                        list(s=PMT$s[PMT$C==1], cond=PMT$cond[PMT$C==1], S=PMT$S[PMT$C==1]), mean))
# str(PMT.RT)
# levels(PMT.RT$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(PMT.RT$S) <- c("PM_Conflict","PM_Nonconflict")
# wsAnova(PMT.RT)
PMT.RT.Means <- mneffects(PMT.RT,list("S","cond",c("S","cond"), digits=3))

PMT.RT.S <- data.frame(rbind(round(PMT.RT.Means$S,3), round(se2(PMT.RT, facs=c("S")),3)))
rownames(PMT.RT.S)<- c("M","SE")
# PMT.RT.S

PMT.RT.Cond <- data.frame(rbind(round(PMT.RT.Means$cond,3), round(se2(PMT.RT, facs=c("cond")),3)))
rownames(PMT.RT.Cond)<- c("M","SE")
# PMT.RT.Cond

# Split up data object for cond comparisons
PMT.RT$S <- as.numeric(factor(PMT.RT$S))
bonf.PMT.RT <- data.frame()
bonf.PMT.RT <- tapply(PMT.RT$y,
                      list(s=PMT.RT$s,
                           cond=factor(PMT.RT$cond)), mean)
# bonf.PMT.RT
tAB <- t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,2], paired=T)
dAB <- cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,2], method="paired")
AB <- data.frame(cbind(t=tAB$statistic, df=tAB$parameter, p=tAB$p.value, Difference=tAB$estimate, d=dAB))

tBC <- t.test (bonf.PMT.RT[,2], bonf.PMT.RT[,3], paired=T)
dBC <- cohensD(x=bonf.PMT.RT[,2], y=bonf.PMT.RT[,3], method="paired")
BC <- data.frame(cbind(t=tBC$statistic, df=tBC$parameter, p=tBC$p.value, Difference=tBC$estimate, d=dBC))

tCD <- t.test (bonf.PMT.RT[,3], bonf.PMT.RT[,4], paired=T)
dCD <- cohensD(x=bonf.PMT.RT[,3], y=bonf.PMT.RT[,4], method="paired")
CD <- data.frame(cbind(t=tCD$statistic, df=tCD$parameter, p=tCD$p.value, Difference=tCD$estimate, d=dCD))

# tAC <- t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,3], paired=T)
# dAC <- cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,3], method="paired")
# tBD <- t.test (bonf.PMT.RT[,2], bonf.PMT.RT[,4], paired=T)
# dBD <- cohensD(x=bonf.PMT.RT[,2], y=bonf.PMT.RT[,4], method="paired")
# tAD <- t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,4], paired=T)
# dAD <- cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,4], method="paired")

T.Contrasts.PMT.RT <- data.frame(rbind(AB, BC, CD))
T.Contrasts.PMT.RT$t <- round(T.Contrasts.PMT.RT$t, digits=2)
T.Contrasts.PMT.RT$Difference <- round(T.Contrasts.PMT.RT$Difference, digits=2)
T.Contrasts.PMT.RT$d <- round(T.Contrasts.PMT.RT$d, digits=2)
T.Contrasts.PMT.RT$p <- round(T.Contrasts.PMT.RT$p, digits=3)
T.Contrasts.PMT.RT$p <- format.pval(T.Contrasts.PMT.RT$p, digits=2, eps= 0.001)
row.names(T.Contrasts.PMT.RT) <- c("A-B","B-C","C-D")
# T.Contrasts.PMT.RT
pander(T.Contrasts.PMT.RT)
write.csv(T.Contrasts.PMT.RT, file="analysis/T.Contrasts.PMT.RT.E1.csv")



# # # ANOVA analysis: Ongoing RT for PM versus non-PM Trials # # #
#
# Ongoing Task RT Object
# str(data.E1)
# any((data.E1$block=="3" & data.E1$S=="pc" & data.E1$R=="C"))
# any((data.E1$block=="3" & data.E1$S=="pn" & data.E1$R=="N"))

CDT.Reactive <- data.E1[(data.E1$block=="3" & ((data.E1$S=="cc" & data.E1$R=="C") | (data.E1$S=="pc" & data.E1$R=="C") |
                                                   (data.E1$S=="nn" & data.E1$R=="N") | (data.E1$S=="pn" & data.E1$R=="N"))),]
CDT.Reactive$block <- factor(CDT.Reactive$block)
CDT.Reactive$R <- factor(CDT.Reactive$R)
# any((CDT.Reactive$block=="3" & CDT.Reactive$S=="pc" & CDT.Reactive$R=="C"))
# any((CDT.Reactive$block=="3" & CDT.Reactive$S=="pn" & CDT.Reactive$R=="N"))
# any(CDT.Reactive$R=="P")
# any(CDT.Reactive$block=="2")
# any(CDT.Reactive$S=="pc")
# str(CDT.Reactive)
# levels(CDT.Reactive$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.Reactive$S) <- c("Conflict","Nonconflict","PM_Conflict","PM_Nonconflict")
CDT.Reactive <- arr2df(tapply(CDT.Reactive$RT,
                        list(s=CDT.Reactive$s,
                             cond=CDT.Reactive$cond, S=CDT.Reactive$S), mean))
# str(CDT.Reactive)
# CDT.Reactive.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=CDT.Reactive)
# save(CDT.Reactive.lmer.E1, file="analysis/CDT.Reactive.lmer.E1.RData")
load("analysis/CDT.Reactive.lmer.E1.RData")
CDT.Reactive.glm.E1 <- Anova(CDT.Reactive.lmer.E1,type="II")
# CDT.Reactive.glm.E1
CDT.Reactive.TABLE.E1 <- data.frame(CDT.Reactive.glm.E1$Chisq, CDT.Reactive.glm.E1$Df, CDT.Reactive.glm.E1$Pr)
CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Chisq <- round(CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Chisq, 2)
CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Pr <- format.pval(CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Pr, digits=2, eps= 0.001)
CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Pr <- gsub("0\\.", ".", CDT.Reactive.TABLE.E1$CDT.Reactive.glm.E1.Pr)
rownames(CDT.Reactive.TABLE.E1) <- c("Stimulus","Time Pressure","Stimulus by Time Pressure")
colnames(CDT.Reactive.TABLE.E1) <- c("Chi-square","df","p")
pander(CDT.Reactive.TABLE.E1)
write.csv(CDT.Reactive.TABLE.E1, file="analysis/LMER.CDT.Reactive.E1.csv")


# wsAnova(CDT.Reactive)
CDT.Reactive.Means <- mneffects(CDT.Reactive,list("S","cond", digits=3))

CDT.Reactive.S <- data.frame(rbind(round(CDT.Reactive.Means$S,3), round(se2(CDT.Reactive, facs=c("S")),3)))
rownames(CDT.Reactive.S)<- c("M","SE")
CDT.Reactive.S

CDT.Reactive.Cond <- data.frame(rbind(round(CDT.Reactive.Means$cond,3), round(se2(CDT.Reactive, facs=c("cond")),3)))
rownames(CDT.Reactive.Cond)<- c("M","SE")
CDT.Reactive.Cond

# Split up data object for S comparisons
CDT.Reactive$S <- as.numeric(factor(CDT.Reactive$S))
bonf.CDT.Reactive <- data.frame()
bonf.CDT.Reactive <- tapply(CDT.Reactive$y,
                      list(s=CDT.Reactive$s,
                           S=factor(CDT.Reactive$S)), mean)
# bonf.CDT.Reactive
tccpc <- t.test (bonf.CDT.Reactive[,1], bonf.CDT.Reactive[,3], paired=T)
dccpc <- cohensD(x=bonf.CDT.Reactive[,1], y=bonf.CDT.Reactive[,3], method="paired")
ccpc <- data.frame(cbind(t=tccpc$statistic, df=tccpc$parameter, p=tccpc$p.value, Difference=tccpc$estimate, d=dccpc))

tnnpn <- t.test (bonf.CDT.Reactive[,2], bonf.CDT.Reactive[,4], paired=T)
dnnpn <- cohensD(x=bonf.CDT.Reactive[,2], y=bonf.CDT.Reactive[,4], method="paired")
nnpn <- data.frame(cbind(t=tnnpn$statistic, df=tnnpn$parameter, p=tnnpn$p.value, Difference=tnnpn$estimate, d=dnnpn))

T.Contrasts.Reactive <- data.frame(rbind(ccpc, nnpn))
T.Contrasts.Reactive$t <- round(T.Contrasts.Reactive$t, digits=2)
T.Contrasts.Reactive$Difference <- round(T.Contrasts.Reactive$Difference, digits=2)
T.Contrasts.Reactive$d <- round(T.Contrasts.Reactive$d, digits=2)
T.Contrasts.Reactive$p <- round(T.Contrasts.Reactive$p, digits=3)
T.Contrasts.Reactive$p <- format.pval(T.Contrasts.Reactive$p, digits=2, eps= 0.001)
row.names(T.Contrasts.Reactive) <- c("CR.minus.CR-to-PMC","NR.minus.NR-to-PMN")
# T.Contrasts.Reactive
pander(T.Contrasts.Reactive)
write.csv(T.Contrasts.Reactive, file="analysis/T.Contrasts.Reactive.E1.csv")




# # # Save off Tables of Means and SDs # # #
#

## Ongoing Task Accuracy
means1 <- CDT.Acc.S
means2 <- CDT.Acc.Block
means3 <- CDT.Acc.Cond

pander(means1)
pander(means2)
pander(means3)

write.csv(means1, file="analysis/CDT.Acc.S.E1.csv")
write.csv(means2, file="analysis/CDT.Acc.Block.E1.csv")
write.csv(means3, file="analysis/CDT.Acc.Cond.E1.csv")


## Ongoing Task RT
means1 <- CDT.RT.S
means2 <- CDT.RT.Block
means3 <- CDT.RT.Cond

pander(means1)
pander(means2)
pander(means3)

write.csv(means1, file="analysis/CDT.RT.S.E1.csv")
write.csv(means2, file="analysis/CDT.RT.Block.E1.csv")
write.csv(means3, file="analysis/CDT.RT.Cond.E1.csv")


## PM Accuracy
means1 <- PMT.Acc.S
means2 <- PMT.Acc.Cond

pander(means1)
pander(means2)

write.csv(means1, file="analysis/PMT.Acc.S.E1.csv")
write.csv(means2, file="analysis/PMT.Acc.Cond.E1.csv")


## PM RT
means1 <- PMT.RT.S
means2 <- PMT.RT.Cond

pander(means1)
pander(means2)

write.csv(means1, file="analysis/PMT.RT.S.E1.csv")
write.csv(means2, file="analysis/PMT.RT.Cond.E1.csv")


## Reactive Control RT Check LMER

# Check that congruent ongoing task responses to PM stimuli (i.e., conflict responses to PM conflicts) are faster than correct ongoing task responses to non-PM stimuli (although this will be partly due to statistical facilitation from the fast PM accumulator on PM trials filtering out ongoing task responses).

means1 <- CDT.Reactive.S

pander(means1)

write.csv(means1, file="analysis/CDT.Reactive.S.E1.csv")



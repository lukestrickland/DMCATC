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
require(lsr)
require(lme4)
require(car)
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
# install.packages(pkgs) #install
sapply(pkgs, require, character.only = T) #load
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

# CDT.Acc.glmer.E1 <- glmer(C ~ S*block*cond+(1|s), data=CDT, family=binomial(link="probit"))
# save(CDT.Acc.glmer.E1, file="CDT.Acc.glmer.E1.RData")
load("CDT.Acc.glmer.E1.RData")
CDT.Acc.glm.E1 <- Anova(CDT.Acc.glmer.E1,type="II")
CDT.Acc.glm.E1


# PMT.Acc.glmer.E1 <- glmer(C ~ S*cond+(1|s), data=PMT, family=binomial(link="probit"))
# save(PMT.Acc.glmer.E1, file="PMT.Acc.glmer.E1.RData")
load("PMT.Acc.glmer.E1.RData")
PMT.Acc.glm.E1 <- Anova(PMT.Acc.glmer.E1,type="II")
PMT.Acc.glm.E1


# # # Prep RT dataframes for analysis (keep correct RTs only) # # #
#
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC/analysis")

# Conflict Detection Task correct response trials only (no PM)
CDT.RT <- arr2df(tapply(CDT$RT[CDT$C==1],
                         list(s=CDT$s[CDT$C==1], block=CDT$block[CDT$C==1],
                              cond=CDT$cond[CDT$C==1], S=CDT$S[CDT$C==1]), mean))
head(CDT.RT)
str(CDT.RT)

# PM Task correct response trials only
PMT.RT <- arr2df(tapply(PMT$RT[PMT$C==1],
                        list(s=PMT$s[PMT$C==1], cond=PMT$cond[PMT$C==1], S=PMT$S[PMT$C==1]), mean))
head(PMT.RT)
str(PMT.RT)


# CDTC.RT.lmer.E1 <- lmer(y ~ S*block*cond+(1|s), data=CDT.RT)
# save(CDTC.RT.lmer.E1, file="CDTC.RT.lmer.E1.RData")
load("CDTC.RT.lmer.E1.RData")
CDT.RT.glm.E1 <- Anova(CDTC.RT.lmer.E1,type="II")
CDT.RT.glm.E1


# PMC.RT.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=PMT.RT)
# save(PMC.RT.lmer.E1, file="PMC.RT.lmer.E1.RData")
load("PMC.RT.lmer.E1.RData")
PM.RT.glm.E1 <- Anova(PMC.RT.lmer.E1,type="II")
PM.RT.glm.E1


# # # ANOVA analysis: Ongoing Task Accuracy # # #
#
#
# Ongoing Task Accuracy Object
CDT.Acc <- arr2df(tapply(CDT$C,list(s=CDT$s, block=CDT$block, cond=CDT$cond, S=CDT$S),mean))

str(CDT.Acc)
levels(CDT.Acc$block) <- c("Control","PM")
levels(CDT.Acc$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.Acc$S) <- c("Conflict","Nonconflict")
wsAnova(CDT.Acc)
# Manifests: Ongoing Task Accuracy
mneffects(CDT.Acc,list("S","block","cond",c("S","block"),c("S","cond"),c("block","cond"), digits=3))
round(se(CDT.Acc, facs=c("S")),3)
round(se(CDT.Acc, facs=c("block")),3)
round(se(CDT.Acc, facs=c("cond")),3)
round(se(CDT.Acc, facs=c("S","block")),3)
round(se(CDT.Acc, facs=c("S","cond")),3)
round(se(CDT.Acc, facs=c("block","cond")),3)
# Split up data object for cond comparisons
CDT.Acc$S <- as.numeric(factor(CDT.Acc$S))
bonf.CDT.Acc <- data.frame()
bonf.CDT.Acc <- tapply(CDT.Acc$y,
                                list(s=CDT.Acc$s,
                                     cond=factor(CDT.Acc$cond)), mean)
# bonf.CDT.Acc
t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,2], paired=T)
cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,2], method="paired")
t.test (bonf.CDT.Acc[,2], bonf.CDT.Acc[,3], paired=T)
cohensD(x=bonf.CDT.Acc[,2], y=bonf.CDT.Acc[,3], method="paired")
t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,3], paired=T)
cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,3], method="paired")

t.test (bonf.CDT.Acc[,3], bonf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.CDT.Acc[,3], y=bonf.CDT.Acc[,4], method="paired")
t.test (bonf.CDT.Acc[,2], bonf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.CDT.Acc[,2], y=bonf.CDT.Acc[,4], method="paired")
t.test (bonf.CDT.Acc[,1], bonf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.CDT.Acc[,1], y=bonf.CDT.Acc[,4], method="paired")
# Split up data object for S*cond comparisons
CDT.Acc$S <- as.numeric(factor(CDT.Acc$S))
bonf.Conflict.CDT.Acc <- data.frame()
bonf.Conflict.CDT.Acc <- tapply(CDT.Acc$y[CDT.Acc$S==1],
                                list(s=CDT.Acc$s[CDT.Acc$S==1],
                                     cond=factor(CDT.Acc$cond[CDT.Acc$S==1])), mean)
# bonf.Conflict.CDT.Acc
t.test (bonf.Conflict.CDT.Acc[,1], bonf.Conflict.CDT.Acc[,2], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,1], y=bonf.Conflict.CDT.Acc[,2], method="paired")
t.test (bonf.Conflict.CDT.Acc[,2], bonf.Conflict.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,2], y=bonf.Conflict.CDT.Acc[,3], method="paired")
t.test (bonf.Conflict.CDT.Acc[,1], bonf.Conflict.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,1], y=bonf.Conflict.CDT.Acc[,3], method="paired")

t.test (bonf.Conflict.CDT.Acc[,3], bonf.Conflict.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,3], y=bonf.Conflict.CDT.Acc[,4], method="paired")
t.test (bonf.Conflict.CDT.Acc[,2], bonf.Conflict.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,2], y=bonf.Conflict.CDT.Acc[,4], method="paired")
t.test (bonf.Conflict.CDT.Acc[,1], bonf.Conflict.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.Acc[,1], y=bonf.Conflict.CDT.Acc[,4], method="paired")

CDT.Acc$S <- as.numeric(factor(CDT.Acc$S))
bonf.Nonconf.CDT.Acc <- data.frame()
bonf.Nonconf.CDT.Acc <- tapply(CDT.Acc$y[CDT.Acc$S==2],
                               list(s=CDT.Acc$s[CDT.Acc$S==2],
                                    cond=factor(CDT.Acc$cond[CDT.Acc$S==2])), mean)
# bonf.Nonconf.CDT.Acc
t.test (bonf.Nonconf.CDT.Acc[,1], bonf.Nonconf.CDT.Acc[,2], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,1], y=bonf.Nonconf.CDT.Acc[,2], method="paired")
t.test (bonf.Nonconf.CDT.Acc[,2], bonf.Nonconf.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,2], y=bonf.Nonconf.CDT.Acc[,3], method="paired")
t.test (bonf.Nonconf.CDT.Acc[,1], bonf.Nonconf.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,1], y=bonf.Nonconf.CDT.Acc[,3], method="paired")

t.test (bonf.Nonconf.CDT.Acc[,3], bonf.Nonconf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,3], y=bonf.Nonconf.CDT.Acc[,4], method="paired")
t.test (bonf.Nonconf.CDT.Acc[,2], bonf.Nonconf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,2], y=bonf.Nonconf.CDT.Acc[,4], method="paired")
t.test (bonf.Nonconf.CDT.Acc[,1], bonf.Nonconf.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.Acc[,1], y=bonf.Nonconf.CDT.Acc[,4], method="paired")

# Split up data object for block*cond comparisons
CDT.Acc$block <- as.numeric(factor(CDT.Acc$block))
bonf.Control.CDT.Acc <- data.frame()
bonf.Control.CDT.Acc <- tapply(CDT.Acc$y[CDT.Acc$block==1],
                                list(s=CDT.Acc$s[CDT.Acc$block==1],
                                     cond=factor(CDT.Acc$cond[CDT.Acc$block==1])), mean)
bonf.Control.CDT.Acc
t.test (bonf.Control.CDT.Acc[,1], bonf.Control.CDT.Acc[,2], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,1], y=bonf.Control.CDT.Acc[,2], method="paired")
t.test (bonf.Control.CDT.Acc[,2], bonf.Control.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,2], y=bonf.Control.CDT.Acc[,3], method="paired")
t.test (bonf.Control.CDT.Acc[,1], bonf.Control.CDT.Acc[,3], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,1], y=bonf.Control.CDT.Acc[,3], method="paired")

t.test (bonf.Control.CDT.Acc[,3], bonf.Control.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,3], y=bonf.Control.CDT.Acc[,4], method="paired")
t.test (bonf.Control.CDT.Acc[,2], bonf.Control.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,2], y=bonf.Control.CDT.Acc[,4], method="paired")
t.test (bonf.Control.CDT.Acc[,1], bonf.Control.CDT.Acc[,4], paired=T)
cohensD(x=bonf.Control.CDT.Acc[,1], y=bonf.Control.CDT.Acc[,4], method="paired")

CDT.Acc$block <- as.numeric(factor(CDT.Acc$block))
bonf.PM.CDT.Acc <- data.frame()
bonf.PM.CDT.Acc <- tapply(CDT.Acc$y[CDT.Acc$block==2],
                               list(s=CDT.Acc$s[CDT.Acc$block==2],
                                    cond=factor(CDT.Acc$cond[CDT.Acc$block==2])), mean)
bonf.PM.CDT.Acc
t.test (bonf.PM.CDT.Acc[,1], bonf.PM.CDT.Acc[,2], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,1], y=bonf.PM.CDT.Acc[,2], method="paired")
t.test (bonf.PM.CDT.Acc[,2], bonf.PM.CDT.Acc[,3], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,2], y=bonf.PM.CDT.Acc[,3], method="paired")
t.test (bonf.PM.CDT.Acc[,1], bonf.PM.CDT.Acc[,3], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,1], y=bonf.PM.CDT.Acc[,3], method="paired")

t.test (bonf.PM.CDT.Acc[,3], bonf.PM.CDT.Acc[,4], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,3], y=bonf.PM.CDT.Acc[,4], method="paired")
t.test (bonf.PM.CDT.Acc[,2], bonf.PM.CDT.Acc[,4], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,2], y=bonf.PM.CDT.Acc[,4], method="paired")
t.test (bonf.PM.CDT.Acc[,1], bonf.PM.CDT.Acc[,4], paired=T)
cohensD(x=bonf.PM.CDT.Acc[,1], y=bonf.PM.CDT.Acc[,4], method="paired")

# # # ANOVA analysis: Ongoing Task RT # # #
#
#
# Ongoing Task RT Object
CDT.RT <- arr2df(tapply(CDT$RT[CDT$C==1],
                        list(s=CDT$s[CDT$C==1], block=CDT$block[CDT$C==1],
                             cond=CDT$cond[CDT$C==1], S=CDT$S[CDT$C==1]), mean))
str(CDT.RT)
levels(CDT.RT$block) <- c("Control","PM")
levels(CDT.RT$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.RT$S) <- c("Conflict","Nonconflict")
wsAnova(CDT.RT)
# Manifests: Ongoing Task RT
mneffects(CDT.RT,list("S","block","cond",c("S","block"),c("S","cond"), digits=3))
round(se(CDT.RT, facs=c("S")),3)
round(se(CDT.RT, facs=c("block")),3)
round(se(CDT.RT, facs=c("cond")),3)
round(se(CDT.RT, facs=c("S","cond")),3)
# Split up data object for cond comparisons
CDT.RT$S <- as.numeric(factor(CDT.RT$S))
bonf.CDT.RT <- data.frame()
bonf.CDT.RT <- tapply(CDT.RT$y,
                       list(s=CDT.RT$s,
                            cond=factor(CDT.RT$cond)), mean)
bonf.CDT.RT
t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,2], paired=T)
cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,2], method="paired")
t.test (bonf.CDT.RT[,2], bonf.CDT.RT[,3], paired=T)
cohensD(x=bonf.CDT.RT[,2], y=bonf.CDT.RT[,3], method="paired")
t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,3], paired=T)
cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,3], method="paired")

t.test (bonf.CDT.RT[,3], bonf.CDT.RT[,4], paired=T)
cohensD(x=bonf.CDT.RT[,3], y=bonf.CDT.RT[,4], method="paired")
t.test (bonf.CDT.RT[,2], bonf.CDT.RT[,4], paired=T)
cohensD(x=bonf.CDT.RT[,2], y=bonf.CDT.RT[,4], method="paired")
t.test (bonf.CDT.RT[,1], bonf.CDT.RT[,4], paired=T)
cohensD(x=bonf.CDT.RT[,1], y=bonf.CDT.RT[,4], method="paired")
# Split up data object for S*cond comparisons
CDT.RT$S <- as.numeric(factor(CDT.RT$S))
bonf.Conflict.CDT.RT <- data.frame()
bonf.Conflict.CDT.RT <- tapply(CDT.RT$y[CDT.RT$S==1],
                                list(s=CDT.RT$s[CDT.RT$S==1],
                                     cond=factor(CDT.RT$cond[CDT.RT$S==1])), mean)

t.test (bonf.Conflict.CDT.RT[,1], bonf.Conflict.CDT.RT[,2], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,1], y=bonf.Conflict.CDT.RT[,2], method="paired")
t.test (bonf.Conflict.CDT.RT[,2], bonf.Conflict.CDT.RT[,3], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,2], y=bonf.Conflict.CDT.RT[,3], method="paired")
t.test (bonf.Conflict.CDT.RT[,1], bonf.Conflict.CDT.RT[,3], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,1], y=bonf.Conflict.CDT.RT[,3], method="paired")

t.test (bonf.Conflict.CDT.RT[,3], bonf.Conflict.CDT.RT[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,3], y=bonf.Conflict.CDT.RT[,4], method="paired")
t.test (bonf.Conflict.CDT.RT[,2], bonf.Conflict.CDT.RT[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,2], y=bonf.Conflict.CDT.RT[,4], method="paired")
t.test (bonf.Conflict.CDT.RT[,1], bonf.Conflict.CDT.RT[,4], paired=T)
cohensD(x=bonf.Conflict.CDT.RT[,1], y=bonf.Conflict.CDT.RT[,4], method="paired")

CDT.RT$S <- as.numeric(factor(CDT.RT$S))
bonf.Nonconf.CDT.RT <- data.frame()
bonf.Nonconf.CDT.RT <- tapply(CDT.RT$y[CDT.RT$S==2],
                               list(s=CDT.RT$s[CDT.RT$S==2],
                                    cond=factor(CDT.RT$cond[CDT.RT$S==2])), mean)

t.test (bonf.Nonconf.CDT.RT[,1], bonf.Nonconf.CDT.RT[,2], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,1], y=bonf.Nonconf.CDT.RT[,2], method="paired")
t.test (bonf.Nonconf.CDT.RT[,2], bonf.Nonconf.CDT.RT[,3], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,2], y=bonf.Nonconf.CDT.RT[,3], method="paired")
t.test (bonf.Nonconf.CDT.RT[,1], bonf.Nonconf.CDT.RT[,3], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,1], y=bonf.Nonconf.CDT.RT[,3], method="paired")

t.test (bonf.Nonconf.CDT.RT[,3], bonf.Nonconf.CDT.RT[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,3], y=bonf.Nonconf.CDT.RT[,4], method="paired")
t.test (bonf.Nonconf.CDT.RT[,2], bonf.Nonconf.CDT.RT[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,2], y=bonf.Nonconf.CDT.RT[,4], method="paired")
t.test (bonf.Nonconf.CDT.RT[,1], bonf.Nonconf.CDT.RT[,4], paired=T)
cohensD(x=bonf.Nonconf.CDT.RT[,1], y=bonf.Nonconf.CDT.RT[,4], method="paired")



# # # ANOVA analysis: PM Accuracy # # #
#
#
# PM Accuracy Object
PMT.Acc <- arr2df(tapply(PMT$C,
                         list(s=PMT$s, cond=PMT$cond, S=PMT$S),mean))
wsAnova(PMT.Acc)
mneffects(PMT.Acc,list("S","cond",c("S","cond"), digits=3))
round(se(PMT.Acc, facs= c("S")),3)
round(se(PMT.Acc, facs= c("cond")),3)
# Split up data object for cond comparisons
PMT.Acc$S <- as.numeric(factor(PMT.Acc$S))
bonf.PMT.Acc <- data.frame()
bonf.PMT.Acc <- tapply(PMT.Acc$y,
                       list(s=PMT.Acc$s,
                            cond=factor(PMT.Acc$cond)), mean)
bonf.PMT.Acc
t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,2], paired=T)
cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,2], method="paired")
t.test (bonf.PMT.Acc[,2], bonf.PMT.Acc[,3], paired=T)
cohensD(x=bonf.PMT.Acc[,2], y=bonf.PMT.Acc[,3], method="paired")
t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,3], paired=T)
cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,3], method="paired")

t.test (bonf.PMT.Acc[,3], bonf.PMT.Acc[,4], paired=T)
cohensD(x=bonf.PMT.Acc[,3], y=bonf.PMT.Acc[,4], method="paired")
t.test (bonf.PMT.Acc[,2], bonf.PMT.Acc[,4], paired=T)
cohensD(x=bonf.PMT.Acc[,2], y=bonf.PMT.Acc[,4], method="paired")
t.test (bonf.PMT.Acc[,1], bonf.PMT.Acc[,4], paired=T)
cohensD(x=bonf.PMT.Acc[,1], y=bonf.PMT.Acc[,4], method="paired")

# # # ANOVA analysis: PM RT # # #
#
#
# PM RT Object
PMT.RT <- arr2df(tapply(PMT$RT[PMT$C==1],
                        list(s=PMT$s[PMT$C==1], cond=PMT$cond[PMT$C==1], S=PMT$S[PMT$C==1]), mean))
wsAnova(PMT.RT)
mneffects(PMT.RT,list("S","cond",c("S","cond"), digits=3))
round(se2(PMT.RT, facs=c("S")),3)
round(se2(PMT.RT, facs=c("cond")),3)
# Split up data object for cond comparisons
PMT.RT$S <- as.numeric(factor(PMT.RT$S))
bonf.PMT.RT <- data.frame()
bonf.PMT.RT <- tapply(PMT.RT$y,
                      list(s=PMT.RT$s,
                           cond=factor(PMT.RT$cond)), mean)
bonf.PMT.RT
t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,2], paired=T)
cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,2], method="paired")
t.test (bonf.PMT.RT[,2], bonf.PMT.RT[,3], paired=T)
cohensD(x=bonf.PMT.RT[,2], y=bonf.PMT.RT[,3], method="paired")
t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,3], paired=T)
cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,3], method="paired")

t.test (bonf.PMT.RT[,3], bonf.PMT.RT[,4], paired=T)
cohensD(x=bonf.PMT.RT[,3], y=bonf.PMT.RT[,4], method="paired")
t.test (bonf.PMT.RT[,2], bonf.PMT.RT[,4], paired=T)
cohensD(x=bonf.PMT.RT[,2], y=bonf.PMT.RT[,4], method="paired")
t.test (bonf.PMT.RT[,1], bonf.PMT.RT[,4], paired=T)
cohensD(x=bonf.PMT.RT[,1], y=bonf.PMT.RT[,4], method="paired")


# # # ANOVA analysis: Ongoing RT for PM versus non-PM Trials # # #
#
#
# Ongoing Task RT Object
head(data.E1)
CDT.Reactive <- data.E1[(data.E1$block=="3" & ((data.E1$S=="cc" & data.E1$R=="C") | (data.E1$S=="pc" & data.E1$R=="C") |
                                                   (data.E1$S=="nn" & data.E1$R=="N") | (data.E1$S=="pn" & data.E1$R=="N"))),]
CDT.Reactive$block <- factor(CDT.Reactive$block)
CDT.Reactive$R <- factor(CDT.Reactive$R)
any(CDT.Reactive$R=="P")
any(CDT.Reactive$block=="2")

str(CDT.Reactive)

CDT.Reactive <- arr2df(tapply(CDT$RT,
                        list(s=CDT$s,
                             cond=CDT$cond, S=CDT$S), mean))
str(CDT.Reactive)
CDT.Reactive.lmer.E1 <- lmer(y ~ S*cond+(1|s), data=CDT.Reactive)
save(CDT.Reactive.lmer.E1, file="CDT.Reactive.lmer.E1.RData")
load("CDT.Reactive.lmer.E1.RData")
CDT.Reactive.glm.E1 <- Anova(CDT.Reactive.lmer.E1,type="II")
CDT.Reactive.glm.E1

levels(CDT.Reactive$cond) <- c("LL.LT","LL.HT","HL.LT","HL.HT")
levels(CDT.Reactive$S) <- c("Conflict","Nonconflict","PM (Conflict)","PM (Nonconflict)")

wsAnova(CDT.Reactive)
mneffects(CDT.Reactive,list("S","cond",c("S","cond"), digits=3))
round(se2(CDT.Reactive, facs=c("S")),3)
round(se2(CDT.Reactive, facs=c("cond")),3)
# Split up data object for S comparisons
CDT.Reactive$S <- as.numeric(factor(CDT.Reactive$S))
bonf.CDT.Reactive <- data.frame()
bonf.CDT.Reactive <- tapply(CDT.Reactive$y,
                      list(s=CDT.Reactive$s,
                           S=factor(CDT.Reactive$S)), mean)
bonf.CDT.Reactive
t.test (bonf.CDT.Reactive[,1], bonf.CDT.Reactive[,2], paired=T)
cohensD(x=bonf.CDT.Reactive[,1], y=bonf.CDT.Reactive[,2], method="paired")
t.test (bonf.CDT.Reactive[,3], bonf.CDT.Reactive[,4], paired=T)
cohensD(x=bonf.CDT.Reactive[,3], y=bonf.CDT.Reactive[,4], method="paired")

t.test (bonf.CDT.Reactive[,1], bonf.CDT.Reactive[,3], paired=T)
cohensD(x=bonf.CDT.Reactive[,1], y=bonf.CDT.Reactive[,3], method="paired")
t.test (bonf.CDT.Reactive[,2], bonf.CDT.Reactive[,4], paired=T)
cohensD(x=bonf.CDT.Reactive[,2], y=bonf.CDT.Reactive[,4], method="paired")





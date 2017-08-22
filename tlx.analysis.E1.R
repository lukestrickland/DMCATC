###################### Russ data analysis template ###############################
# Running this top to bottom should correspond to my data analysis from the
# PMDC manuscript.
rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
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


# # # Prep dataframes for analysis # # #
#
tlx <- read.csv("data/E1_TLX.csv", header = TRUE)
tlx$s <- factor(tlx$s)
str(tlx)
head(tlx)

# tlx.PM <- tlx[ tlx$block=="PM",]
# tlx.PM$block <- factor(tlx.PM$block)
# str(tlx.PM)

pairs(tlx[ tlx$block=="PM",-c(1,2,3,4,6,8,10, 16,17,18,19)])

cor.test(tlx.PM$effort, tlx.PM$mean_v.P)

#
TLX.lmer.E1 <- lmer(effort ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1

tAB <- t.test (tlx[ tlx$cond=="A", "effort" ], tlx[ tlx$cond=="B", "effort" ], paired = TRUE)
dAB <- cohensD(tlx[ tlx$cond=="A", "effort" ], tlx[ tlx$cond=="B", "effort" ], method = "paired")
AB <- data.frame(cbind(t=tAB$statistic, df=tAB$parameter, p=tAB$p.value, Difference=tAB$estimate, d=dAB))

tBC <- t.test (tlx[ tlx$cond=="B", "effort" ], tlx[ tlx$cond=="C", "effort" ], paired = TRUE)
dBC <- cohensD(tlx[ tlx$cond=="B", "effort" ], tlx[ tlx$cond=="C", "effort" ], method = "paired")
BC <- data.frame(cbind(t=tBC$statistic, df=tBC$parameter, p=tBC$p.value, Difference=tBC$estimate, d=dBC))

tCD <- t.test (tlx[ tlx$cond=="C", "effort" ], tlx[ tlx$cond=="D", "effort" ], paired = TRUE)
dCD <- cohensD(tlx[ tlx$cond=="C", "effort" ], tlx[ tlx$cond=="D", "effort" ], method = "paired")
CD <- data.frame(cbind(t=tCD$statistic, df=tCD$parameter, p=tCD$p.value, Difference=tCD$estimate, d=dCD))

T.Contrasts.TLX.Effort.TP <- data.frame(rbind(AB, BC, CD))
T.Contrasts.TLX.Effort.TP$t <- round(T.Contrasts.TLX.Effort.TP$t, digits=2)
T.Contrasts.TLX.Effort.TP$Difference <- round(T.Contrasts.TLX.Effort.TP$Difference, digits=2)
T.Contrasts.TLX.Effort.TP$d <- round(T.Contrasts.TLX.Effort.TP$d, digits=2)
T.Contrasts.TLX.Effort.TP$p <- round(T.Contrasts.TLX.Effort.TP$p, digits=3)
T.Contrasts.TLX.Effort.TP$p <- format.pval(T.Contrasts.TLX.Effort.TP$p, digits=2, eps= 0.001)
row.names(T.Contrasts.TLX.Effort.TP) <- c("A-B","B-C","C-D")
# T.Contrasts.TLX.Effort.TP
pander(T.Contrasts.TLX.Effort.TP)
write.csv(T.Contrasts.TLX.Effort.TP, file="analysis/T.Contrasts.TLX.Effort.TP.E1.csv")


ggplot(tlx, aes(x=cond, y=effort, fill=block)) +
    geom_boxplot()


#
TLX.lmer.E1 <- lmer(effort ~ mean_v.C+mean_v.N+mean_v.Cerr+mean_v.Nerr+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1










TLX.lmer.E1 <- lmer(effort ~ mean_v.Cerr+mean_v.Nerr+mean_v.P+(1|s), data=tlx.PM)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1

TLX.lmer.E1 <- lmer(mean_v.P ~ mental+temporal+effort+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1

ggplot(tlx, aes(x=effort, y=mean_v.P)) +
    geom_point(shape=1,color="darkblue") +
    geom_smooth(method=lm,color="darkblue") +
    geom_point(aes(x=effort, y=mean_v.C),shape=1,color="forestgreen") +
    geom_smooth(aes(x=effort, y=mean_v.C),method=lm,color="forestgreen") +
    geom_point(aes(x=effort, y=mean_v.N),shape=1,color="forestgreen") +
    geom_smooth(aes(x=effort, y=mean_v.N),method=lm,color="forestgreen") +
    geom_point(aes(x=effort, y=mean_v.Cerr),shape=1,color="red") +
    geom_smooth(aes(x=effort, y=mean_v.Cerr),method=lm,color="red") +
    geom_point(aes(x=effort, y=mean_v.Nerr),shape=1,color="red") +
    geom_smooth(aes(x=effort, y=mean_v.Nerr),method=lm,color="red") +
    xlab("TLX Effort Rating") +
    ylab("Drift Rate") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 15, vjust=0.3))




ggplot(tlx, aes(x=cond, y=mental, fill=block)) +
    geom_boxplot()

ggplot(tlx, aes(x=cond, y=temporal, fill=block)) +
    geom_boxplot()

ggplot(tlx, aes(x=cond, y=physical, fill=block)) +
    geom_boxplot()

ggplot(tlx, aes(x=cond, y=performance, fill=block)) +
    geom_boxplot()

ggplot(tlx, aes(x=cond, y=frustration, fill=block)) +
    geom_boxplot()


TLX.lmer.E1 <- lmer(mental ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1
TLX.lmer.E1 <- lmer(physical ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1
TLX.lmer.E1 <- lmer(temporal ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1
TLX.lmer.E1 <- lmer(performance ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1
TLX.lmer.E1 <- lmer(effort ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1
TLX.lmer.E1 <- lmer(frustration ~ block*cond+(1|s), data=tlx)
TLX.lm.E1 <- Anova(TLX.lmer.E1,type="II")
TLX.lm.E1






# # # Calculate and save off averaged parameter objects for each subject for analysis # # #
#
# # #
# subj.meanthetas <- function (samples){
#     samps <- lapply(samples, function(x) x["theta"])
#     ## thetas into big array for apply
#     samps2<- unlist(samps)
#     dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
#     dim(samps2) <- dim3
#     samps3<- apply(samps2, c(4,2), mean)
#     ## back to a theta list after applied
#     colnames(samps3)<- colnames(samps[[1]]$theta)
#     df <- cbind(names(samples), data.frame(samps3))
#     names(df)[1] <- "s"
#     df
# }
# # # #
#
# subj.thetas.E1 <- subj.meanthetas(samples)
# save(subj.thetas.E1, file="data/after_sampling/subj.thetas.E1.RData")
#
# rm(samples)
#
# #






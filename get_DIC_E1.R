rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")
load_model ("LBA","lbaN_B.R")
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
# install.packages(pkgs) #install
sapply(pkgs, require, character.only = T) #load

# # # DIC Model Selection # # #
#
#     E1 - Neutral Emphasis
#

load("data/samples/E1.block.B.V_cond.B.V.PMFA.samples.RData")  # Top Model
E1.PMFA <- h.IC.dmc(E1.block.B.V_cond.B.V.PMFA.samples,DIC=TRUE)
E1.PMFA.npars <- E1.block.B.V_cond.B.V.PMFA.samples$p1$n.pars
rm(E1.block.B.V_cond.B.V.PMFA.samples)

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Selected Model
E1.PMV <- h.IC.dmc(E1.block.B.V_cond.B.V.PMV.samples,DIC=TRUE)
E1.PMV.npars <- E1.block.B.V_cond.B.V.PMV.samples$p1$n.pars
rm(E1.block.B.V_cond.B.V.PMV.samples)

load("data/samples/E1.blockVonly.samples.RData")
E1.blockVonly <- h.IC.dmc(E1.blockVonly.samples,DIC=TRUE)
E1.blockVonly.npars <- E1.blockVonly.samples$p1$n.pars
rm(E1.blockVonly.samples)

load("data/samples/E1.blockBonly.samples.RData")
E1.blockBonly <- h.IC.dmc(E1.blockBonly.samples,DIC=TRUE)
E1.blockBonly.npars <- E1.blockBonly.samples$p1$n.pars
rm(E1.blockBonly.samples)

load("data/samples/E1.condVonly.samples.RData")
E1.condVonly <- h.IC.dmc(E1.condVonly.samples,DIC=TRUE)
E1.condVonly.npars <- E1.condVonly.samples$p1$n.pars
rm(E1.condVonly.samples)

load("data/samples/E1.condBonly.samples.RData")
E1.condBonly <- h.IC.dmc(E1.condBonly.samples,DIC=TRUE)
E1.condBonly.npars <- E1.condBonly.samples$p1$n.pars
rm(E1.condBonly.samples)



sum(E1.PMFA)
sum(E1.PMV)
sum(E1.blockVonly)
sum(E1.blockBonly)
sum(E1.condVonly)
sum(E1.condBonly)

DIC.TABLE.E1 <- data.frame(cbind(Model=c("Top Model","Selected Model",
                                         "Selected Model with B fixed over PM Block",
                                         "Selected Model with V fixed over PM Block",
                                         "Selected Model with B fixed over Time Pressure",
                                         "Selected Model with V fixed over Time Pressure"
                                         ),
                                 n.pars=c(E1.PMFA.npars,
                                          E1.PMV.npars,
                                          E1.blockVonly.npars,
                                          E1.blockBonly.npars,
                                          E1.condVonly.npars,
                                          E1.condBonly.npars
                                          ),
                                 DIC=round(c(sum(E1.PMFA),
                                             sum(E1.PMV),
                                             sum(E1.blockVonly),
                                             sum(E1.blockBonly),
                                             sum(E1.condVonly),
                                             sum(E1.condBonly)
                                             ),digits = 3)))
colnames(DIC.TABLE.E1) <- c("Model","Number of Parameters","DIC")
DIC.TABLE.E1
write.csv(DIC.TABLE.E1, file="analysis/DIC.TABLE.E1.csv")


# # # Comparison of summed DICs
#
#     DIC difference > 10 rules out model with larger DIC
#

# # # Group DIC function # # #

# h.IC.dmc <- function(hsamples,DIC=FALSE,fast=TRUE,use.pd=NA)
#   # Applies IC.dmc each suubject, prints sum.
#   # Invisibly returns results for each subject.
# {
#   if ( any(names(hsamples)=="theta") )
#     stop("For a single subject use Dstats.dmc")
#   ds <- lapply(hsamples,Dstats.dmc,fast=fast,save=FALSE)
#   pds <- lapply(ds,pd.dmc)
#   if ( is.na(use.pd) ) {
#     pd <- numeric(length(pds))
#     for (i in 1:length(pds)) {
#       if (ds[[i]]$minD < ds[[i]]$Dmean)
#         pd[i] <- pds[[i]]$Pmin else
#           pd[i] <- pds[[i]]$Pmean
#     }
#   } else {
#     if ( use.pd %in% names(pds[[1]]) )
#       pd <- unlist(lapply(pds,function(x){x[[use.pd]]})) else
#         stop(paste("use.pd must be one of:",paste(names(pds[[1]]),collapse=",")))
#   }
#   ICs <- numeric(length(ds))
#   for (i in 1:length(ds)) if (DIC)
#     ICs[i] <- ds[[i]]$meanD+pd[i] else
#       ICs[i] <- ds[[i]]$meanD+2*pd[i]
#   if (DIC) cat("Summed DIC\n") else cat("Summed BPIC\n")
#   print(sum(ICs))
#   invisible(ICs)
# }


rm(list=ls())
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")
load_model ("LBA","lbaN_B.R")
# source("LSAnova.R")
require("gridExtra")
require("lme4")
require("car")
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
sapply(pkgs, require, character.only = T) #load

# # # Load samples object # # #
#
load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.E1 <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)

samples.E1 <- samples.E1[names(samples.E1) != "p17"]  # Exclude p17 E1 due to no PM responses


# # # Load functions # # #
#
group.inference.dist <- function (hsamples, fun) {
    inference <- list()
    for (i in 1:length(hsamples)) {
        thetas <- hsamples[[i]]$theta
        inference [[i]] <- fun (thetas)
    }
    inf2 <- unlist(inference)
    dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
    dim(inf2) <- dim3
    apply(inf2, c(1,2,3), mean)
}

minp <- function (effect) min(ecdf(effect)(0), 1-ecdf(effect)(0))

zandp <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), "(", round(p,3), ")", sep="")
}

mean.sd <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    M <- mean(effect)
    SD <- sd(effect)
    data.frame(M, SD)
}

# # # Reactive Control # # #
#
Reactive.ccAC <- function (thetas) thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.pcA3C",, drop=F]
Reactive.ccBC <- function (thetas) thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.pcB3C",, drop=F]
Reactive.ccCC <- function (thetas) thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.pcC3C",, drop=F]
Reactive.ccDC <- function (thetas) thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.pcD3C",, drop=F]

Reactive.nnAN <- function (thetas) thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.pnA3N",, drop=F]
Reactive.nnBN <- function (thetas) thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.pnB3N",, drop=F]
Reactive.nnCN <- function (thetas) thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.pnC3N",, drop=F]
Reactive.nnDN <- function (thetas) thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.pnD3N",, drop=F]

Reactive.nnAC <- function (thetas) thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.pnA3C",, drop=F]
Reactive.nnBC <- function (thetas) thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F]
Reactive.nnCC <- function (thetas) thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F]
Reactive.nnDC <- function (thetas) thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F]

Reactive.ccAN <- function (thetas) thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.pcA3N",, drop=F]
Reactive.ccBN <- function (thetas) thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F]
Reactive.ccCN <- function (thetas) thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F]
Reactive.ccDN <- function (thetas) thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F]

Reactive.Table <- data.frame(rbind(
    A.NonPM.minus.PM=c(zandp(samples.E1, Reactive.ccAC),
                       zandp(samples.E1, Reactive.nnAN),
                       zandp(samples.E1, Reactive.nnAC),
                       zandp(samples.E1, Reactive.ccAN)
    ),

    B.NonPM.minus.PM=c(zandp(samples.E1, Reactive.ccBC),
                       zandp(samples.E1, Reactive.nnBN),
                       zandp(samples.E1, Reactive.nnBC),
                       zandp(samples.E1, Reactive.ccBN)
    ),

    C.NonPM.minus.PM=c(zandp(samples.E1, Reactive.ccCC),
                       zandp(samples.E1, Reactive.nnCN),
                       zandp(samples.E1, Reactive.nnCC),
                       zandp(samples.E1, Reactive.ccCN)
    ),

    D.NonPM.minus.PM=c(zandp(samples.E1, Reactive.ccDC),
                       zandp(samples.E1, Reactive.nnDN),
                       zandp(samples.E1, Reactive.nnDC),
                       zandp(samples.E1, Reactive.ccDN)
    )
))
colnames(Reactive.Table) <- c("Conflict","Nonconflict","Conflict (Error)","Nonconflict (Error)")
Reactive.Table


# # # Reactive Control - Inhibition of Correct minus Error Stimuli # # #
#
Reactive.AC <- function (thetas) (thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.pcA3C",, drop=F]) -
    (thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.pnA3C",, drop=F])
Reactive.BC <- function (thetas) (thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.pcB3C",, drop=F]) -
    (thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F])
Reactive.CC <- function (thetas) (thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.pcC3C",, drop=F]) -
    (thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F])
Reactive.DC <- function (thetas) (thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.pcD3C",, drop=F]) -
    (thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F])

Reactive.AN <- function (thetas) (thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.pnA3N",, drop=F]) -
    (thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.pcA3N",, drop=F])
Reactive.BN <- function (thetas) (thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.pnB3N",, drop=F]) -
    (thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F])
Reactive.CN <- function (thetas) (thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.pnC3N",, drop=F]) -
    (thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F])
Reactive.DN <- function (thetas) (thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.pnD3N",, drop=F]) -
    (thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F])

Reactive.Cor.minus.Err.Table <- data.frame(rbind(
    A.Cor.minus.Err=c(zandp(samples.E1, Reactive.AC),
                      zandp(samples.E1, Reactive.AN)
    ),

    B.Cor.minus.Err=c(zandp(samples.E1, Reactive.BC),
                      zandp(samples.E1, Reactive.BN)
    ),

    C.Cor.minus.Err=c(zandp(samples.E1, Reactive.CC),
                      zandp(samples.E1, Reactive.CN)
    ),

    D.Cor.minus.Err=c(zandp(samples.E1, Reactive.DC),
                      zandp(samples.E1, Reactive.DN)
    )
))
colnames(Reactive.Cor.minus.Err.Table) <- c("Conflict Correct-Error Diff","Nonconflict Correct-Error Diff")
Reactive.Cor.minus.Err.Table

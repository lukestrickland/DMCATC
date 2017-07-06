
rm(list=ls())
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")
# source("LSAnova.R")
require(gridExtra)
require("lme4")
require(car)
require(plyr)
require(dplyr)
require("pander")

load("C:/Users/Russell Boag/Documents/GitHub/DMCATC/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object

samples.E1 <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)


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


## Reactive Control

Reactive.ccAC <- function (thetas) thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.pcA3C",, drop=F]
Reactive.ccBC <- function (thetas) thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.pcB3C",, drop=F]
Reactive.ccCC <- function (thetas) thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.pcC3C",, drop=F]
Reactive.ccDC <- function (thetas) thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.pcD3C",, drop=F]

Reactive.nnAN <- function (thetas) thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.pnA3N",, drop=F]
Reactive.nnBN <- function (thetas) thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.pnB3N",, drop=F]
Reactive.nnCN <- function (thetas) thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.pnC3N",, drop=F]
Reactive.nnDN <- function (thetas) thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.pnD3N",, drop=F]

Reactive.ccAN <- function (thetas) thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.pcA3N",, drop=F]
Reactive.ccBN <- function (thetas) thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F]
Reactive.ccCN <- function (thetas) thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F]
Reactive.ccDN <- function (thetas) thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F]

Reactive.nnAC <- function (thetas) thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.pnA3C",, drop=F]
Reactive.nnBC <- function (thetas) thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F]
Reactive.nnCC <- function (thetas) thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F]
Reactive.nnDC <- function (thetas) thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F]



HFreacpN <- function (thetas) (thetas[,"mean_v.FwN",, drop=F] - thetas[,"mean_v.FpN",, drop=F]) - (thetas[,"mean_v.HwN",, drop=F] - thetas[,"mean_v.HpN",, drop=F])


# HFreacpW <- function (thetas) (thetas[,"mean_v.FwW",, drop=F] - thetas[,"mean_v.FpW",, drop=F]) - (thetas[,"mean_v.HwW",, drop=F] - thetas[,"mean_v.HpW",, drop=F])
#
# HFreacpP <- function (thetas) (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.HpP",, drop=F]) - (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.FpP",, drop=F])
# HFreacpP <- function (thetas) (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.HpP",, drop=F]) - (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.FpP",, drop=F])


zandp <- function(samples, fun){
	effect<- group.inference.dist(samples, fun)
	Z <- mean(effect)/sd(effect)
	p <- minp(effect)
	paste(round(Z,2), "(", round(p,3), ")", sep="")
}


Reactive.Table <- data.frame(rbind(
    A=c(zandp(samples.E1, Reactive.ccAC),
        zandp(samples.E1, Reactive.nnAN),
        zandp(samples.E1, Reactive.nnAC),
        zandp(samples.E1, Reactive.ccAN)
        ),

    B=c(zandp(samples.E1, Reactive.ccBC),
        zandp(samples.E1, Reactive.nnBN),
        zandp(samples.E1, Reactive.nnBC),
        zandp(samples.E1, Reactive.ccBN)
        ),

    C=c(zandp(samples.E1, Reactive.ccCC),
        zandp(samples.E1, Reactive.nnCN),
        zandp(samples.E1, Reactive.nnCC),
        zandp(samples.E1, Reactive.ccCN)
        ),

    D=c(zandp(samples.E1, Reactive.ccDC),
        zandp(samples.E1, Reactive.nnDN),
        zandp(samples.E1, Reactive.nnDC),
        zandp(samples.E1, Reactive.ccDN)
        )
))

colnames(Reactive.Table) <- c("Conflict","Nonconflict","Conflict (FA)","Nonconflict (FA)")

Reactive.Table



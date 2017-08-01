
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

# # # Load samples object # # #
#
load("C:/Users/Russell Boag/Documents/GitHub/DMCATC/samples/A4.block.B.V_cond.B.V.PMV.samples.RData")  # Samples object
samples.A4 <- A4.block.B.V_cond.B.V.PMV.samples
rm(A4.block.B.V_cond.B.V.PMV.samples)

# samples.A4[[1]]$theta

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

Reactive.ccAN <- function (thetas) thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.pcA3N",, drop=F]
Reactive.ccBN <- function (thetas) thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F]
Reactive.ccCN <- function (thetas) thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F]
Reactive.ccDN <- function (thetas) thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F]

Reactive.nnAC <- function (thetas) thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.pnA3C",, drop=F]
Reactive.nnBC <- function (thetas) thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F]
Reactive.nnCC <- function (thetas) thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F]
Reactive.nnDC <- function (thetas) thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F]



# HFreacpN <- function (thetas) (thetas[,"mean_v.FwN",, drop=F] - thetas[,"mean_v.FpN",, drop=F]) - (thetas[,"mean_v.HwN",, drop=F] - thetas[,"mean_v.HpN",, drop=F])


# HFreacpW <- function (thetas) (thetas[,"mean_v.FwW",, drop=F] - thetas[,"mean_v.FpW",, drop=F]) - (thetas[,"mean_v.HwW",, drop=F] - thetas[,"mean_v.HpW",, drop=F])
#
# HFreacpP <- function (thetas) (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.HpP",, drop=F]) - (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.FpP",, drop=F])
# HFreacpP <- function (thetas) (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.HpP",, drop=F]) - (thetas[,"mean_v.fa",, drop=F] - thetas[,"mean_v.FpP",, drop=F])


# # # Proactive Control # # #
#
Proactive.AC <- function (thetas) thetas[,"B.A3C",, drop=F] - thetas[,"B.A2C",, drop=F]
Proactive.BC <- function (thetas) thetas[,"B.B3C",, drop=F] - thetas[,"B.B2C",, drop=F]
Proactive.CC <- function (thetas) thetas[,"B.C3C",, drop=F] - thetas[,"B.C2C",, drop=F]
Proactive.DC <- function (thetas) thetas[,"B.D3C",, drop=F] - thetas[,"B.D2C",, drop=F]

Proactive.AN <- function (thetas) thetas[,"B.A3N",, drop=F] - thetas[,"B.A2N",, drop=F]
Proactive.BN <- function (thetas) thetas[,"B.B3N",, drop=F] - thetas[,"B.B2N",, drop=F]
Proactive.CN <- function (thetas) thetas[,"B.C3N",, drop=F] - thetas[,"B.C2N",, drop=F]
Proactive.DN <- function (thetas) thetas[,"B.D3N",, drop=F] - thetas[,"B.D2N",, drop=F]


# # # Capacity (Difference in Drift between Control and PM blocks) # # #
#
Capacity.ccAC <- function (thetas) thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.ccA2C",, drop=F]
Capacity.ccBC <- function (thetas) thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.ccB2C",, drop=F]
Capacity.ccCC <- function (thetas) thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.ccC2C",, drop=F]
Capacity.ccDC <- function (thetas) thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.ccD2C",, drop=F]

Capacity.nnAN <- function (thetas) thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.nnA2N",, drop=F]
Capacity.nnBN <- function (thetas) thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.nnB2N",, drop=F]
Capacity.nnCN <- function (thetas) thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.nnC2N",, drop=F]
Capacity.nnDN <- function (thetas) thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.nnD2N",, drop=F]

Capacity.ccAN <- function (thetas) thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.ccA2N",, drop=F]
Capacity.ccBN <- function (thetas) thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.ccB2N",, drop=F]
Capacity.ccCN <- function (thetas) thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.ccC2N",, drop=F]
Capacity.ccDN <- function (thetas) thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.ccD2N",, drop=F]

Capacity.nnAC <- function (thetas) thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.nnA2C",, drop=F]
Capacity.nnBC <- function (thetas) thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.nnB2C",, drop=F]
Capacity.nnCC <- function (thetas) thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.nnC2C",, drop=F]
Capacity.nnDC <- function (thetas) thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.nnD2C",, drop=F]


# # # Proactive Control (PM Cost) by Time Pressure # # #
#
ProactiveTP.ABC <- function (thetas) (thetas[,"B.A3C",, drop=F] - thetas[,"B.A2C",, drop=F]) -
    (thetas[,"B.B3C",, drop=F] - thetas[,"B.B2C",, drop=F])
ProactiveTP.BCC <- function (thetas) (thetas[,"B.B3C",, drop=F] - thetas[,"B.B2C",, drop=F]) -
    (thetas[,"B.C3C",, drop=F] - thetas[,"B.C2C",, drop=F])
ProactiveTP.CDC <- function (thetas) (thetas[,"B.C3C",, drop=F] - thetas[,"B.C2C",, drop=F]) -
    (thetas[,"B.D3C",, drop=F] - thetas[,"B.D2C",, drop=F])

ProactiveTP.ABN <- function (thetas) (thetas[,"B.A3N",, drop=F] - thetas[,"B.A2N",, drop=F]) -
    (thetas[,"B.B3N",, drop=F] - thetas[,"B.B2N",, drop=F])
ProactiveTP.BCN <- function (thetas) (thetas[,"B.B3N",, drop=F] - thetas[,"B.B2N",, drop=F]) -
    (thetas[,"B.C3N",, drop=F] - thetas[,"B.C2N",, drop=F])
ProactiveTP.CDN <- function (thetas) (thetas[,"B.C3N",, drop=F] - thetas[,"B.C2N",, drop=F]) -
    (thetas[,"B.D3N",, drop=F] - thetas[,"B.D2N",, drop=F])


# # # Effort/Arousal by Time Pressure # # #
#
Effort.TPccCABC <- function (thetas) thetas[,"mean_v.ccA2C",, drop=F] - thetas[,"mean_v.ccB2C",, drop=F]
Effort.TPccCBCC <- function (thetas) thetas[,"mean_v.ccB2C",, drop=F] - thetas[,"mean_v.ccC2C",, drop=F]
Effort.TPccCCDC <- function (thetas) thetas[,"mean_v.ccC2C",, drop=F] - thetas[,"mean_v.ccD2C",, drop=F]

Effort.TPccPMABC <- function (thetas) thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.ccB3C",, drop=F]
Effort.TPccPMBCC <- function (thetas) thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.ccC3C",, drop=F]
Effort.TPccPMCDC <- function (thetas) thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.ccD3C",, drop=F]

Effort.TPnnCABN <- function (thetas) thetas[,"mean_v.nnA2N",, drop=F] - thetas[,"mean_v.nnB2N",, drop=F]
Effort.TPnnCBCN <- function (thetas) thetas[,"mean_v.nnB2N",, drop=F] - thetas[,"mean_v.nnC2N",, drop=F]
Effort.TPnnCCDN <- function (thetas) thetas[,"mean_v.nnC2N",, drop=F] - thetas[,"mean_v.nnD2N",, drop=F]

Effort.TPnnPMABN <- function (thetas) thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.nnB3N",, drop=F]
Effort.TPnnPMBCN <- function (thetas) thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.nnC3N",, drop=F]
Effort.TPnnPMCDN <- function (thetas) thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.nnD3N",, drop=F]


#
Effort.TPnnCABC <- function (thetas) thetas[,"mean_v.nnA2C",, drop=F] - thetas[,"mean_v.nnB2C",, drop=F]
Effort.TPnnCBCC <- function (thetas) thetas[,"mean_v.nnB2C",, drop=F] - thetas[,"mean_v.nnC2C",, drop=F]
Effort.TPnnCCDC <- function (thetas) thetas[,"mean_v.nnC2C",, drop=F] - thetas[,"mean_v.nnD2C",, drop=F]

Effort.TPnnPMABC <- function (thetas) thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.nnB3C",, drop=F]
Effort.TPnnPMBCC <- function (thetas) thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.nnC3C",, drop=F]
Effort.TPnnPMCDC <- function (thetas) thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.nnD3C",, drop=F]

Effort.TPccCABN <- function (thetas) thetas[,"mean_v.ccA2N",, drop=F] - thetas[,"mean_v.ccB2N",, drop=F]
Effort.TPccCBCN <- function (thetas) thetas[,"mean_v.ccB2N",, drop=F] - thetas[,"mean_v.ccC2N",, drop=F]
Effort.TPccCCDN <- function (thetas) thetas[,"mean_v.ccC2N",, drop=F] - thetas[,"mean_v.ccD2N",, drop=F]

Effort.TPccPMABN <- function (thetas) thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.ccB3N",, drop=F]
Effort.TPccPMBCN <- function (thetas) thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.ccC3N",, drop=F]
Effort.TPccPMCDN <- function (thetas) thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.ccD3N",, drop=F]

#
Effort.TPppPMABP <- function (thetas) thetas[,"mean_v.ppA3P",, drop=F] - thetas[,"mean_v.ppB3P",, drop=F]
Effort.TPppPMBCP <- function (thetas) thetas[,"mean_v.ppB3P",, drop=F] - thetas[,"mean_v.ppC3P",, drop=F]
Effort.TPppPMCDP <- function (thetas) thetas[,"mean_v.ppC3P",, drop=F] - thetas[,"mean_v.ppD3P",, drop=F]

Effort.TPpcPMABC <- function (thetas) thetas[,"mean_v.pcA3C",, drop=F] - thetas[,"mean_v.pcB3C",, drop=F]
Effort.TPpcPMBCC <- function (thetas) thetas[,"mean_v.pcB3C",, drop=F] - thetas[,"mean_v.pcC3C",, drop=F]
Effort.TPpcPMCDC <- function (thetas) thetas[,"mean_v.pcC3C",, drop=F] - thetas[,"mean_v.pcD3C",, drop=F]

Effort.TPpnPMABN <- function (thetas) thetas[,"mean_v.pnA3N",, drop=F] - thetas[,"mean_v.pnB3N",, drop=F]
Effort.TPpnPMBCN <- function (thetas) thetas[,"mean_v.pnB3N",, drop=F] - thetas[,"mean_v.pnC3N",, drop=F]
Effort.TPpnPMCDN <- function (thetas) thetas[,"mean_v.pnC3N",, drop=F] - thetas[,"mean_v.pnD3N",, drop=F]

Effort.TPpnPMABC <- function (thetas) thetas[,"mean_v.pnA3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F]
Effort.TPpnPMBCC <- function (thetas) thetas[,"mean_v.pnB3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F]
Effort.TPpnPMCDC <- function (thetas) thetas[,"mean_v.pnC3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F]

Effort.TPpcPMABN <- function (thetas) thetas[,"mean_v.pcA3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F]
Effort.TPpcPMBCN <- function (thetas) thetas[,"mean_v.pcB3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F]
Effort.TPpcPMCDN <- function (thetas) thetas[,"mean_v.pcC3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F]


# # # Z-score and P-value Tables # # #
#
Reactive.Table <- data.frame(rbind(
    A.NonPM.minus.PM=c(zandp(samples.A4, Reactive.ccAC),
        zandp(samples.A4, Reactive.nnAN),
        zandp(samples.A4, Reactive.nnAC),
        zandp(samples.A4, Reactive.ccAN)
        ),

    B.NonPM.minus.PM=c(zandp(samples.A4, Reactive.ccBC),
        zandp(samples.A4, Reactive.nnBN),
        zandp(samples.A4, Reactive.nnBC),
        zandp(samples.A4, Reactive.ccBN)
        ),

    C.NonPM.minus.PM=c(zandp(samples.A4, Reactive.ccCC),
        zandp(samples.A4, Reactive.nnCN),
        zandp(samples.A4, Reactive.nnCC),
        zandp(samples.A4, Reactive.ccCN)
        ),

    D.NonPM.minus.PM=c(zandp(samples.A4, Reactive.ccDC),
        zandp(samples.A4, Reactive.nnDN),
        zandp(samples.A4, Reactive.nnDC),
        zandp(samples.A4, Reactive.ccDN)
        )
))
colnames(Reactive.Table) <- c("Conflict","Nonconflict","Conflict (FA)","Nonconflict (FA)")
# Reactive.Table

#
Proactive.Table <- data.frame(rbind(
    A.PM.minus.Control=c(zandp(samples.A4, Proactive.AC),
        zandp(samples.A4, Proactive.AN)
    ),

    B.PM.minus.Control=c(zandp(samples.A4, Proactive.BC),
        zandp(samples.A4, Proactive.BN)
    ),

    C.PM.minus.Control=c(zandp(samples.A4, Proactive.CC),
        zandp(samples.A4, Proactive.CN)
    ),

    D.PM.minus.Control=c(zandp(samples.A4, Proactive.DC),
        zandp(samples.A4, Proactive.DN)
    )
))
colnames(Proactive.Table) <- c("Conflict","Nonconflict")
# Proactive.Table

#
Capacity.Table <- data.frame(rbind(
    A.PM.minus.Control=c(zandp(samples.A4, Capacity.ccAC),
                         zandp(samples.A4, Capacity.nnAN),
                         zandp(samples.A4, Capacity.nnAC),
                         zandp(samples.A4, Capacity.ccAN)
    ),

    B.PM.minus.Control=c(zandp(samples.A4, Capacity.ccBC),
                         zandp(samples.A4, Capacity.nnBN),
                         zandp(samples.A4, Capacity.nnBC),
                         zandp(samples.A4, Capacity.ccBN)
    ),

    C.PM.minus.Control=c(zandp(samples.A4, Capacity.ccCC),
                         zandp(samples.A4, Capacity.nnCN),
                         zandp(samples.A4, Capacity.nnCC),
                         zandp(samples.A4, Capacity.ccCN)
    ),

    D.PM.minus.Control=c(zandp(samples.A4, Capacity.ccDC),
                         zandp(samples.A4, Capacity.nnDN),
                         zandp(samples.A4, Capacity.nnDC),
                         zandp(samples.A4, Capacity.ccDN)
    )
))
colnames(Capacity.Table) <- c("Conflict","Nonconflict","Conflict (FA)","Nonconflict (FA)")
# Capacity.Table

#
ProactiveTP.Table <- data.frame(rbind(
    PM.Cost.A.minus.B=c(zandp(samples.A4, ProactiveTP.ABC),zandp(samples.A4, ProactiveTP.ABN)
    ),

    PM.Cost.B.minus.C=c(zandp(samples.A4, ProactiveTP.BCC),zandp(samples.A4, ProactiveTP.BCN)
    ),

    PM.Cost.C.minus.D=c(zandp(samples.A4, ProactiveTP.CDC),zandp(samples.A4, ProactiveTP.CDN)
    )
))
colnames(ProactiveTP.Table) <- c("Conflict","Nonconflict")
# ProactiveTP.Table

#
Effort.by.TP.Table.Wide <- data.frame(rbind(
    Conflict=c(zandp(samples.A4, Effort.TPccCABC),
               zandp(samples.A4, Effort.TPccCBCC),
               zandp(samples.A4, Effort.TPccCCDC),

               zandp(samples.A4, Effort.TPccPMABC),
               zandp(samples.A4, Effort.TPccPMBCC),
               zandp(samples.A4, Effort.TPccPMCDC)
    ),

    Nonconflict=c(zandp(samples.A4, Effort.TPnnCABN),
                  zandp(samples.A4, Effort.TPnnCBCN),
                  zandp(samples.A4, Effort.TPnnCCDN),

                  zandp(samples.A4, Effort.TPnnPMABN),
                  zandp(samples.A4, Effort.TPnnPMBCN),
                  zandp(samples.A4, Effort.TPnnPMCDN)
    ),

    Conflict.FA=c(zandp(samples.A4, Effort.TPnnCABC),
                  zandp(samples.A4, Effort.TPnnCBCC),
                  zandp(samples.A4, Effort.TPnnCCDC),

                  zandp(samples.A4, Effort.TPnnPMABC),
                  zandp(samples.A4, Effort.TPnnPMBCC),
                  zandp(samples.A4, Effort.TPnnPMCDC)
    ),

    Nonconflict.FA=c(zandp(samples.A4, Effort.TPccCABN),
                     zandp(samples.A4, Effort.TPccCBCN),
                     zandp(samples.A4, Effort.TPccCCDN),

                     zandp(samples.A4, Effort.TPccPMABN),
                     zandp(samples.A4, Effort.TPccPMBCN),
                     zandp(samples.A4, Effort.TPccPMCDN)
    ),
    PM.to.PM=c("-","-","-",
         zandp(samples.A4, Effort.TPppPMABP),
         zandp(samples.A4, Effort.TPppPMBCP),
         zandp(samples.A4, Effort.TPppPMCDP)
    ),
    Conflict.to.PMC=c("-","-","-",
                      zandp(samples.A4, Effort.TPpcPMABC),
                      zandp(samples.A4, Effort.TPpcPMBCC),
                      zandp(samples.A4, Effort.TPpcPMCDC)
    ),
    Nonconflict.to.PMN=c("-","-","-",
                      zandp(samples.A4, Effort.TPpnPMABN),
                      zandp(samples.A4, Effort.TPpnPMBCN),
                      zandp(samples.A4, Effort.TPpnPMCDN)
    ),
    Conflict.to.PMN=c("-","-","-",
                      zandp(samples.A4, Effort.TPpnPMABC),
                      zandp(samples.A4, Effort.TPpnPMBCC),
                      zandp(samples.A4, Effort.TPpnPMCDC)
    ),
    Nonconflict.to.PMC=c("-","-","-",
                        zandp(samples.A4, Effort.TPpcPMABN),
                        zandp(samples.A4, Effort.TPpcPMBCN),
                        zandp(samples.A4, Effort.TPpcPMCDN)
    )
))
colnames(Effort.by.TP.Table.Wide) <- c("A-B (Control)","B-C (Control)","C-D (Control)","A-B (PM)","B-C (PM)","C-D (PM)")
# Effort.by.TP.Table.Wide

#
# Effort.by.TP.Table.Long <- data.frame(rbind(
#     Control.A.minus.B=c(zandp(samples.A4, Effort.TPccCABC),
#                         zandp(samples.A4, Effort.TPnnCABN),
#                         zandp(samples.A4, Effort.TPnnCABC),
#                         zandp(samples.A4, Effort.TPccCABN),
#                         "-","-","-","-","-"
#                         ),
#     Control.B.minus.C=c(zandp(samples.A4, Effort.TPccCBCC),
#                         zandp(samples.A4, Effort.TPnnCBCN),
#                         zandp(samples.A4, Effort.TPnnCBCC),
#                         zandp(samples.A4, Effort.TPccCBCN),
#                         "-","-","-","-","-"
#     ),
#     Control.C.minus.D=c(zandp(samples.A4, Effort.TPccCCDC),
#                         zandp(samples.A4, Effort.TPnnCCDN),
#                         zandp(samples.A4, Effort.TPnnCCDC),
#                         zandp(samples.A4, Effort.TPccCCDN),
#                         "-","-","-","-","-"
#     ),
#     PM.A.minus.B=c(zandp(samples.A4, Effort.TPccPMABC),
#                    zandp(samples.A4, Effort.TPnnPMABN),
#                    zandp(samples.A4, Effort.TPnnPMABC),
#                    zandp(samples.A4, Effort.TPccPMABN),
#
#                    zandp(samples.A4, Effort.TPppPMABP),
#
#                    zandp(samples.A4, Effort.TPpcPMABC),
#                    zandp(samples.A4, Effort.TPpnPMABN),
#                    zandp(samples.A4, Effort.TPpnPMABC),
#                    zandp(samples.A4, Effort.TPpcPMABN)
#
#     ),
#     PM.B.minus.C=c(zandp(samples.A4, Effort.TPccPMBCC),
#                    zandp(samples.A4, Effort.TPnnPMBCN),
#                    zandp(samples.A4, Effort.TPnnPMBCC),
#                    zandp(samples.A4, Effort.TPccPMBCN),
#
#                    zandp(samples.A4, Effort.TPppPMBCP),
#
#                    zandp(samples.A4, Effort.TPpcPMBCC),
#                    zandp(samples.A4, Effort.TPpnPMBCN),
#                    zandp(samples.A4, Effort.TPpnPMBCC),
#                    zandp(samples.A4, Effort.TPpcPMBCN)
#     ),
#     PM.C.minus.D=c(zandp(samples.A4, Effort.TPccPMCDC),
#                    zandp(samples.A4, Effort.TPnnPMCDN),
#                    zandp(samples.A4, Effort.TPnnPMCDC),
#                    zandp(samples.A4, Effort.TPccPMCDN),
#
#                    zandp(samples.A4, Effort.TPppPMCDP),
#
#                    zandp(samples.A4, Effort.TPpcPMCDC),
#                    zandp(samples.A4, Effort.TPpnPMCDN),
#                    zandp(samples.A4, Effort.TPpnPMCDC),
#                    zandp(samples.A4, Effort.TPpcPMCDN)
#     )
# ))
#
# colnames(Effort.by.TP.Table.Long) <- c("Conflict","Nonconflict","Conflict (FA)","Nonconflict (FA)","PM","Conflict to PMC","Nonconflict to PMC","Conflict to PMN","Nonconflict to PMN")
# Effort.by.TP.Table.Long

#
# Effort.by.TP.Table.PM <- data.frame(rbind(
#     PM.to.PM=c(zandp(samples.A4, Effort.TPppPMABP),
#                zandp(samples.A4, Effort.TPppPMBCP),
#                zandp(samples.A4, Effort.TPppPMCDP),
#                "-","-","-"
#     ),
#
#     Conflict.to.PM=c(zandp(samples.A4, Effort.TPpcPMABC),
#                      zandp(samples.A4, Effort.TPpcPMBCC),
#                      zandp(samples.A4, Effort.TPpcPMCDC),
#
#                      zandp(samples.A4, Effort.TPpcPMABN),
#                      zandp(samples.A4, Effort.TPpcPMBCN),
#                      zandp(samples.A4, Effort.TPpcPMCDN)
#     ),
#
#     Nonconflict.to.PM=c(zandp(samples.A4, Effort.TPpnPMABN),
#                         zandp(samples.A4, Effort.TPpnPMBCN),
#                         zandp(samples.A4, Effort.TPpnPMCDN),
#
#                         zandp(samples.A4, Effort.TPpnPMABC),
#                         zandp(samples.A4, Effort.TPpnPMBCC),
#                         zandp(samples.A4, Effort.TPpnPMCDC)
#
#     )
# ))
# colnames(Effort.by.TP.Table.PM) <- c("A-B (Correct)","B-C (Correct)","C-D (Correct)","A-B (FA)","B-C (FA)","C-D (FA)")
# Effort.by.TP.Table.PM


Reactive.Table
Proactive.Table
Capacity.Table
ProactiveTP.Table
Effort.by.TP.Table.Wide
# Effort.by.TP.Table.Long
# Effort.by.TP.Table.PM

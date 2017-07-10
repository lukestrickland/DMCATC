# PMDC manuscript.
rm(list=ls())
# setwd("~/russ_model_analyses")
setwd("D:/Software/DMC_ATCPMDC")
# setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")


load("data/samples/A1.block.B.V_cond.B.V.PMV.samples.RData")
samples.A1 <- A1.block.B.V_cond.B.V.PMV.samples 
rm(A1.block.B.V_cond.B.V.PMV.samples)
load("data/samples/A2.block.B.V_cond.B.V.PMV.samples.RData")
samples.A2 <- A2.block.B.V_cond.B.V.PMV.samples 
rm(A2.block.B.V_cond.B.V.PMV.samples)
load("data/samples/A3.block.B.V_cond.B.V.PMV.samples.RData")
samples.A3 <- A3.block.B.V_cond.B.V.PMV.samples 
rm(A3.block.B.V_cond.B.V.PMV.samples)
load("data/samples/A4.block.B.V_cond.B.V.PMV.samples.RData")
samples.A4 <- A4.block.B.V_cond.B.V.PMV.samples 
rm(A4.block.B.V_cond.B.V.PMV.samples)

av.B.3P <- function (thetas) (thetas[,"B.A3P",, drop=F] + thetas[,"B.B3P",, drop=F]
                              + thetas[,"B.C3P",, drop=F] + thetas[,"B.D3P",, drop=F])/4
av.B.diff.C <- function (thetas) ((thetas[,"B.A3C",, drop=F] + thetas[,"B.B3C",, drop=F]
                                   + thetas[,"B.C3C",, drop=F] + thetas[,"B.D3C",, drop=F])/4 - 
                                    (thetas[,"B.A2C",, drop=F] + thetas[,"B.B2C",, drop=F] 
                                     + thetas[,"B.C2C",, drop=F] + thetas[,"B.D2C",, drop=F])/4)
av.B.diff.N <- function (thetas) ((thetas[,"B.A3N",, drop=F] + thetas[,"B.B3N",, drop=F] + 
                                     thetas[,"B.C3N",, drop=F] + thetas[,"B.D3N",, drop=F])/4 - 
                                    (thetas[,"B.A2N",, drop=F]
                                  + thetas[,"B.B2N",, drop=F] + thetas[,"B.C2N",, drop=F] +
                                    thetas[,"B.D2N",, drop=F])/4)


#Russ enjoy sorting through that mate 
Z.p.acrossexp(samples.A1, samples.A2, fun= av.B.3P)
Z.p.acrossexp(samples.A2, samples.A3, fun= av.B.3P)
Z.p.acrossexp(samples.A1, samples.A3, fun= av.B.3P)
Z.p.acrossexp(samples.A1, samples.A4, fun= av.B.3P)
Z.p.acrossexp(samples.A1, samples.A2, fun= av.B.diff.C )
Z.p.acrossexp(samples.A1, samples.A3, fun= av.B.diff.C )
Z.p.acrossexp(samples.A2, samples.A3, fun= av.B.diff.C )
Z.p.acrossexp(samples.A1, samples.A4, fun= av.B.diff.C )
Z.p.acrossexp(samples.A2, samples.A4, fun= av.B.diff.C )
Z.p.acrossexp(samples.A1, samples.A2, fun= av.B.diff.N )
Z.p.acrossexp(samples.A2, samples.A3, fun= av.B.diff.N )
Z.p.acrossexp(samples.A1, samples.A3, fun= av.B.diff.N )
Z.p.acrossexp(samples.A1, samples.A4, fun= av.B.diff.N )

#Interactions next? cond x proactive ? 4x as many comparisons as above
# have fun mate :) (nah probably start with a graph and if it looks like nothing
#leave it)



# # # Reactive Control # # #

Reactive.ccC <- function (thetas) ((thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.pcA3C",, drop=F]) +
                                     (thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.pcB3C",, drop=F]) +
                                     (thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.pcC3C",, drop=F]) +
                                     (thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.pcD3C",, drop=F]))/4
Reactive.nnN <- function (thetas) ((thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.pnA3N",, drop=F]) +
                                     (thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.pnB3N",, drop=F]) +
                                     (thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.pnC3N",, drop=F]) +
                                     (thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.pnD3N",, drop=F]))/4
Reactive.ccN <- function (thetas) ((thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.pcA3N",, drop=F]) +
                                     (thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.pcB3N",, drop=F]) +
                                     (thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.pcC3N",, drop=F]) +
                                     (thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.pcD3N",, drop=F]))/4
Reactive.nnC <- function (thetas) ((thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.pnA3C",, drop=F]) +
                                     (thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.pnB3C",, drop=F]) +
                                     (thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.pnC3C",, drop=F]) +
                                     (thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.pnD3C",, drop=F]))/4
#Some more fun for you Russ enjoy!
Z.p.acrossexp(samples.A1, samples.A2, fun= Reactive.ccC)
Z.p.acrossexp(samples.A2, samples.A3, fun= Reactive.ccC)
Z.p.acrossexp(samples.A1, samples.A3, fun= Reactive.ccC)
Z.p.acrossexp(samples.A1, samples.A4, fun= Reactive.ccC)
Z.p.acrossexp(samples.A1, samples.A2, fun= Reactive.nnN )
Z.p.acrossexp(samples.A1, samples.A3, fun= Reactive.nnN )
Z.p.acrossexp(samples.A2, samples.A3, fun= Reactive.nnN )
Z.p.acrossexp(samples.A1, samples.A4, fun= Reactive.nnN )
Z.p.acrossexp(samples.A2, samples.A4, fun= Reactive.nnN )
Z.p.acrossexp(samples.A1, samples.A2, fun= Reactive.ccN )
Z.p.acrossexp(samples.A2, samples.A3, fun= Reactive.ccN )
Z.p.acrossexp(samples.A1, samples.A3, fun= Reactive.ccN )
Z.p.acrossexp(samples.A1, samples.A4, fun= Reactive.ccN )
Z.p.acrossexp(samples.A1, samples.A2, fun= Reactive.nnC )
Z.p.acrossexp(samples.A2, samples.A3, fun= Reactive.nnC )
Z.p.acrossexp(samples.A1, samples.A3, fun= Reactive.nnC )
Z.p.acrossexp(samples.A1, samples.A4, fun= Reactive.nnC )

#Interactions next? cond x Reactive ? 4x as many comparisons as above
# have fun mate :) (nah probably start with a graph and if it looks like nothing
#leave it)

rm(list=ls())
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")

CDT.Acc.S <- read.csv("analysis/CDT.Acc.S.E1.csv")
CDT.Acc.Block <- read.csv("analysis/CDT.Acc.Block.E1.csv")
CDT.Acc.Cond <- read.csv("analysis/CDT.Acc.Cond.E1.csv")

CDT.Acc.S
CDT.Acc.Cond

CDT.RT.S <- read.csv("analysis/CDT.RT.S.E1.csv")
CDT.RT.Block <- read.csv("analysis/CDT.RT.Block.E1.csv")
CDT.RT.Cond <- read.csv("analysis/CDT.RT.Cond.E1.csv")

CDT.RT.S

CDT.Reactive.S <- read.csv("analysis/CDT.Reactive.S.E1.csv")
CDT.Reactive.TABLE <- read.csv("analysis/CDT.Reactive.TABLE.E1.csv")

CDT.Reactive.S
CDT.Reactive.TABLE



sdvs <- read.csv("analysis/sdvs.E1.csv")
sdvs <- data.frame(sdv.table)
colnames(sdvs) <- c("Parameter", "Mean", "SD")

sdv.table <- rbind(
    c("", "Conflict", "Nonconflict", "PM (Conflict)", "PM (Nonconflict)"),
    c("Conflict",
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.ccC" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.ccC" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.nnC" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.nnC" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.pcC" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.pcC" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.pnC" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.pnC" ],2), ")", sep = "")),

    c("Nonconflict",
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.ccN" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.ccN" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.nnN" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.nnN" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.pcN" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.pcN" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.pnN" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.pnN" ],2), ")", sep = "")),

    c("PM", "Fixed at 0.5", "Fixed at 0.5",
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.ppP" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.ppP" ],2), ")", sep = ""),
      paste(round(sdvs$Mean[ sdvs$Parameter=="sdv.ppP" ],2), " (",
            round(sdvs$SD[ sdvs$Parameter=="sdv.ppP" ],2), ")", sep = ""))
)
sdv.table <- data.frame(sdv.table)
colnames(sdv.table) <- c("Accumulator", "", "", "Stimulus Type", "")
sdv.table

pander(sdv.table)

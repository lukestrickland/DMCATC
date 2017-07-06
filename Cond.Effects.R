rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")


load("data/after_sampling/E1PP.RData")

sim <- do.call(rbind, E1PP)
currentsim <- sim
# data <- lapply(E1PP, function(x) attr(x, "data"))
# data <- do.call(rbind, data)

cond.effects <- function (currentsim) {

    RTccCA <- NA;RTccNA <- NA;RTnnNA <- NA;RTnnCA <- NA
    RTccCB <- NA;RTccNB <- NA;RTnnNB <- NA;RTnnCB <- NA
    RTccCC <- NA;RTccNC <- NA;RTnnNC <- NA;RTnnCC <- NA
    RTccCD <- NA;RTccND <- NA;RTnnND <- NA;RTnnCD <- NA

    RTdiffccCAB = NA;RTdiffccNAB = NA
    RTdiffnnNAB = NA;RTdiffnnCAB = NA

    RTdiffccCBC = NA;RTdiffccNBC = NA
    RTdiffnnNBC = NA;RTdiffnnCBC = NA

    RTdiffccCCD = NA;RTdiffccNCD = NA
    RTdiffnnNCD = NA;RTdiffnnCCD = NA

    accCA <- NA;accCB <- NA;accCC <- NA;accCD <- NA
    accNA <- NA;accNB <- NA;accNC <- NA;accND <- NA

    accdiffCAB = NA; accdiffNAB = NA
    accdiffCBC = NA; accdiffNBC = NA
    accdiffCCD = NA; accdiffNCD = NA

    pmaccA <- NA; pmaccB <- NA; pmaccC <- NA; pmaccD <- NA
    pmaccdiffAB <- NA; pmaccdiffBC <- NA; pmaccdiffCD <- NA

    pmaccA <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="A"])/
        length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="A"])
    pmaccB <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="B"])/
        length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="B"])
    pmaccC <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="C"])/
        length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="C"])
    pmaccD <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="D"])/
        length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="D"])

    pmaccdiffAB <- pmaccA - pmaccB
    pmaccdiffBC <- pmaccB - pmaccC
    pmaccdiffCD <- pmaccC - pmaccD

    # RT for each stim by response by cond
    RTccCA <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])
    RTccNA <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="A"])
    RTnnNA <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])
    RTnnCA <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="A"])

    RTccCB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])
    RTccNB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"])
    RTnnNB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])
    RTnnCB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"])

    RTccCC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])
    RTccNC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"])
    RTnnNC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])
    RTnnCC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"])

    RTccCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])
    RTccND <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="D"])
    RTnnND <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])
    RTnnCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="D"])

    # RT differences between time pressure conditions for each stim by response
    RTdiffccCAB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])
    RTdiffccNAB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="A"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"])
    RTdiffnnNAB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])
    RTdiffnnCAB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="A"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"])

    RTdiffccCBC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])
    RTdiffccNBC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"])
    RTdiffnnNBC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])
    RTdiffnnCBC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"])

    RTdiffccCCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])
    RTdiffccNCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="D"])
    RTdiffnnNCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])
    RTdiffnnCCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="D"])

    # Accuracy by stimulus by condition
    accCA <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="A"])
    accCB <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"])
    accCC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"])
    accCD <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="D"])

    accNA <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="A"])
    accNB <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"])
    accNC <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"])
    accND <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="D"])

    # Accuracy differences between time pressure conditions for each stimulus
    accdiffCAB <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="A"]) -
        length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"])
    accdiffNAB <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="A"]) -
        length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"])

    accdiffCBC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"]) -
        length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"])
    accdiffNBC <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"]) -
        length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"])

    accdiffCCD <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"]) -
        length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="D"])
    accdiffNCD <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"]) -
        length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="D"])

    out <- c(pmaccA,pmaccB,pmaccC,pmaccD,

             pmaccdiffAB,pmaccdiffBC,pmaccdiffCD,

             RTccCA,RTccCB,RTccCC,RTccCD,
             RTnnNA,RTnnNB,RTnnNC,RTnnND,
             RTnnCA,RTnnCB,RTnnCC,RTnnCD,
             RTccNA,RTccNB,RTccNC,RTccND,

             RTdiffccCAB,RTdiffccCBC,RTdiffccCCD,
             RTdiffnnNAB,RTdiffnnNBC,RTdiffnnNCD,
             RTdiffnnCAB,RTdiffnnCBC,RTdiffnnCCD,
             RTdiffccNAB,RTdiffccNBC,RTdiffccNCD,

             accCA,accCB,accCC,accCD,
             accNA,accNB,accNC,accND,

             accdiffCAB,accdiffCBC,accdiffCCD,
             accdiffNAB,accdiffNBC,accdiffNCD

             # pmrtdiff,

             )

    names(out) <- c("PM Accuracy A","PM Accuracy B","PM Accuracy C","PM Accuracy D",
                    "PM Acc Diff A-B","PM Acc Diff B-C","PM Acc Diff C-D",

                    "RT Conflict A","RT Conflict B","RT Conflict C","RT Conflict D",
                    "RT Nonconflict A","RT Nonconflict B","RT Nonconflict C","RT Nonconflict D",
                    "RT Conflict (FA) A","RT Conflict (FA) B","RT Conflict (FA) C","RT Conflict (FA) D",
                    "RT Nonconflict (FA) A","RT Nonconflict (FA) B","RT Nonconflict (FA) C","RT Nonconflict (FA) D",

                    "RT Diff Conflict A-B","RT Diff Conflict B-C","RT Diff Conflict C-D",
                    "RT Diff Nonconflict A-B","RT Diff Nonconflict B-C","RT Diff Nonconflict C-D",
                    "RT Diff Conflict (FA) A-B","RT Diff Conflict (FA) B-C","RT Diff Conflict (FA) C-D",
                    "RT Diff Nonconflict (FA) A-B","RT Diff Nonconflict (FA) B-C","RT Diff Nonconflict (FA) C-D",

                    "Accuracy Conflict A","Accuracy Conflict B","Accuracy Conflict C","Accuracy Conflict D",
                    "Accuracy Nonconflict A","Accuracy Nonconflict B","Accuracy Nonconflict C","Accuracy Nonconflict D",

                    "Acc Diff Conflict A-B","Acc Diff Conflict B-C","Acc Diff Conflict C-D",
                    "Acc Diff Nonconflict A-B","Acc Diff Nonconflict B-C","Acc Diff Nonconflict C-D"

                    )
    out

}

Cond.Effects <- cond.effects(currentsim)
Cond.Effects <- data.frame(Cond.Effects)
Cond.Effects


#
load("data/after_sampling/A1PP.RData")
load("data/after_sampling/A2PP.RData")
load("data/after_sampling/A3PP.RData")
load("data/after_sampling/A4PP.RData")

sim.A1 <- do.call(rbind, A1PP)
sim.A2 <- do.call(rbind, A2PP)
sim.A3 <- do.call(rbind, A3PP)
sim.A4 <- do.call(rbind, A4PP)

Cond.Effects.A1 <- cond.effects(sim.A1)
Cond.Effects.A2 <- cond.effects(sim.A2)
Cond.Effects.A3 <- cond.effects(sim.A3)
Cond.Effects.A4 <- cond.effects(sim.A4)

data.frame(Cond.Effects.A1, Cond.Effects.A2, Cond.Effects.A3, Cond.Effects.A4)

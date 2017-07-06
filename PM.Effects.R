

rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")




# E1PP <- h.post.predict.dmc(samples, save.simulation=
#                                TRUE)
# save(E1PP, file="data/after_sampling/E1PP.RData")
load("data/after_sampling/E1PP.RData")

sim <- do.call(rbind, E1PP)

currentsim <- sim

x1pmeffects <- function (currentsim) {

    costccC = NA;costccN = NA
    costnnN = NA;costnnC = NA
    accC = NA; accN = NA
    # nonaccC = NA; nonaccN = NA
    pmacc <- NA

    pmacc <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P"])/
        length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn")])

    costccC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])

    costccN <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="3"]) -
        mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="2"])

    costnnN <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])

    costnnC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="3"]) -
        mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="2"])

    accC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="3"]) -
        length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])/
        length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="2"])

    accN <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="3"]) -
        length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])/
        length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="2"])


    out <- c(pmacc,
             # pmrtdiff,
             costccC,costnnC,
             costnnN,costccN,
             # noncostccC,noncostccN,
             # noncostnnN,noncostnnC,
             accC,
             # nonaccC,
             accN
             # nonaccN
             )

    names(out) <- c("PM Accuracy",
                    # "pmRTdiff",
                    "RT Cost Conflict","RT Cost Conflict (FA)",
                    "RT Cost Nonconflict","RT Cost Nonconflict (FA)",
                    # "noncostccC","noncostccN",
                    # "noncostnnN","noncostnnC",
                    "Accuracy Cost Conflict",
                    # "nonaccC",
                    "Accuracy Cost Nonconflict"
                    # "nonaccN"
                    )
    out

}

PM.Effects <- x1pmeffects(currentsim)
PM.Effects <- data.frame(PM.Effects)
PM.Effects


#
load("data/after_sampling/A1PP.RData")
load("data/after_sampling/A2PP.RData")
load("data/after_sampling/A3PP.RData")
load("data/after_sampling/A4PP.RData")

sim.A1 <- do.call(rbind, A1PP)
sim.A2 <- do.call(rbind, A2PP)
sim.A3 <- do.call(rbind, A3PP)
sim.A4 <- do.call(rbind, A4PP)

PM.Effects.A1 <- x1pmeffects(sim.A1)
PM.Effects.A2 <- x1pmeffects(sim.A2)
PM.Effects.A3 <- x1pmeffects(sim.A3)
PM.Effects.A4 <- x1pmeffects(sim.A4)

data.frame(PM.Effects.A1, PM.Effects.A2, PM.Effects.A3, PM.Effects.A4)

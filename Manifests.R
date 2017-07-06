rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

get.accuracy <- function (data) {

    acc.C <- NA; acc.N <- NA
    acc.P <- NA

    acc.C.control <- NA;acc.C.pm <- NA;
    acc.N.control <- NA;acc.N.pm <- NA;

    acc.C.A <- NA;acc.C.B <- NA;acc.C.C <- NA;acc.C.D <- NA;
    acc.N.A <- NA;acc.N.B <- NA;acc.N.C <- NA;acc.N.D <- NA;
    acc.P.A <- NA;acc.P.B <- NA;acc.P.C <- NA;acc.P.D <- NA;

    acc.C <- length(data$RT[data$S=="cc" & data$R=="C"])/
        length(data$RT[data$S=="cc"])
    acc.N <- length(data$RT[data$S=="nn" & data$R=="N"])/
        length(data$RT[data$S=="nn"])
    acc.P <- length(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P"])/
        length(data$RT[(data$S=="pc"|data$S=="pn")])

    acc.C.control <- length(data$RT[data$S=="cc" & data$R=="C" & data$block=="2"])/
        length(data$RT[data$S=="cc" & data$block=="2"])
    acc.C.pm <- length(data$RT[data$S=="cc" & data$R=="C" & data$block=="3"])/
        length(data$RT[data$S=="cc" & data$block=="3"])

    acc.N.control <- length(data$RT[data$S=="nn" & data$R=="N" & data$block=="2"])/
        length(data$RT[data$S=="nn" & data$block=="2"])
    acc.N.pm <- length(data$RT[data$S=="nn" & data$R=="N" & data$block=="3"])/
        length(data$RT[data$S=="nn" & data$block=="3"])

    acc.C.A <- length(data$RT[data$S=="cc" & data$R=="C" & data$cond=="A"])/
        length(data$RT[data$S=="cc" & data$cond=="A"])
    acc.C.B <- length(data$RT[data$S=="cc" & data$R=="C" & data$cond=="B"])/
        length(data$RT[data$S=="cc" & data$cond=="B"])
    acc.C.C <- length(data$RT[data$S=="cc" & data$R=="C" & data$cond=="C"])/
        length(data$RT[data$S=="cc" & data$cond=="C"])
    acc.C.D <- length(data$RT[data$S=="cc" & data$R=="C" & data$cond=="D"])/
        length(data$RT[data$S=="cc" & data$cond=="D"])

    acc.N.A <- length(data$RT[data$S=="nn" & data$R=="N" & data$cond=="A"])/
        length(data$RT[data$S=="nn" & data$cond=="A"])
    acc.N.B <- length(data$RT[data$S=="nn" & data$R=="N" & data$cond=="B"])/
        length(data$RT[data$S=="nn" & data$cond=="B"])
    acc.N.C <- length(data$RT[data$S=="nn" & data$R=="N" & data$cond=="C"])/
        length(data$RT[data$S=="nn" & data$cond=="C"])
    acc.N.D <- length(data$RT[data$S=="nn" & data$R=="N" & data$cond=="D"])/
        length(data$RT[data$S=="nn" & data$cond=="D"])

    acc.P.A <- length(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="A"])/
        length(data$RT[(data$S=="pc"|data$S=="pn") & data$cond=="A"])
    acc.P.B <- length(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="B"])/
        length(data$RT[(data$S=="pc"|data$S=="pn") & data$cond=="B"])
    acc.P.C <- length(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="C"])/
        length(data$RT[(data$S=="pc"|data$S=="pn") & data$cond=="C"])
    acc.P.D <- length(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="D"])/
        length(data$RT[(data$S=="pc"|data$S=="pn") & data$cond=="D"])

    out <- c(acc.C, acc.N,
             acc.P,

             acc.C.control,acc.C.pm,
             acc.N.control,acc.N.pm,

             acc.C.A,acc.C.B,acc.C.C,acc.C.D,
             acc.N.A,acc.N.B,acc.N.C,acc.N.D,
             acc.P.A,acc.P.B,acc.P.C,acc.P.D

             )

    names(out) <- c("Conflict","Nonconflict",
                    "PM",

                    "Conflict Control","Conflict PM",
                    "Nonconflict Control","Nonconflict PM",

                    "Conflict A","Conflict B","Conflict C","Conflict D",
                    "Nonconflict A","Nonconflict B","Nonconflict C","Nonconflict D",
                    "PM A","PM B","PM C","PM D"
                    )
    out

}


get.mean.RT <- function (data) {

    RT.C <- NA; RT.N <- NA
    RT.P <- NA

    RT.C.control <- NA;RT.C.pm <- NA;
    RT.N.control <- NA;RT.N.pm <- NA;

    RT.C.A <- NA;RT.C.B <- NA;RT.C.C <- NA;RT.C.D <- NA;
    RT.N.A <- NA;RT.N.B <- NA;RT.N.C <- NA;RT.N.D <- NA;
    RT.P.A <- NA;RT.P.B <- NA;RT.P.C <- NA;RT.P.D <- NA;

    RT.C <- mean(data$RT[data$S=="cc" & data$R=="C"])
    RT.N <- mean(data$RT[data$S=="nn" & data$R=="N"])
    RT.P <- mean(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P"])

    RT.C.control <- mean(data$RT[data$S=="cc" & data$R=="C" & data$block=="2"])
    RT.C.pm <- mean(data$RT[data$S=="cc" & data$R=="C" & data$block=="3"])

    RT.N.control <- mean(data$RT[data$S=="cc" & data$R=="N" & data$block=="2"])
    RT.N.pm <- mean(data$RT[data$S=="cc" & data$R=="N" & data$block=="3"])

    RT.C.A <- mean(data$RT[data$S=="cc" & data$R=="C" & data$cond=="A"])
    RT.C.B <- mean(data$RT[data$S=="cc" & data$R=="C" & data$cond=="B"])
    RT.C.C <- mean(data$RT[data$S=="cc" & data$R=="C" & data$cond=="C"])
    RT.C.D <- mean(data$RT[data$S=="cc" & data$R=="C" & data$cond=="D"])

    RT.N.A <- mean(data$RT[data$S=="nn" & data$R=="N" & data$cond=="A"])
    RT.N.B <- mean(data$RT[data$S=="nn" & data$R=="N" & data$cond=="B"])
    RT.N.C <- mean(data$RT[data$S=="nn" & data$R=="N" & data$cond=="C"])
    RT.N.D <- mean(data$RT[data$S=="nn" & data$R=="N" & data$cond=="D"])

    RT.P.A <- mean(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="A"])
    RT.P.B <- mean(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="B"])
    RT.P.C <- mean(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="C"])
    RT.P.D <- mean(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="D"])


    out <- c(RT.C,RT.N,
             RT.P,

             RT.C.control,RT.C.pm,
             RT.N.control,RT.N.pm,

             RT.C.A,RT.C.B,RT.C.C,RT.C.D,
             RT.N.A,RT.N.B,RT.N.C,RT.N.D,
             RT.P.A,RT.P.B,RT.P.C,RT.P.D
    )

    names(out) <- c("Conflict","Nonconflict",
                    "PM",

                    "Conflict Control","Conflict PM",
                    "Nonconflict Control","Nonconflict PM",

                    "Conflict A","Conflict B","Conflict C","Conflict D",
                    "Nonconflict A","Nonconflict B","Nonconflict C","Nonconflict D",
                    "PM A","PM B","PM C","PM D"
    )
    out

}


get.sd.RT <- function (data) {

    RT.C <- NA; RT.N <- NA
    RT.P <- NA

    RT.C.control <- NA;RT.C.pm <- NA;
    RT.N.control <- NA;RT.N.pm <- NA;

    RT.C.A <- NA;RT.C.B <- NA;RT.C.C <- NA;RT.C.D <- NA;
    RT.N.A <- NA;RT.N.B <- NA;RT.N.C <- NA;RT.N.D <- NA;
    RT.P.A <- NA;RT.P.B <- NA;RT.P.C <- NA;RT.P.D <- NA;

    RT.C <- sd(data$RT[data$S=="cc" & data$R=="C"])
    RT.N <- sd(data$RT[data$S=="nn" & data$R=="N"])
    RT.P <- sd(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P"])

    RT.C.control <- sd(data$RT[data$S=="cc" & data$R=="C" & data$block=="2"])
    RT.C.pm <- sd(data$RT[data$S=="cc" & data$R=="C" & data$block=="3"])

    RT.N.control <- sd(data$RT[data$S=="cc" & data$R=="N" & data$block=="2"])
    RT.N.pm <- sd(data$RT[data$S=="cc" & data$R=="N" & data$block=="3"])

    RT.C.A <- sd(data$RT[data$S=="cc" & data$R=="C" & data$cond=="A"])
    RT.C.B <- sd(data$RT[data$S=="cc" & data$R=="C" & data$cond=="B"])
    RT.C.C <- sd(data$RT[data$S=="cc" & data$R=="C" & data$cond=="C"])
    RT.C.D <- sd(data$RT[data$S=="cc" & data$R=="C" & data$cond=="D"])

    RT.N.A <- sd(data$RT[data$S=="nn" & data$R=="N" & data$cond=="A"])
    RT.N.B <- sd(data$RT[data$S=="nn" & data$R=="N" & data$cond=="B"])
    RT.N.C <- sd(data$RT[data$S=="nn" & data$R=="N" & data$cond=="C"])
    RT.N.D <- sd(data$RT[data$S=="nn" & data$R=="N" & data$cond=="D"])

    RT.P.A <- sd(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="A"])
    RT.P.B <- sd(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="B"])
    RT.P.C <- sd(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="C"])
    RT.P.D <- sd(data$RT[(data$S=="pc"|data$S=="pn") & data$R=="P" & data$cond=="D"])


    out <- c(RT.C,RT.N,
             RT.P,

             RT.C.control,RT.C.pm,
             RT.N.control,RT.N.pm,

             RT.C.A,RT.C.B,RT.C.C,RT.C.D,
             RT.N.A,RT.N.B,RT.N.C,RT.N.D,
             RT.P.A,RT.P.B,RT.P.C,RT.P.D
    )

    names(out) <- c("Conflict","Nonconflict",
                    "PM",

                    "Conflict Control","Conflict PM",
                    "Nonconflict Control","Nonconflict PM",

                    "Conflict A","Conflict B","Conflict C","Conflict D",
                    "Nonconflict A","Nonconflict B","Nonconflict C","Nonconflict D",
                    "PM A","PM B","PM C","PM D"
    )
    out

}


#
load("data/after_sampling/E1PP.RData")

# sim <- do.call(rbind, E1PP)
data.E1 <- lapply(E1PP, function(x) attr(x, "data"))
data.E1 <- do.call(rbind, data.E1)

Accuracy.E1 <- get.accuracy(data.E1)
Accuracy.E1 <- data.frame(Accuracy.E1)
Accuracy.E1

RT.E1 <- get.mean.RT(data.E1)
RT.E1 <- data.frame(RT.E1)
RT.E1

SD.E1 <- get.sd.RT(data.E1)
SD.E1 <- data.frame(SD.E1)
SD.E1
cbind(RT.E1,SD.E1)

ddply(data.E1, ~S*R, summarise, M=mean(RT), SD=sd(RT))
ddply(data.E1, ~S*R*cond, summarise, M=mean(RT), SD=sd(RT))
ddply(data.E1, ~S*R*block, summarise, M=mean(RT), SD=sd(RT))


#
load("data/after_sampling/A1PP.RData")
load("data/after_sampling/A2PP.RData")
load("data/after_sampling/A3PP.RData")
load("data/after_sampling/A4PP.RData")

# sim.A1 <- do.call(rbind, A1PP)
# sim.A2 <- do.call(rbind, A2PP)
# sim.A3 <- do.call(rbind, A3PP)
# sim.A4 <- do.call(rbind, A4PP)

data.A1 <- lapply(A1PP, function(x) attr(x, "data"))
data.A1 <- do.call(rbind, data.A1)
data.A2 <- lapply(A2PP, function(x) attr(x, "data"))
data.A2 <- do.call(rbind, data.A2)
data.A3 <- lapply(A3PP, function(x) attr(x, "data"))
data.A3 <- do.call(rbind, data.A3)
data.A4 <- lapply(A4PP, function(x) attr(x, "data"))
data.A4 <- do.call(rbind, data.A4)

Accuracy.A1 <- get.accuracy(data.A1)
Accuracy.A2 <- get.accuracy(data.A2)
Accuracy.A3 <- get.accuracy(data.A3)
Accuracy.A4 <- get.accuracy(data.A4)

data.frame(Accuracy.A1, Accuracy.A2, Accuracy.A3, Accuracy.A4)



rm(list=ls())
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

# # # Load samples objects
load("~/GitHub/DMCATC/data/samples/A1.block.B.V_cond.B.V.PMV.samples.RData")  # A1 samples object
samples.A1 <- A1.block.B.V_cond.B.V.PMV.samples
rm(A1.block.B.V_cond.B.V.PMV.samples)
rm(samples.A1)

load("~/GitHub/DMCATC/data/samples/A2.block.B.V_cond.B.V.PMV.samples.RData")  # A2 samples object
samples.A2 <- A2.block.B.V_cond.B.V.PMV.samples
rm(A2.block.B.V_cond.B.V.PMV.samples)
rm(samples.A2)

load("~/GitHub/DMCATC/data/samples/A3.block.B.V_cond.B.V.PMV.samples.RData")  # A3 samples object
samples.A3 <- A3.block.B.V_cond.B.V.PMV.samples
rm(A3.block.B.V_cond.B.V.PMV.samples)
rm(samples.A3)

load("~/GitHub/DMCATC/data/samples/A4.block.B.V_cond.B.V.PMV.samples.RData")  # A4 samples object
samples.A4 <- A4.block.B.V_cond.B.V.PMV.samples
rm(A4.block.B.V_cond.B.V.PMV.samples)
rm(samples.A4)
#

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

# # # Functions to get thresholds averaged over condition # # #
#

av.B.2C <- function (thetas) (thetas[,"B.A2C",, drop=F] + thetas[,"B.B2C",, drop=F] + thetas[,"B.C2C",, drop=F] + thetas[,"B.D2C",, drop=F])/4

av.B.2N <- function (thetas) (thetas[,"B.A2N",, drop=F] + thetas[,"B.B2N",, drop=F] + thetas[,"B.C2N",, drop=F] + thetas[,"B.D2N",, drop=F])/4

av.B.3C <- function (thetas) (thetas[,"B.A3C",, drop=F] + thetas[,"B.B3C",, drop=F] + thetas[,"B.C3C",, drop=F] + thetas[,"B.D3C",, drop=F])/4

av.B.3N <- function (thetas) (thetas[,"B.A3N",, drop=F] + thetas[,"B.B3N",, drop=F] + thetas[,"B.C3N",, drop=F] + thetas[,"B.D3N",, drop=F])/4

av.B.3P <- function (thetas) (thetas[,"B.A3P",, drop=F] + thetas[,"B.B3P",, drop=F] + thetas[,"B.C3P",, drop=F] + thetas[,"B.D3P",, drop=F])/4

av.B.diff.C <- function (thetas) ((thetas[,"B.A3C",, drop=F] + thetas[,"B.B3C",, drop=F] + thetas[,"B.C3C",, drop=F] + thetas[,"B.D3C",, drop=F])/4 - (thetas[,"B.A2C",, drop=F] + thetas[,"B.B2C",, drop=F] + thetas[,"B.C2C",, drop=F] + thetas[,"B.D2C",, drop=F])/4)

av.B.diff.N <- function (thetas) ((thetas[,"B.A3N",, drop=F] + thetas[,"B.B3N",, drop=F] + thetas[,"B.C3N",, drop=F] + thetas[,"B.D3N",, drop=F])/4 - (thetas[,"B.A2N",, drop=F] + thetas[,"B.B2N",, drop=F] + thetas[,"B.C2N",, drop=F] + thetas[,"B.D2N",, drop=F])/4)
#

B.diff.C.A1 <- mean.sd(samples.A1,av.B.diff.C)
B.diff.N.A1 <- mean.sd(samples.A1,av.B.diff.N)
B.3P.A1 <- mean.sd(samples.A1,av.B.3P)
av.B.A1 <- data.frame(rbind(B.diff.C.A1=B.diff.C.A1,
                            B.diff.N.A1=B.diff.N.A1,
                            B.3P.A1=B.3P.A1))

B.diff.C.A2 <- mean.sd(samples.A2,av.B.diff.C)
B.diff.N.A2 <- mean.sd(samples.A2,av.B.diff.N)
B.3P.A2 <- mean.sd(samples.A2,av.B.3P)
av.B.A2 <- data.frame(rbind(B.diff.C.A2=B.diff.C.A2,
                            B.diff.N.A2=B.diff.N.A2,
                            B.3P.A2=B.3P.A2))

B.diff.C.A3 <- mean.sd(samples.A3,av.B.diff.C)
B.diff.N.A3 <- mean.sd(samples.A3,av.B.diff.N)
B.3P.A3 <- mean.sd(samples.A3,av.B.3P)
av.B.A3 <- data.frame(rbind(B.diff.C.A3=B.diff.C.A3,
                            B.diff.N.A3=B.diff.N.A3,
                            B.3P.A3=B.3P.A3))

B.diff.C.A4 <- mean.sd(samples.A4,av.B.diff.C)
B.diff.N.A4 <- mean.sd(samples.A4,av.B.diff.N)
B.3P.A4 <- mean.sd(samples.A4,av.B.3P)
av.B.A4 <- data.frame(rbind(B.diff.C.A4=B.diff.C.A4,
                            B.diff.N.A4=B.diff.N.A4,
                            B.3P.A4=B.3P.A4))


Bs <- data.frame(rbind(av.B.A1,av.B.A2,av.B.A3,av.B.A4))
Bs

# # #
# Thresholds data frame for ggplot
Bs$R <- NA; Bs$Emphasis <- NA
Bs$R[grep (".C", rownames(Bs))] <- "Conflict";
Bs$R[grep (".N", rownames(Bs))] <- "Nonconflict";
Bs$R[grep (".P", rownames(Bs))] <- "PM";

Bs$Emphasis[grep ("A1", rownames(Bs))] <- "Neut.";
Bs$Emphasis[grep ("A2", rownames(Bs))] <- "Ongoing";
Bs$Emphasis[grep ("A3", rownames(Bs))] <- "Speed";
Bs$Emphasis[grep ("A4", rownames(Bs))] <- "PM";



plot.df <- Bs
plot.df
plot.df$Emphasis <- factor(plot.df$Emphasis, levels=c("Neut.","Ongoing","Speed","PM"))

B.Ongoing <- ggplot(plot.df[ plot.df$R=="Conflict" | plot.df$R=="Nonconflict",],
            aes(factor(Emphasis),M)) +
    geom_line(aes(group=R, y=M), linetype=2) +
    geom_point(stat = "identity",aes(), size=3) +
    geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
    xlab("Task Emphasis") + ylab("B") +
    scale_shape_discrete("PM Block:") +
    ylim(0.85,1.25) +
    theme(text = element_text()) +
    theme(
        axis.line.x = element_line(),
        axis.line.y = element_line()
    ) +
    ggtitle("Difference in Control and PM Ongoing Task Thresholds by Task Emphasis") +
    facet_grid(.~ R)
B.Ongoing

B.PM <- ggplot(plot.df[ plot.df$R=="PM",],
                    aes(factor(Emphasis),M)) +
    geom_line(aes(group=R, y=M), linetype=2) +
    geom_point(stat = "identity",aes(), size=3) +
    geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
    xlab("Task Emphasis") + ylab("B") +
    scale_shape_discrete("PM Block:") +
    ylim(1.85,2.1) +
    theme(text = element_text()) +
    theme(
        axis.line.x = element_line(),
        axis.line.y = element_line()
    ) +
    ggtitle("PM Threshold by Task Emphasis") +
    facet_grid(.~ R)
B.PM

# # #

# # # Reactive Control # # #
#
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



# > mean(rccc)
# [1] 0.5632861
# > mean(rnnc)
# [1] 0.5342158

# mean(rccc)
# [1] 0.5320913
# > mean(rnnc)
# [1] 0.6097464
#
# mean(rccc)
# [1] 0.5081637
# > mean(rnnc)
# [1] 0.4436646


mean.sd <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    M <- mean(effect)
    SD <- sd(effect)
    data.frame(M, SD)
}

Reactive.ccC.A1 <- mean.sd(samples.A1,Reactive.ccC)
Reactive.nnN.A1 <- mean.sd(samples.A1,Reactive.nnN)
Reactive.nnC.A1 <- mean.sd(samples.A1,Reactive.nnC)
Reactive.ccN.A1 <- mean.sd(samples.A1,Reactive.ccN)
Reac.A1 <- rbind(A1.ccC=Reactive.ccC.A1,A1.nnN=Reactive.nnN.A1,A1.nnC=Reactive.nnC.A1,A1.ccN=Reactive.ccN.A1)
Reac.A1

Reactive.ccC.A2 <- mean.sd(samples.A2,Reactive.ccC)
Reactive.nnN.A2 <- mean.sd(samples.A2,Reactive.nnN)
Reactive.nnC.A2 <- mean.sd(samples.A2,Reactive.nnC)
Reactive.ccN.A2 <- mean.sd(samples.A2,Reactive.ccN)
Reac.A2 <- rbind(A2.ccC=Reactive.ccC.A2,A2.nnN=Reactive.nnN.A2,A2.nnC=Reactive.nnC.A2,A2.ccN=Reactive.ccN.A2)
Reac.A2

Reactive.ccC.A3 <- mean.sd(samples.A3,Reactive.ccC)
Reactive.nnN.A3 <- mean.sd(samples.A3,Reactive.nnN)
Reactive.nnC.A3 <- mean.sd(samples.A3,Reactive.nnC)
Reactive.ccN.A3 <- mean.sd(samples.A3,Reactive.ccN)
Reac.A3 <- rbind(A3.ccC=Reactive.ccC.A3,A3.nnN=Reactive.nnN.A3,A3.nnC=Reactive.nnC.A3,A3.ccN=Reactive.ccN.A3)
Reac.A3

Reactive.ccC.A4 <- mean.sd(samples.A4,Reactive.ccC)
Reactive.nnN.A4 <- mean.sd(samples.A4,Reactive.nnN)
Reactive.nnC.A4 <- mean.sd(samples.A4,Reactive.nnC)
Reactive.ccN.A4 <- mean.sd(samples.A4,Reactive.ccN)
Reac.A4 <- rbind(A4.ccC=Reactive.ccC.A4,A4.nnN=Reactive.nnN.A4,A4.nnC=Reactive.nnC.A4,A4.ccN=Reactive.ccN.A4)
Reac.A4

Mean.Reactive <- data.frame(rbind(Reac.A1,Reac.A2,Reac.A3,Reac.A4))
Mean.Reactive

Mean.Reactive$S <- NA; Mean.Reactive$Emphasis <- NA; Mean.Reactive$R <- NA;

Mean.Reactive$S[grep ("cc", rownames(Mean.Reactive))] <- "Conflict"
Mean.Reactive$S[grep ("nn", rownames(Mean.Reactive))] <- "Nonconflict"
Mean.Reactive$S[grep ("pc", rownames(Mean.Reactive))] <- "PM (Conflict)"
Mean.Reactive$S[grep ("pn", rownames(Mean.Reactive))] <- "PM (Nonconflict)"

Mean.Reactive$R[grep ("ccC", rownames(Mean.Reactive))] <- "Correct"
Mean.Reactive$R[grep ("nnN", rownames(Mean.Reactive))] <- "Correct"
Mean.Reactive$R[grep ("nnC", rownames(Mean.Reactive))] <- "Error"
Mean.Reactive$R[grep ("ccN", rownames(Mean.Reactive))] <- "Error"

Mean.Reactive$Emphasis[grep ("A1", rownames(Mean.Reactive))] <- "Neut."
Mean.Reactive$Emphasis[grep ("A2", rownames(Mean.Reactive))] <- "Ongoing"
Mean.Reactive$Emphasis[grep ("A3", rownames(Mean.Reactive))] <- "Speed"
Mean.Reactive$Emphasis[grep ("A4", rownames(Mean.Reactive))] <- "PM"


plot.reactive <- Mean.Reactive
plot.reactive$Emphasis <- factor(plot.reactive$Emphasis, levels=c("Neut.","Ongoing","Speed","PM"))

ggplot(plot.reactive, aes(factor(Emphasis),M)) +
    geom_line(aes(group=R, y=M), linetype=2) +
    geom_point(stat = "identity", aes(shape=R), size=3) +
    geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
    xlab("Task Emphasis") + ylab("V") +
    scale_shape_discrete("Response:") +
    ylim(0.4,1) +
    theme(text = element_text()) +
    theme(
        axis.line.x = element_line(),
        axis.line.y = element_line()
    ) + ggtitle("Reactive Control by Task Emphasis") +
    facet_grid(. ~ S)

#
# # # Nondecision Time # # #

av.t0 <- function (thetas) thetas[,"t0",, drop=F]

t0.A1 <- mean.sd(samples.A1,av.t0)
t0.A2 <- mean.sd(samples.A2,av.t0)
t0.A3 <- mean.sd(samples.A3,av.t0)
t0.A4 <- mean.sd(samples.A4,av.t0)

t0.A1
t0 <- rbind(A1.t0=t0.A1,A2.t0=t0.A2,A3.t0=t0.A3,A4.t0=t0.A4)
t0

t0$Emphasis <- NA;

t0$Emphasis[grep ("A1", rownames(t0))] <- "Neutral"
t0$Emphasis[grep ("A2", rownames(t0))] <- "Ongoing"
t0$Emphasis[grep ("A3", rownames(t0))] <- "Speed"
t0$Emphasis[grep ("A4", rownames(t0))] <- "PM"

plot.t0 <- t0
plot.t0$Emphasis <- factor(plot.t0$Emphasis, levels=c("Neutral","Ongoing","Speed","PM"))

ggplot(plot.t0, aes(factor(Emphasis),M)) +
    geom_line(aes(y=M), linetype=2) +
    geom_point(stat = "identity", aes(), size=3) +
    geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
    xlab("Task Emphasis") + ylab("t0") +
    scale_shape_discrete("Response:") +
    ylim(0.25,0.4) +
    theme(text = element_text()) +
    theme(
        axis.line.x = element_line(),
        axis.line.y = element_line()
    ) + ggtitle("Nondecision Time by Task Emphasis")

###########


# # # Capacity (Difference in Drift between Control and PM blocks) # # #
#
Capacity.ccC <- function (thetas) ((thetas[,"mean_v.ccA3C",, drop=F] - thetas[,"mean_v.ccA2C",, drop=F]) +
                                       (thetas[,"mean_v.ccB3C",, drop=F] - thetas[,"mean_v.ccB2C",, drop=F]) +
                                       (thetas[,"mean_v.ccC3C",, drop=F] - thetas[,"mean_v.ccC2C",, drop=F]) +
                                       (thetas[,"mean_v.ccD3C",, drop=F] - thetas[,"mean_v.ccD2C",, drop=F]))/4

Capacity.nnN <- function (thetas) ((thetas[,"mean_v.nnA3N",, drop=F] - thetas[,"mean_v.nnA2N",, drop=F]) +
                                       (thetas[,"mean_v.nnB3N",, drop=F] - thetas[,"mean_v.nnB2N",, drop=F]) +
                                       (thetas[,"mean_v.nnC3N",, drop=F] - thetas[,"mean_v.nnC2N",, drop=F]) +
                                       (thetas[,"mean_v.nnD3N",, drop=F] - thetas[,"mean_v.nnD2N",, drop=F]))/4

Capacity.nnC <- function (thetas) ((thetas[,"mean_v.nnA3C",, drop=F] - thetas[,"mean_v.nnA2C",, drop=F]) +
                                       (thetas[,"mean_v.nnB3C",, drop=F] - thetas[,"mean_v.nnB2C",, drop=F]) +
                                       (thetas[,"mean_v.nnC3C",, drop=F] - thetas[,"mean_v.nnC2C",, drop=F]) +
                                       (thetas[,"mean_v.nnD3C",, drop=F] - thetas[,"mean_v.nnD2C",, drop=F]))/4

Capacity.ccN <- function (thetas) ((thetas[,"mean_v.ccA3N",, drop=F] - thetas[,"mean_v.ccA2N",, drop=F]) +
                                       (thetas[,"mean_v.ccB3N",, drop=F] - thetas[,"mean_v.ccB2N",, drop=F]) +
                                       (thetas[,"mean_v.ccC3N",, drop=F] - thetas[,"mean_v.ccC2N",, drop=F]) +
                                       (thetas[,"mean_v.ccD3N",, drop=F] - thetas[,"mean_v.ccD2N",, drop=F]))/4


cap.ccC.A1 <- mean.sd(samples.A1,Capacity.ccC)
cap.nnN.A1 <- mean.sd(samples.A1,Capacity.nnN)
cap.nnC.A1 <- mean.sd(samples.A1,Capacity.nnC)
cap.ccN.A1 <- mean.sd(samples.A1,Capacity.ccN)

cap.ccC.A2 <- mean.sd(samples.A2,Capacity.ccC)
cap.nnN.A2 <- mean.sd(samples.A2,Capacity.nnN)
cap.nnC.A2 <- mean.sd(samples.A2,Capacity.nnC)
cap.ccN.A2 <- mean.sd(samples.A2,Capacity.ccN)

cap.ccC.A3 <- mean.sd(samples.A3,Capacity.ccC)
cap.nnN.A3 <- mean.sd(samples.A3,Capacity.nnN)
cap.nnC.A3 <- mean.sd(samples.A3,Capacity.nnC)
cap.ccN.A3 <- mean.sd(samples.A3,Capacity.ccN)

cap.ccC.A4 <- mean.sd(samples.A4,Capacity.ccC)
cap.nnN.A4 <- mean.sd(samples.A4,Capacity.nnN)
cap.nnC.A4 <- mean.sd(samples.A4,Capacity.nnC)
cap.ccN.A4 <- mean.sd(samples.A4,Capacity.ccN)

capacity <- rbind(A1.cap.ccC=cap.ccC.A1,
                  A1.cap.nnN=cap.nnN.A1,
                  A1.cap.nnC=cap.nnC.A1,
                  A1.cap.ccN=cap.ccN.A1,

                  A2.cap.ccC=cap.ccC.A2,
                  A2.cap.nnN=cap.nnN.A2,
                  A2.cap.nnC=cap.nnC.A2,
                  A2.cap.ccN=cap.ccN.A2,

                  A3.cap.ccC=cap.ccC.A3,
                  A3.cap.nnN=cap.nnN.A3,
                  A3.cap.nnC=cap.nnC.A3,
                  A3.cap.ccN=cap.ccN.A3,

                  A4.cap.ccC=cap.ccC.A4,
                  A4.cap.nnN=cap.nnN.A4,
                  A4.cap.nnC=cap.nnC.A4,
                  A4.cap.ccN=cap.ccN.A4)
capacity

capacity$S <- NA; capacity$R <- NA; capacity$Emphasis <- NA;

capacity$S[grep ("cc", rownames(capacity))] <- "Conflict"
capacity$S[grep ("nn", rownames(capacity))] <- "Nonconflict"

capacity$R[grep ("ccC", rownames(capacity))] <- "Correct"
capacity$R[grep ("nnN", rownames(capacity))] <- "Correct"
capacity$R[grep ("nnC", rownames(capacity))] <- "Error"
capacity$R[grep ("ccN", rownames(capacity))] <- "Error"

capacity$Emphasis[grep ("A1", rownames(capacity))] <- "Neutral"
capacity$Emphasis[grep ("A2", rownames(capacity))] <- "Ongoing"
capacity$Emphasis[grep ("A3", rownames(capacity))] <- "Speed"
capacity$Emphasis[grep ("A4", rownames(capacity))] <- "PM"


plot.capacity <- capacity
plot.capacity$Emphasis <- factor(plot.capacity$Emphasis, levels=c("Neutral","Ongoing","Speed","PM"))

ggplot(plot.capacity, aes(factor(Emphasis),M)) +
    geom_line(aes(group=R, y=M), linetype=2) +
    geom_point(stat = "identity", aes(shape=R), size=3) +
    geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
    xlab("Task Emphasis") + ylab("V") +
    scale_shape_discrete("Response:") +
    ylim(0,0.5) +
    theme(text = element_text()) +
    theme(
        axis.line.x = element_line(),
        axis.line.y = element_line()
    ) + ggtitle("PM-Control Drift Rate Differences (Capacity) by Task Emphasis") +
    facet_grid(. ~ S)

A1=c(zandp(samples.A1, Capacity.ccC),
     zandp(samples.A1, Capacity.nnN),
     zandp(samples.A1, Capacity.nnC),
     zandp(samples.A1, Capacity.ccN))

A2=c(zandp(samples.A2, Capacity.ccC),
     zandp(samples.A2, Capacity.nnN),
     zandp(samples.A2, Capacity.nnC),
     zandp(samples.A2, Capacity.ccN))

A3=c(zandp(samples.A3, Capacity.ccC),
     zandp(samples.A3, Capacity.nnN),
     zandp(samples.A3, Capacity.nnC),
     zandp(samples.A3, Capacity.ccN))

A4=c(zandp(samples.A4, Capacity.ccC),
     zandp(samples.A4, Capacity.nnN),
     zandp(samples.A4, Capacity.nnC),
     zandp(samples.A4, Capacity.ccN))

Capacity.Table <- data.frame(rbind(A1,A2,A3,A4))
colnames(Capacity.Table) <- c("Conflict","Nonconflict","Conflict (Error)","Nonconflict (Error)")
Capacity.Table

# > Capacity.Table
#    Conflict Nonconflict Conflict (Error) Nonconflict (Error)
# A1 17.96(0)    13.95(0)         13.16(0)            17.16(0)
# A2 14.25(0)     8.97(0)         11.62(0)            14.82(0)
# A3 11.57(0)      9.4(0)         13.41(0)            16.11(0)
# A4 18.65(0)     16.9(0)         16.33(0)            21.54(0)


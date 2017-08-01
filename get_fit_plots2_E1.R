

###################### Russ data analysis template ###############################
# Running this top to bottom should correspond to my data analysis from the
# PMDC manuscript.
rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
sapply(pkgs, require, character.only = T) #load

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)

samples <- samples[names(samples) != "p17"]  # Exclude p17 E1 due to no PM responses

# Check how many runs it took to converge
# If any say "FAIL" then it didn't converge

# for(i in 1:length(samples)) {
#   print(attr(samples[[i]], "auto"))
# }


# Check your samples are all the same length for every participant
# for(i in 1:length(samples)) {
#   print(samples[[i]]$nmc)
# }

# Check the gelman diags. <1.1 is enforced by the sampling algorithm
# gelman.diag.E1 <- gelman.diag.dmc(samples)
# save(gelman.diag.E1, file="gelman.diag.E1.RData")
# load("gelman.diag.E1.RData")


# How many parameters? How many chains?
# str(samples[[1]])
# names(samples[[1]])
# samples[[1]]$theta

# Checkout some trace plots for participants and confirm that the samples look
# like mcmc samples should look.
# plot.dmc(samples[[1]])


# Plot the posterior parameter distributions against the priors to make sure
# that you don't have huge prior influence
# plot.dmc(samples[[1]], p.prior= samples[[1]]$p.prior)

#
# E1PP <- h.post.predict.dmc(samples, save.simulation=TRUE, n.post = 100)
# save(E1PP, file="data/after_sampling/E1PP.RData")

load("data/after_sampling/E1PP.RData")

# E1PP <- h.post.predict.dmc(samples, save.simulation=TRUE, n.post = 250)
# save(E1PP, file="data/after_sampling/E1PP.250.RData")
#
# load("data/after_sampling/E1PP.250.RData")


# rm(samples)


sim <- do.call(rbind, E1PP)
# Do the same for the data
data <- lapply(E1PP, function(x) attr(x, "data"))
data <- do.call(rbind, data)

# Levels for graphing
lev.S <- c("Conflict", "Nonconflict", "PM (Conflict)", "PM (Nonconflict)")
lev.PM <- c("Control", "PM")
lev.cond <- c("LL.LT", "LL.HT",
              "HL.LT", "HL.HT")
lev.R <- c("CR", "NR", "PMR")

# Get a gglist as stored in the PP object.
GGLIST <- get.fitgglist.dmc(sim,data)

# clean out the NAs
GGLIST <- lapply(GGLIST, function(x)  x[is.finite(x$data),])

# Re-order to swap S and PM for desired panel order
# (side by side control vs PM, top and bottom for stim type
# xaxis for cond)
pp.obj <- GGLIST[[1]]
pp.obj <- pp.obj[,c(1,3,2,4,5,6,7,8)]

# Re-level pp.obj for plotting purposes (split cond into TP and load factors)
pp.obj$TP <- NA; pp.obj$load <- NA
pp.obj
pp.obj$TP[grep ("A", pp.obj$cond)] <- "Low"
pp.obj$TP[grep ("B", pp.obj$cond)] <- "High"
pp.obj$TP[grep ("C", pp.obj$cond)] <- "Low"
pp.obj$TP[grep ("D", pp.obj$cond)] <- "High"
pp.obj$TP <- factor(pp.obj$TP)
pp.obj$TP <- factor(pp.obj$TP, levels=c("Low","High"))

pp.obj$load[grep ("A", pp.obj$cond)] <- "2 decisions per trial"
pp.obj$load[grep ("B", pp.obj$cond)] <- "2 decisions per trial"
pp.obj$load[grep ("C", pp.obj$cond)] <- "5 decisions per trial"
pp.obj$load[grep ("D", pp.obj$cond)] <- "5 decisions per trial"
pp.obj$load <- factor(pp.obj$load)

pp.obj$cond <- NULL
pp.obj <- pp.obj[ , c(1,8,2,9,3,4,5,6,7)]
pp.obj
str(pp.obj)


# Better factor names
levels(pp.obj$S) <- lev.S
levels(pp.obj$block) <- lev.PM
# levels(pp.obj$cond) <- lev.cond
levels(pp.obj$R) <- lev.R

## Take only the ONGOING accuracies and drop the R column.
ongoing.acc.obj.2 <- pp.obj[pp.obj$load=="2 decisions per trial" & ((pp.obj$S=="Conflict" & pp.obj$R=="CR")|
                            (pp.obj$S=="Nonconflict" & pp.obj$R=="NR")),
                           !(names(pp.obj) %in% "R")]

ongoing.acc.plots.2 <- ggplot.RP.dmc(ongoing.acc.obj.2, xaxis='TP') +
    ylim(0.4,1) +
    ylab("Accuracy") +
    xlab("Time Pressure") +
    ggtitle("Model Fits (Low Trial Load): \nOngoing Task Accuracy by PM Block and Time Pressure")
ongoing.acc.plots.2

ongoing.acc.obj.5 <- pp.obj[pp.obj$load=="5 decisions per trial" & ((pp.obj$S=="Conflict" & pp.obj$R=="CR")|
                                                    (pp.obj$S=="Nonconflict" & pp.obj$R=="NR")),
                            !(names(pp.obj) %in% "R")]

ongoing.acc.plots.5 <- ggplot.RP.dmc(ongoing.acc.obj.5, xaxis='TP') +
    ylim(0.4,1) +
    ylab("Accuracy") +
    xlab("Time Pressure") +
    ggtitle("Model Fits (High Trial Load): \nOngoing Task Accuracy by PM Block and Time Pressure")
ongoing.acc.plots.5

# ggsave("E1.Fits.Acc.Ongoing.png", plot = last_plot())

## Take only the PM accuracies and drop the R column.
PM.acc.obj <- pp.obj[(pp.obj$S=="PM (Conflict)" & pp.obj$R=="PMR")|
                            (pp.obj$S=="PM (Nonconflict)" & pp.obj$R=="PMR"),
                          !(names(pp.obj) %in% "R")]

PM.acc.plots <- ggplot.RP.dmc(PM.acc.obj,panels.ppage=4, xaxis='TP') +
    ylim(0.4,1) +
    ylab("Accuracy") +
    xlab("Time Pressure") +
    ggtitle("Model Fits: \nPM Accuracy by Time Pressure and Trial Load")
PM.acc.plots

# ggsave("E1.Fits.Acc.PM.png", plot = last_plot())

# The grid arrange will depend on your design
E1.acc.plots.Ongoing <- grid.arrange(ongoing.acc.plots.2, ongoing.acc.plots.5,layout_matrix = cbind(
    c(1,1,2,2), c(1,1,2,2)))
E1.acc.plots.PM <- PM.acc.plots

# ggsave("E1.Fits.Acc.png", plot = E1.acc.plots, width = 9, height = 12)

## stolp ycaruccA

# RT plots
# Ongoing Task RT - dump PM, dump false alarms (so rare)
oRT.obj <- GGLIST[[2]][(GGLIST[[2]]$S=="nn"|
                            GGLIST[[2]]$S=="cc") & GGLIST[[2]]$R!="P",]

oRT.obj <- oRT.obj[,c(1,3,2,4,5,6,7,8,9)]
head(oRT.obj)

# Re-level oRT.obj for plotting purposes (split cond into TP and load factors)
oRT.obj$TP <- NA; oRT.obj$load <- NA
# oRT.obj
oRT.obj$TP[grep ("A", oRT.obj$cond)] <- "Low"
oRT.obj$TP[grep ("B", oRT.obj$cond)] <- "High"
oRT.obj$TP[grep ("C", oRT.obj$cond)] <- "Low"
oRT.obj$TP[grep ("D", oRT.obj$cond)] <- "High"
oRT.obj$TP <- factor(oRT.obj$TP)
oRT.obj$TP <- factor(oRT.obj$TP, levels=c("Low","High"))

oRT.obj$load[grep ("A", oRT.obj$cond)] <- "2 decisions per trial"
oRT.obj$load[grep ("B", oRT.obj$cond)] <- "2 decisions per trial"
oRT.obj$load[grep ("C", oRT.obj$cond)] <- "5 decisions per trial"
oRT.obj$load[grep ("D", oRT.obj$cond)] <- "5 decisions per trial"
oRT.obj$load <- factor(oRT.obj$load)

oRT.obj$cond <- NULL
oRT.obj <- oRT.obj[ , c(1,9,2,10,3,4,5,6,7,8)]
oRT.obj
str(oRT.obj)

#
levels(oRT.obj$S) <- lev.S
levels(oRT.obj$block) <- lev.PM
# levels(oRT.obj$cond) <- lev.cond
levels(oRT.obj$R) <- lev.R

#
corr.RT.2 <- oRT.obj[oRT.obj$load=="2 decisions per trial" & ((oRT.obj$S=="Conflict" & oRT.obj$R=="CR")|
                                                                  (oRT.obj$S=="Nonconflict" & oRT.obj$R=="NR")),  ]
# corr.RT.2
corr.RT.2 <- corr.RT.2[,-5]
# corr.RT.2

corr.RT.plots.2 <- ggplot.RT.dmc(corr.RT.2, panels.ppage=4, xaxis="TP") +
    ylim(0,7) +
    xlab("Time Pressure") +
    ylab("RT (s)") +
    ggtitle("Model Fits (Low Trial Load): \nOngoing Task Correct RTs by PM Block and Time Pressure")
corr.RT.plots.2

corr.RT.5 <- oRT.obj[oRT.obj$load=="5 decisions per trial" & ((oRT.obj$S=="Conflict" & oRT.obj$R=="CR")|
                                                                  (oRT.obj$S=="Nonconflict" & oRT.obj$R=="NR")),  ]
# corr.RT.5
corr.RT.5 <- corr.RT.5[,-5]
# corr.RT.5

corr.RT.plots.5 <- ggplot.RT.dmc(corr.RT.5, panels.ppage=4, xaxis="TP") +
    ylim(0,7) +
    xlab("Time Pressure") +
    ylab("RT (s)") +
    ggtitle("Model Fits (High Trial Load): \nOngoing Task Correct RTs by PM Block and Time Pressure")
corr.RT.plots.5

# ggsave("E1.Fits.RT.Correct.png", plot = last_plot())

#
err.RT.2 <- oRT.obj[oRT.obj$load=="2 decisions per trial" & ((oRT.obj$S=="Conflict" & oRT.obj$R=="NR")|
                        (oRT.obj$S=="Nonconflict" & oRT.obj$R=="CR")),  ]
err.RT.2 <- err.RT.2[,-5]

err.RT.plots.2 <- ggplot.RT.dmc(err.RT.2, panels.ppage=4, xaxis="TP") +
    ylim(0,7) +
    xlab("Time Pressure") +
    ylab ("RT (s)") +
    ggtitle("Model Fits (Low Trial Load): \nOngoing Task Error RTs by PM Block and Time Pressure")
err.RT.plots.2

#
err.RT.5 <- oRT.obj[oRT.obj$load=="5 decisions per trial" & ((oRT.obj$S=="Conflict" & oRT.obj$R=="NR")|
                                                                 (oRT.obj$S=="Nonconflict" & oRT.obj$R=="CR")),  ]
err.RT.5 <- err.RT.5[,-5]

err.RT.plots.5 <- ggplot.RT.dmc(err.RT.5, panels.ppage=4, xaxis="TP") +
    ylim(0,7) +
    xlab("Time Pressure") +
    ylab ("RT (s)") +
    ggtitle("Model Fits (High Trial Load): \nOngoing Task Error RTs by PM Block and Time Pressure")
err.RT.plots.5

# ggsave("E1.Fits.RT.Error.png", plot = last_plot())

E1.RT.plots.Ongoing.2 <- grid.arrange(corr.RT.plots.2, err.RT.plots.2, layout_matrix = cbind(
    c(1,1,2,2), c(1,1,2,2)))

E1.RT.plots.Ongoing.5 <- grid.arrange(corr.RT.plots.5, err.RT.plots.5, layout_matrix = cbind(
    c(1,1,2,2), c(1,1,2,2)))

# # #  PM RT Plots # # #
# There is not that much RT data for PM trials.
# so aggregate differently
# and show fit to TOTAL RT rather than separate correct/error.

pRT.obj <- get.fitgglist.dmc(sim, data, noR=T)[[2]]
pRT.obj <- pRT.obj[is.finite(pRT.obj$data),c(1,3,4,5,6,7,8)]
head(pRT.obj)

# Re-level pRT.obj for plotting purposes (split cond into TP and load factors)
pRT.obj$TP <- NA; pRT.obj$load <- NA
# pRT.obj
pRT.obj$TP[grep ("A", pRT.obj$cond)] <- "Low"
pRT.obj$TP[grep ("B", pRT.obj$cond)] <- "High"
pRT.obj$TP[grep ("C", pRT.obj$cond)] <- "Low"
pRT.obj$TP[grep ("D", pRT.obj$cond)] <- "High"
pRT.obj$TP <- factor(pRT.obj$TP)
pRT.obj$TP <- factor(pRT.obj$TP, levels=c("Low","High"))

pRT.obj$load[grep ("A", pRT.obj$cond)] <- "2 decisions per trial"
pRT.obj$load[grep ("B", pRT.obj$cond)] <- "2 decisions per trial"
pRT.obj$load[grep ("C", pRT.obj$cond)] <- "5 decisions per trial"
pRT.obj$load[grep ("D", pRT.obj$cond)] <- "5 decisions per trial"
pRT.obj$load <- factor(pRT.obj$load)

pRT.obj$cond <- NULL
pRT.obj <- pRT.obj[ , c(1,7,8,2,3,4,5,6)]
pRT.obj
str(pRT.obj)

levels(pRT.obj$S) <- lev.S
# levels(pRT.obj$cond) <- lev.cond

pRT.obj <- pRT.obj[pRT.obj$S=="PM (Nonconflict)"|pRT.obj$S=="PM (Conflict)",]

PM.RT.plots <- ggplot.RT.dmc(pRT.obj, xaxis="TP") +
  ylim(0,5) +
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") +
  ggtitle("Model Fits: \nPM RT by Time Pressure and Trial Load")
PM.RT.plots

# ggsave("E1.Fits.RT.PM.png", plot = last_plot())

E1.RT.plots.PM <- PM.RT.plots
# ggsave("E1.Fits.RT.Ongoing.png", plot = E1.RT.plots.Ongoing, width = 9, height = 12)
# ggsave("E1.Fits.RT.PM.png", plot = PM.RTgraph, width = 9, height = 6)
#

# # # Nonresponse Predictions # # #
#
# Do out of sample predictions of non-responses to see whether they are
# consistent with the model.
load("data/exp_data/okdats.E1.NR.RData")

# Exclude p17 E1 due to no PM responses
okdats <- okdats[ okdats$s!="p17",]
okdats$s <- factor(okdats$s)
# str(okdats)

names(okdats)[length(okdats)] <- "trial.pos"
levels(okdats$block)<- c("2", "3")

for (i in 1:length(samples)) {
  data <- okdats[okdats$s==levels(okdats$s)[i],]
  data <-data[,c(2,3,4,5,6,7)]
  attr(samples[[i]], "NRdata") <- data
}

# h.ordermatched.sims <- h.post.predict.dmc.MATCHORDER(samples)
# save(h.ordermatched.sims,file="data/after_sampling/E1NRPP.RData")
load("data/after_sampling/E1NRPP.RData")

NRsim <- do.call(rbind, h.ordermatched.sims)
getNRdata <- lapply(h.ordermatched.sims, function(x) attr(x, "data"))
NRdata <- do.call(rbind, getNRdata)
NRdata <- add.trial.cumsum.data(NRdata)
NRsim <- add.trial.cumsum.sim (NRsim, NRdata)
###

NR.df <- get.NRs.ATCDMC(NRsim, NRdata, miss_fun=get.trials.missed.E1_A4)
NR.df
# Re-level NR.df for plotting purposes (split cond into TP and load factors)
NR.df$TP <- NA; NR.df$load <- NA
# NR.df
NR.df$TP[grep ("A", NR.df$cond)] <- "Low"
NR.df$TP[grep ("B", NR.df$cond)] <- "High"
NR.df$TP[grep ("C", NR.df$cond)] <- "Low"
NR.df$TP[grep ("D", NR.df$cond)] <- "High"
NR.df$TP <- factor(NR.df$TP)
NR.df$TP <- factor(NR.df$TP, levels=c("Low","High"))

NR.df$load[grep ("A", NR.df$cond)] <- "2 decisions per trial"
NR.df$load[grep ("B", NR.df$cond)] <- "2 decisions per trial"
NR.df$load[grep ("C", NR.df$cond)] <- "5 decisions per trial"
NR.df$load[grep ("D", NR.df$cond)] <- "5 decisions per trial"
NR.df$load <- factor(NR.df$load)

NR.df$cond <- NULL
NR.df <- NR.df[ , c(5,6,1,2,3,4)]
# NR.df
# str(NR.df)

# levels(NR.df$cond) <- lev.cond

NR.plot <- ggplot(NR.df, aes(TP, mean))
# NR.plot
NR.plot <- NR.plot + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower),
                                                        width=0.2) +
    geom_point(aes(TP, data), pch=21, size=3,
               colour="black")+geom_line(aes(group = 1, y=data), linetype=2) +
    ylab("Probability of Nonresponse") +
    ylim(0,0.2) +
    xlab("Time Pressure") +
    facet_wrap(~ load, nrow = NULL, ncol = NULL, scales = "fixed") +
    ggtitle("Model Fits: \nPredicted Probability of Nonresponse by Time Pressure and Trial Load")
NR.plot
# ggsave(NR.plot, file="E1.NR.Prob.png", width = 9, height = 4.5)


# ?facet_wrap



## Model Fits to Response Accuracy
plot <- E1.acc.plots.Ongoing
ggsave("figures/E1/E1.Fits.Accuracy.Ongoing.png", plot = plot, height = 13, width = 9)
# plot(plot)

plot <- E1.acc.plots.PM
ggsave("figures/E1/E1.Fits.Accuracy.PM.png", plot = plot, height = 6.5, width = 9)
# plot(plot)


## Model Fits to RT (Ongoing)
plot <- E1.RT.plots.Ongoing.2
ggsave("figures/E1/E1.Fits.RT.Ongoing.2.png", plot = plot, height = 13, width = 9)
# plot(plot)

plot <- E1.RT.plots.Ongoing.5
ggsave("figures/E1/E1.Fits.RT.Ongoing.5.png", plot = plot, height = 13, width = 9)
# plot(plot)


## Model Fits to RT (PM)
plot <- E1.RT.plots.PM
ggsave("figures/E1/E1.Fits.RT.PM.png", plot = plot, height = 6.5, width = 9)
# plot(plot)

## Out of Sample Predicted Nonresponse Proportions
plot <- NR.plot
ggsave("figures/E1/E1.Fits.NR.png", plot = plot, height = 4, width = 9)
# plot(plot)




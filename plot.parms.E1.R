###################### Russ data analysis template ###############################
# Running this top to bottom should correspond to my data analysis from the
# PMDC manuscript.
rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")


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
# E1PP <- h.post.predict.dmc(samples, save.simulation=
#                              TRUE)
# save(E1PP, file="data/after_sampling/E1PP.RData")
# # <<<<<<< HEAD
load("data/after_sampling/E1PP.RData")
# =======
# >>>>>>> df826efb76acd8ad0e1cb10663f8f981493a9556

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
pp.obj<- GGLIST[[1]]
pp.obj <- pp.obj[,c(1,3,2,4,5,6,7,8)]

# Better factor names
levels(pp.obj$S) <- lev.S
levels(pp.obj$block) <- lev.PM
levels(pp.obj$cond) <- lev.cond
levels(pp.obj$R) <- lev.R

## Take only the ONGOING accuracies and drop the R column.
ongoing.acc.obj <- pp.obj[(pp.obj$S=="Conflict" & pp.obj$R=="CR")|
                            (pp.obj$S=="Nonconflict" & pp.obj$R=="NR"),
                           !(names(pp.obj) %in% "R")]

ongoing.acc.plots <- ggplot.RP.dmc(ongoing.acc.obj, xaxis='cond') +
    ylim(0.4,1) +
    ylab("Accuracy") +
    xlab("Time Pressure/Traffic Load") +
    ggtitle("Ongoing Task Accuracy by PM Block and Time Pressure")
ongoing.acc.plots

# ggsave("E1.Fits.Acc.Ongoing.png", plot = last_plot())

## Take only the PM accuracies and drop the R column.
PM.acc.obj <- pp.obj[(pp.obj$S=="PM (Conflict)" & pp.obj$R=="PMR")|
                            (pp.obj$S=="PM (Nonconflict)" & pp.obj$R=="PMR"),
                          !(names(pp.obj) %in% "R")]

PM.acc.plots <- ggplot.RP.dmc(PM.acc.obj,panels.ppage=4, xaxis='cond') +
    ylim(0.4,1) +
    ylab("Accuracy") +
    xlab("Time Pressure/Traffic Load") +
    ggtitle("PM Task Accuracy by Time Pressure")
PM.acc.plots

# ggsave("E1.Fits.Acc.PM.png", plot = last_plot())

# The grid arrange will depend on your design
E1.acc.plots <- grid.arrange(ongoing.acc.plots, PM.acc.plots,layout_matrix = cbind(
  c(1,1,2), c(1,1,2)))

# ggsave("E1.Fits.Acc.png", plot = E1.acc.plots, width = 9, height = 12)

## stolp ycaruccA

# RT plots
# Ongoing Task RT - dump PM, dump false alarms (so rare)
oRT.obj <- GGLIST[[2]][(GGLIST[[2]]$S=="nn"|
                            GGLIST[[2]]$S=="cc") & GGLIST[[2]]$R!="P",]

oRT.obj <- oRT.obj[,c(1,3,2,4,5,6,7,8,9)]
# head(oRT.obj)

levels(oRT.obj$S) <- lev.S
levels(oRT.obj$block) <- lev.PM
levels(oRT.obj$cond) <- lev.cond
levels(oRT.obj$R) <- lev.R

#
corr.RTs <- oRT.obj[(oRT.obj$S=="Conflict" & oRT.obj$R=="CR")|
                  (oRT.obj$S=="Nonconflict" & oRT.obj$R=="NR"),  ]
corr.RTs <- corr.RTs[,-4]
# corr.RTs

corr.RTgraph <- ggplot.RT.dmc(corr.RTs, panels.ppage=4, xaxis="cond") +
  ylim(0,7) +
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") +
  ggtitle("Ongoing Task Correct RTs by PM Block and Time Pressure")
corr.RTgraph

# ggsave("E1.Fits.RT.Correct.png", plot = last_plot())


err.RTs <- oRT.obj[(oRT.obj$S=="Conflict" & oRT.obj$R=="NR")|
                     (oRT.obj$S=="Nonconflict" & oRT.obj$R=="CR"),  ]
err.RTs <- err.RTs[,-4]

err.RTgraph <- ggplot.RT.dmc(err.RTs, panels.ppage=4, xaxis="cond") +
  ylim(0,7) +
  xlab("Time Pressure/Traffic Load") +
  ylab ("RT (s)") +
  ggtitle("Ongoing Task Error RTs by PM Block and Time Pressure")
err.RTgraph

# ggsave("E1.Fits.RT.Error.png", plot = last_plot())

# There is not that much RT data for PM trials.
# so aggregate differently
# and show fit to TOTAL RT rather than separate correct/error.

pRT.obj <-get.fitgglist.dmc(sim, data, noR=T)[[2]]
pRT.obj <- pRT.obj[is.finite(pRT.obj$data),c(1,3,4,5,6,7,8)]
# head(pRT.obj)

levels(pRT.obj$S) <- lev.S
levels(pRT.obj$cond) <- lev.cond

pRT.obj <- pRT.obj[pRT.obj$S=="PM (Nonconflict)"|pRT.obj$S=="PM (Conflict)",]

PM.RTgraph <- ggplot.RT.dmc(pRT.obj, xaxis="cond") +
  ylim(0,5) +
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") +
  ggtitle("PM RT by Time Pressure")
PM.RTgraph

# ggsave("E1.Fits.RT.PM.png", plot = last_plot())



# fix xaxis names
# corr.RTgraph <- corr.RTgraph + xlab("")
# err.RTgraph <- err.RTgraph + xlab("")
# PM.RTgraph <- PM.RTgraph + xlab("")

E1.RT.plots.ongoing <- grid.arrange(corr.RTgraph, err.RTgraph, layout_matrix = cbind(
  c(1,1,2,2), c(1,1,2,2)))

# ggsave("E1.Fits.RT.Ongoing.png", plot = E1.RT.plots.ongoing, width = 9, height = 12)
# ggsave("E1.Fits.RT.PM.png", plot = PM.RTgraph, width = 9, height = 6)
#


# Do out of sample predictions of non-responses to see whether they are
# consistent with the model.
load("data/exp_data/okdats.E1.NR.RData")
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

levels(NR.df$cond) <- lev.cond

NR.plot <- ggplot(NR.df, aes(cond, mean))
NR.plot <- NR.plot + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower),
                                                        width=0.2) +
    geom_point(aes(cond, data), pch=21, size=4,
               colour="black")+geom_line(aes(group = 1, y=data), linetype=2) +
    ylab("Probability of non-response") +
    xlab("Time Pressure/ Traffic Load") +
    ggtitle("Predicted Probability of Nonresponse by Time Pressure")
NR.plot
# ggsave(NR.plot, file="E1.NR.Prob.png", width = 9, height = 4.5)


# # # Parameter Plots # # #
#
#
# #

# # #
fixedeffects.meanthetas <- function(samples){
  ## bring longer thetas down to min nmc by sampling
  nmcs<- sapply(samples, function(x) x$nmc)
  nmc <- min(nmcs)

  for (i in 1:length(samples)) if (nmcs[i] > nmc) samples[[i]]$theta <-
      samples[[i]]$theta[,,sample(1:dim(samples[[i]]$theta)[3], nmc)]
  samps <- lapply(samples, function(x) x["theta"])
  ##
  ## thetas into big array for apply
  samps2 <- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  samps3<- apply(samps2, c(1,2,3), mean)
  ## back to a theta list after applied
  colnames(samps3)<- colnames(samps[[1]]$theta)
  samps5<- list(samps3)
  attributes(samps5)<- attributes(samps[[1]])
  samps5
}
# # #
#
# # setwd("data/after_sampling")
# av.thetas.E1 <- fixedeffects.meanthetas(samples)[[1]]
# save(av.thetas.E1, file="data/after_sampling/av.thetas.E1.PMV.RData")
load("data/after_sampling/av.thetas.E1.RData")

msds <- cbind(apply(av.thetas.E1, 2, mean), apply(av.thetas.E1, 2, sd))
colnames(msds) <- c("M", "SD")
msds <- data.frame(msds)
# msds

# # # Nondecision Time # # #
t0 <- msds[grep("t0", rownames(msds)),]
t0

# # # Plot Thresholds # # #
Bs <- msds[grep("B\\.", rownames(msds)),]
Bs

# Data frame for ggplot
Bs$PM <- NA; Bs$R <- NA; Bs$cond <- NA
Bs$PM[grep ("2", rownames(Bs))] <- "Control"; Bs$PM[grep ("3", rownames(Bs))] <- "PM"
Bs$R[grep ("2C", rownames(Bs))] <- "Conflict"; Bs$R[grep ("3C", rownames(Bs))] <- "Conflict";
Bs$R[grep ("2N", rownames(Bs))] <- "Nonconflict"; Bs$R[grep ("3N", rownames(Bs))] <- "Nonconflict";
Bs$R[grep ("P", rownames(Bs))] <- "PM"
Bs$cond[grep ("\\.A", rownames(Bs))] <- "LL.LT";
Bs$cond[grep ("\\.B", rownames(Bs))] <- "LL.HT";
Bs$cond[grep ("\\.C", rownames(Bs))] <- "HL.LT";
Bs$cond[grep ("\\.D", rownames(Bs))] <- "HL.HT"

plot.df <- Bs
plot.df$cond <- factor(plot.df$cond)
plot.df$cond <- factor(plot.df$cond, levels=c("LL.LT","LL.HT","HL.LT","HL.HT"))
levels(plot.df$cond)

ggplot(plot.df, aes(factor(cond),M)) +
  geom_point(stat = "identity",aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressure/Traffic Load") + ylab("B") +
  scale_shape_discrete("PM Block:") +
  ylim(0.8,3.5) +
  theme(text = element_text()) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Response Threshold by Time Pressure") +
  facet_grid(. ~ R)

# setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
# ggsave("E1.Thresholds.png", plot = last_plot(), width = 9, height = 6)

# # # Plot Rates # # #

mvs <- msds[grep ("mean_v", rownames(msds)),]
mvs <- data.frame(mvs)

# Data frame for ggplot
mvs$S <- NA; mvs$PM <- NA; mvs$cond <- NA; mvs$R <- NA

mvs$S[grep ("nn", rownames(mvs))] <- "Nonconflict"
mvs$S[grep ("cc", rownames(mvs))] <- "Conflict"
mvs$S[grep ("pc", rownames(mvs))] <- "PM (Conflict)"
mvs$S[grep ("pn", rownames(mvs))] <- "PM (Nonconflict)"
mvs$S[grep ("pp", rownames(mvs))] <- "PM"
mvs

mvs$PM[grep ("2", rownames(mvs))] <- "Control"
mvs$PM[grep ("3", rownames(mvs))] <- "PM"
mvs

mvs$cond[grep ("A2", rownames(mvs))] <- "LL.LT"
mvs$cond[grep ("A3", rownames(mvs))] <- "LL.LT"
mvs$cond[grep ("B2", rownames(mvs))] <- "LL.HT"
mvs$cond[grep ("B3", rownames(mvs))] <- "LL.HT"
mvs$cond[grep ("C2", rownames(mvs))] <- "HL.LT"
mvs$cond[grep ("C3", rownames(mvs))] <- "HL.LT"
mvs$cond[grep ("D2", rownames(mvs))] <- "HL.HT"
mvs$cond[grep ("D3", rownames(mvs))] <- "HL.HT"
mvs

mvs$R[grep ("2N", rownames(mvs))] <- "NR"
mvs$R[grep ("3N", rownames(mvs))] <- "NR"
mvs$R[grep ("2C", rownames(mvs))] <- "CR"
mvs$R[grep ("3C", rownames(mvs))] <- "CR"
mvs$R[grep ("P", rownames(mvs))] <- "PMR"
mvs

#
plot.df <- mvs
plot.df$S <- factor(plot.df$S)
plot.df$PM <- factor(plot.df$PM)
plot.df$cond <- factor(plot.df$cond)
plot.df$R <- factor(plot.df$R)

dim(plot.df)
plot.df <- plot.df[-53,]
str(plot.df)
plot.df
levels(plot.df$cond)
plot.df$cond <- factor(plot.df$cond, levels=c("LL.LT","LL.HT","HL.LT","HL.HT"))


plot.corr <- plot.df[((plot.df$S=="Conflict" & plot.df$R=="CR")|(plot.df$S=="Nonconflict" & plot.df$R=="NR")|(plot.df$S=="PM" & plot.df$R=="PMR")),]
plot.corr
plot.corr.noPM <- plot.df[((plot.df$S=="Conflict" & plot.df$R=="CR")|(plot.df$S=="Nonconflict" & plot.df$R=="NR")),]
plot.corr.noPM
plot.FA <- plot.df[((plot.df$S=="Conflict" & plot.df$R=="NR")|(plot.df$S=="Nonconflict" & plot.df$R=="CR")),]
plot.FA
plot.reactive <- plot.df[((plot.df$S=="Conflict" & plot.df$R=="CR")|(plot.df$S=="Nonconflict" & plot.df$R=="NR")|
                              (plot.df$S=="PM (Conflict)" & plot.df$R=="CR")|(plot.df$S=="PM (Nonconflict)" & plot.df$R=="NR")),]
plot.reactive
plot.ongoing <- plot.df[((plot.df$S=="Conflict" & plot.df$R=="CR")|(plot.df$S=="Nonconflict" & plot.df$R=="NR")|
                             (plot.df$S=="Conflict" & plot.df$R=="NR")|(plot.df$S=="Nonconflict" & plot.df$R=="CR")),]
plot.ongoing

V.corr.graph <- ggplot(plot.corr, aes(factor(cond),M)) +
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressure/Traffic Load") + ylab("V") +
  scale_shape_discrete("PM Block:") +
  ylim(0.7,2.7) +
  theme(text = element_text()) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Drift Rates for Ongoing Task and PM Correct Responses by Time Pressure") +
  facet_grid(. ~ S + R)

V.corr.graph

# ggsave("E1.Rates.Correct.png", plot = last_plot())

# V.ongoing.corr.graph <- ggplot(plot.corr.noPM, aes(factor(cond),M)) +
#   geom_point(stat = "identity", aes(shape=PM), size=3) +
#   geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
#   xlab("Time Pressure/Traffic Load") + ylab("V") +
#   scale_shape_discrete("PM Block:") +
#   ylim(0.9,2.1) +
#   theme(text = element_text()) +
#   theme(
#     axis.line.x = element_line(),
#     axis.line.y = element_line()
#   ) + geom_line(aes(group=PM, y=M), linetype=2) +
#   ggtitle("Drift Rates for Ongoing Task Correct Responses by Time Pressure") +
#   facet_grid(. ~ S + R)
#
# V.ongoing.corr.graph
#
# ggsave("E1.Rates.Ongoing.Correct.png", plot = last_plot())

V.ongoing.FA.graph <- ggplot(plot.FA, aes(factor(cond),M)) +
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressure/Traffic Load") + ylab("V") +
  scale_shape_discrete("PM Block:") +
  ylim(-1,2) +
  theme(text = element_text()) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Drift Rates for Ongoing Task False Alarms by Time Pressure") +
  facet_grid(. ~ S + R)

V.ongoing.FA.graph

# ggsave("E1.Rates.Ongoing.FA.png", plot = last_plot())

V.reactive.graph <- ggplot(plot.reactive, aes(factor(cond),M)) +
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressure/Traffic Load") + ylab("V") +
  scale_shape_discrete("PM Block:") +
  ylim(0.2,2.2) +
  theme(text = element_text()) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Reactive Inhibition of Ongoing Task Drift Rates by Time Pressure") +
  facet_grid(. ~ S + R)

V.reactive.graph

# ggsave("E1.Rates.Reactive.Inhibition.png", plot = last_plot())


V.plots <- grid.arrange(V.corr.graph, V.ongoing.FA.graph,V.reactive.graph, layout_matrix = cbind(
    c(1,1,1,1,2,2,2,3,3,3,3), c(1,1,1,1,2,2,2,3,3,3,3)))

ggsave("E1.Rates.png", plot = V.plots, width = 9, height = 12)


#
# # # Posterior Exploration # # #
# # Below script requires the main samples object called camples object and assumes
# that there is a PP object named PP

rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)

samples <- samples[names(samples) != "p17"]  # Exclude p17 E1 due to no PM responses

load("data/after_sampling/E1PP.RData")
PP <- E1PP

#Set this to TRUE after running all the sims once
#and then it will6 load them rather than re-run.
run.before=T
# <<<<<<< HEAD
# =======
Enam="E1"
# >>>>>>> a23a0198bd3449e18d4d112a6de53a3f6a23b83a
source("generate_postexploration_E1.R")
##The above g

###bLOCK EFFECTS with params averaged


plot1 <- ggplot(all_effects_predictives.block, aes(S, mean))
plot1 <- plot1 + facet_grid(DV  ~  model, scales = "free", space = "fixed") + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes(S, data), pch=21, size=4, colour="black") + xlab("Stimulus")+
  theme(text = element_text()) + ylab("Cost (PM - Control)") +
    ggtitle("Model Exploration: \nObserved vs Predicted Ongoing Task Accuracy and RT (PM-Control Contrasts)")
plot1


#cond effects with params averaged
all_effects_predictives.cond$contrast
all_effects_predictives.cond$contrast[ is.na(all_effects_predictives.cond$contrast) ] <- "C-D"

plotC <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Conflict",], aes(contrast, mean))
plotC <- plotC+
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +
  theme(text = element_text()) + ylab("Difference") + xlab("") + labs(title = "Conflict Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted Conflict Accuracy and RT (Time Pressure Contrasts)")
plotC

plotN <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Nonconflict",], aes(contrast, mean))
plotN <- plotN+
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +xlab("Contrast")+
  theme(text = element_text()) +ylab("Difference")+ labs(title = "Nonconflict Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted Nonconflict Accuracy and RT (Time Pressure Contrasts)")
plotN

plotP <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="PM",], aes(contrast, mean))
plotP <- plotP+
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +
  theme(text = element_text()) +ylab("Difference")+xlab("")+ labs(title = "PM Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted PM Accuracy and RT (Time Pressure Contrasts)")
plotP

#some possible grid arrangements

grid.arrange(plotC, plotN, plotP, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                        c(1,1,1,2,2,2),
                                                        c(1,1,1,2,2,2),
                                                        c(3,3,3,NA,NA,NA),
                                                        c(3,3,3,NA,NA,NA),
                                                        c(3,3,3,NA,NA,NA)

                                                        ))

grid.arrange(plotC, plotN, plotP, layout_matrix = cbind(c(1,1,1,1,2,2,2,2,3,3,3),
                                                        c(1,1,1,1,2,2,2,2,3,3,3),
                                                        c(1,1,1,1,2,2,2,2,3,3,3)

))

grid.arrange(plotC, plotN, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2)


))
plotP
##Turning off PM control mechanisms systematically and examining PM accuracy/
#PM RT
PM_control_mechanisms$data

plot2 <- ggplot(PM_control_mechanisms, aes(model, mean))
plot2 <- plot2 + facet_grid(DV  ~ ., scales = "free", space = "fixed") +
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_line(aes(group=1, y=data), linetype=2) +
    xlab("Included Control Mechanisms") +
    ylab("PM Trial Performance") +
    scale_x_discrete(labels=c("Proactive and Reactive",
                              "Only Reactive", "Only Proactive", "No Control")) +
    theme(text = element_text()) +
    ggtitle("Model Exploration: \nObserved vs Predicted PM Accuracy and RT by Cognitive Control Mechanism")
plot2


grid.arrange(plotP, plot2, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2)


))


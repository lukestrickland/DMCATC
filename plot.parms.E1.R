###################### Russ data analysis template ############################### 
# Running this top to bottom should correspond to my data analysis from the 
# PMDC manuscript. 
rm(list=ls())
# setwd("~/russ_model_analyses")
setwd("D:/Software/DMC_ATCPMDC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples

# Check how many runs it took to converge
# If any say "FAIL" then it didn't converge

for(i in 1:length(samples)) {
  print(attr(samples[[i]], "auto"))
}


# Check your samples are all the same length for every participant
for(i in 1:length(samples)) {
  print(samples[[i]]$nmc)
}

# Check the gelman diags. <1.1 is enforced by the sampling algorithm
# gelman.diag.E1 <- gelman.diag.dmc(samples)
# save(gelman.diag.E1, file="gelman.diag.E1.RData")
load("gelman.diag.E1.RData")


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

#Generate PPs that respect the fact some RTs will be truncated by the trial
#deadline. 

### Do some stuff to index the non-responses and get the simulated data as if
# the sim included trial deadlines

load("data/exp_data/okdats.E1.NR.RData")
names(okdats)[length(okdats)] <- "trial.pos"
levels(okdats$block)<- c("2", "3")

for (i in 1:length(samples)) {
  data <- okdats[okdats$s==levels(okdats$s)[i],]
  data <-data[,c(2,3,4,5,6,7)]
  attr(samples[[i]], "NRdata") <- data
}

h.matched.samples <- h.post.predict.dmc.MATCHORDER(samples)

sim <- do.call(rbind, h.matched.samples)
getdata <- lapply(h.matched.samples, function(x) attr(x, "data"))
data <- do.call(rbind, getdata)
data<-add.trial.cumsum.data(data)
sim<-add.trial.cumsum.sim (sim, data)

###

#Plot non-resposes

#levels for graphing
lev.S <- c("Conflict", "Nonconflict", "PM (Conflict)", "PM (Nonconflict)")
lev.PM <- c("Control", "PM")
lev.cond <- c("LL.LT", "LL.HT",
              "HL.LT", "HL.HT")
lev.R <- c("CR", "NR", "PMR")

NR.df <- get.NRs.ATCDMC(sim, data, miss_fun=get.trials.missed.E1_A4)

levels(NR.df$cond)<- lev.cond

plot <- ggplot(NR.df, aes(cond, mean)) 
plot + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes(cond, data), pch=21, size=4, 
            colour="black")+geom_line(aes(group = 1, y=data), linetype=2)+
  ylab("Probability of non-response") + xlab("Time Pressure/ Traffic Load")

mark_sim_nrs <- sim$cond=="A" & sim$cumsum>12 |sim$cond=="B" & sim$cumsum>8 |  
  sim$cond=="C" & sim$cumsum>20 |sim$cond=="D" & sim$cumsum>10


#Take the non-responses out and take out the factors we needed only to calc them
sim.noNRs <- sim[!mark_sim_nrs,]
data.withs <- data
data.withs$s <- okdats$s

data<-data[!(data$block=="2" & data$R=="P"),]

data.noNRs <- data.withs[data.withs$R!="M",]
data.noNRs$R <- factor(as.character(data.noNRs$R ))

sim.noNRs <-
  sim.noNRs[,!(names(sim.noNRs) %in% c("trial", "trial.pos", "cumsum"))]
data.noNRs <-
  data.noNRs[,!(names(data.noNRs) %in% c("trial", "trial.pos", "cumsum"))]


data.noNRs <- clean(okdats)

#We removed <0.2s from the data and >3 IQR as we did for the data we estimated
#from
data.noNRs <- data.noNRs[,-length(data.noNRs)]

# Desired level names for graphs of factors
checkdata <- get.hdata.dmc(E1.block.B.V_cond.B.V.PMV.samples)

# Conflict R, Nonconf R, PMR
GGLIST <- get.fitgglist.dmc(sim.noNRs, data.noNRs)
check<- h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples)

ggplot.RT.dmc(check, xaxis="cond")
ggplot.RP.dmc(GGLIST[[1]], xaxis="cond")
ggplot.RT.dmc(GGLIST[[2]], xaxis="cond")


## Accuracy plots
acc.obj <- ggplot.obj[[1]]
acc.obj
str(acc.obj)
# re-order to swap S and PM for desired panel order 
# (side by side control vs PM, top and bottom for stim type
# xaxis for cond)
acc.obj <- acc.obj[,c(1,3,2,4,5,6,7,8)]

# Better factor names
levels(acc.obj$S) <- lev.S
levels(acc.obj$block) <- lev.PM
levels(acc.obj$cond) <- lev.cond
levels(acc.obj$R) <- lev.R

ongoing.acc.plots <- ggplot.RA.dmc(acc.obj[acc.obj$S=="Nonconflict"|
acc.obj$S=="Conflict",], ylim=c(0.5, 1)) + 
  ylab("Accuracy") + 
  xlab("Time Pressure/Traffic Load") + 
  ggtitle("Ongoing Task Accuracy by PM Block and Time Pressure") 
ongoing.acc.plots

ggsave("E1.Fits.Acc.Ongoing.png", plot = last_plot())


PM.acc.plots <- ggplot.RA.dmc(acc.obj[acc.obj$S=="PM (Conflict)"|
                                acc.obj$S=="PM (Nonconflict)",],
              acc.fun=function(x){x$R=="PMR"},           
              panels.ppage=4, ylim=c(0.5, 1)) + 
  ylab("Accuracy") + 
  xlab("Time Pressure/Traffic Load") + 
  ggtitle("PM Task Accuracy by Time Pressure") 
PM.acc.plots

ggsave("E1.Fits.Acc.PM.png", plot = last_plot())

# The grid arrange will depend on your design
# grid.arrange(ongoing.acc.plots, PM.acc.plots,layout_matrix = cbind(
#   c(1,1,2), c(1,1,2)))

## stolp ycaruccA

# RT plots
# Ongoing Task RT - dump PM, dump false alarms (so rare)
oRT.obj <- ggplot.obj[[2]][(ggplot.obj[[2]]$S=="nn"|
                            ggplot.obj[[2]]$S=="cc") & ggplot.obj[[2]]$R!="P",]

oRT.obj <- oRT.obj[,c(1,3,2,4,5,6,7,8,9)]
oRT.obj

levels(oRT.obj$S) <- lev.S
levels(oRT.obj$block) <- lev.PM
levels(oRT.obj$cond) <- lev.cond
levels(oRT.obj$R) <- lev.R

#
corr.RTs <- oRT.obj[(oRT.obj$S=="Conflict" & oRT.obj$R=="CR")|
                  (oRT.obj$S=="Nonconflict" & oRT.obj$R=="NR"),  ]
corr.RTs <- corr.RTs[,-4]
corr.RTs

corr.RTgraph <- ggplot.RT.dmc(corr.RTs, panels.ppage=4, xaxis="cond") + 
  ylim(0.2,6) + 
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") + 
  ggtitle("Ongoing Task Correct RTs by PM Block and Time Pressure") 
corr.RTgraph

ggsave("E1.Fits.RT.Correct.png", plot = last_plot())


err.RTs <- oRT.obj[(oRT.obj$S=="Conflict" & oRT.obj$R=="NR")|
                     (oRT.obj$S=="Nonconflict" & oRT.obj$R=="CR"),  ]
err.RTs <- err.RTs[,-4]

err.RTgraph <- ggplot.RT.dmc(err.RTs, panels.ppage=4, xaxis="cond") + 
  ylim(0.2,7) + 
  xlab("Time Pressure/Traffic Load") +
  ylab ("RT (s)") + 
  ggtitle("Ongoing Task Error RTs by PM Block and Time Pressure") 
err.RTgraph

ggsave("E1.Fits.RT.Error.png", plot = last_plot())

# There is not that much RT data for PM trials.
# so aggregate differently
# and show fit to TOTAL RT rather than separate correct/error.

pRT.obj <-get.fitgglist.dmc(sim, data, noR=T)[[2]]
pRT.obj <- pRT.obj[,c(1,3,4,5,6,7,8)]
pRT.obj

levels(pRT.obj$S) <- lev.S
levels(pRT.obj$cond) <- lev.cond

pRT.obj <- pRT.obj[pRT.obj$S=="PM (Nonconflict)"|pRT.obj$S=="PM (Conflict)",]

PM.RTgraph <- ggplot.RT.dmc(pRT.obj, xaxis="cond") + 
  ylim(0.2,4.5) + 
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") + 
  ggtitle("PM RT by Time Pressure") 
PM.RTgraph

ggsave("E1.Fits.RT.PM.png", plot = last_plot())


#
# fix xaxis names
# corr.RTgraph <- corr.RTgraph + xlab("")
# err.RTgraph <- err.RTgraph + xlab("")
# PM.RTgraph <- PM.RTgraph + xlab("")

# grid.arrange(corr.RTgraph, err.RTgraph, PM.RTgraph, layout_matrix = cbind(
#   c(1,1,2,2,3), c(1,1,2,2,3)))
# 


# # # Parameter Plots # # #
#
#
#

# # #
fixedeffects.meanthetas <- function(samples){
  ## bring longer thetas down to min nmc by sampling
  nmcs<- sapply(samples, function(x) x$nmc)
  nmc <- min(nmcs)
  
  for (i in 1:length(samples)) if (nmcs[i] > nmc) samples[[i]]$theta <- samples[[i]]$theta[,,sample(1:dim(samples[[i]]$theta)[3], nmc)]
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

setwd("~/Modelling")
av.thetas.E1 <- fixedeffects.meanthetas(samples)[[1]]
save(av.thetas.E1, file="av.thetas.E1.PMV.RData")
load("av.thetas.E1.PMV.RData")

msds <- cbind(apply(av.thetas.E1, 2, mean), apply(av.thetas.E1, 2, sd))
colnames(msds) <- c("M", "SD")
msds <- data.frame(msds)
msds

# # # Nondecision Time # # # 
t0 <- msds[grep("t0", rownames(msds)),]
t0

# # # Thresholds # # #
Bs <- msds[grep("B\\.", rownames(msds)),]
Bs

# Data frame for ggplot
Bs$PM <- NA; Bs$R <- NA; Bs$cond <- NA
Bs$PM[grep ("2", rownames(Bs))] <- "Control"; Bs$PM[grep ("3", rownames(Bs))] <- "PM"
Bs$R[grep ("2C", rownames(Bs))] <- "Conflict"; Bs$R[grep ("3C", rownames(Bs))] <- "Conflict"; 
Bs$R[grep ("2N", rownames(Bs))] <- "Nonconflict"; Bs$R[grep ("3N", rownames(Bs))] <- "Nonconflict"; 
Bs$R[grep ("P", rownames(Bs))] <- "PM"
Bs$cond[grep ("\\.A", rownames(Bs))] <- "A"; 
Bs$cond[grep ("\\.B", rownames(Bs))] <- "B"; 
Bs$cond[grep ("\\.C", rownames(Bs))] <- "C"; 
Bs$cond[grep ("\\.D", rownames(Bs))] <- "D"

plot.df <- Bs
plot.df

ggplot(plot.df, aes(factor(cond),M)) + 
  geom_point(stat = "identity",aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressure/Traffic Load") + ylab("B") + 
  scale_shape_discrete("PM Block:") + 
  ylim(0.8,3.5) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) + 
  ggtitle("Response Threshold by Time Pressure") +
  facet_grid(. ~ R)

ggsave("E1.Thresholds.png", plot = last_plot())

# # # Plot Rates # # #

mvs <- msds[grep ("mean_v", rownames(msds)),]
mvs <- data.frame(mvs)
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

mvs$cond[grep ("A2", rownames(mvs))] <- "A"
mvs$cond[grep ("A3", rownames(mvs))] <- "A" 
mvs$cond[grep ("B2", rownames(mvs))] <- "B"
mvs$cond[grep ("B3", rownames(mvs))] <- "B" 
mvs$cond[grep ("C2", rownames(mvs))] <- "C"
mvs$cond[grep ("C3", rownames(mvs))] <- "C"
mvs$cond[grep ("D2", rownames(mvs))] <- "D"
mvs$cond[grep ("D3", rownames(mvs))] <- "D"
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

ggplot(plot.corr, aes(factor(cond),M)) + 
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) + 
  xlab("Time Pressue/Traffic Load") + ylab("V") + 
  scale_shape_discrete("PM Block:") + 
  ylim(0.9,2.7) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Correct Response Drift Rates by Time Pressure") +
  facet_grid(. ~ S + R)

ggsave("E1.Rates.Correct.png", plot = last_plot())

ggplot(plot.corr.noPM, aes(factor(cond),M)) + 
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) + 
  xlab("Time Pressue/Traffic Load") + ylab("V") + 
  scale_shape_discrete("PM Block:") + 
  ylim(0.9,2.1) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Ongoing Task Correct Response Drift Rates by Time Pressure") +
  facet_grid(. ~ S + R)

ggsave("E1.Rates.Ongoing.Correct.png", plot = last_plot())

ggplot(plot.FA, aes(factor(cond),M)) + 
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) + 
  xlab("Time Pressue/Traffic Load") + ylab("V") + 
  scale_shape_discrete("PM Block:") + 
  ylim(-0.5,1.5) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Ongoing Task False Alarm Drift Rates by Time Pressure") +
  facet_grid(. ~ S + R)

ggsave("E1.Rates.Ongoing.FA.png", plot = last_plot())

ggplot(plot.reactive, aes(factor(cond),M)) + 
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) +
  xlab("Time Pressue/Traffic Load") + ylab("V") + 
  scale_shape_discrete("PM Block:") + 
  ylim(0.4,2.1) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  ggtitle("Reactive Inhibition of Ongoing Task Drift Rates by Time Pressure") +
  facet_grid(. ~ S + R)

ggsave("E1.Rates.Reactive.Inhibition.png", plot = last_plot())

ggplot(plot.ongoing, aes(factor(cond),M)) + 
  geom_point(stat = "identity", aes(shape=PM), size=3) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.2)) + 
  xlab("Time Pressue Condition") + ylab("V") + 
  scale_shape_discrete("PM Block:") + 
  ylim(-0.5,2.2) +
  theme(text = element_text(size=20)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  ) + geom_line(aes(group=PM, y=M), linetype=2) +
  facet_grid(. ~ S + R) +
  ggtitle("Ongoing Task Drift Rates by Time Pressure") 

ggsave("E1.Rates.Ongoing.Task.png", plot = last_plot())




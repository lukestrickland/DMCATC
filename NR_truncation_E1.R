###################### Russ data analysis template ###############################
# Running this top to bottom should correspond to my data analysis from the
# PMDC manuscript.
rm(list=ls())
# setwd("~/russ_model_analyses")
setwd("D:/Software/DMC_ATCPMDC")
# setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")


#try a few diffrent models with out of sample fitting to see
# if we can detect any over-fitting of the chosen model this way
load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples

load("data/samples/E1.sdv.samples.RData")
samples <- E1.1sdv.samples

load("data/samples/block.B.V_cond.B.V_ABC_D.PMFA.samples.RData")
samples <- block.B.V_cond.B.V_ABC_D.PMFA.samples

# # Do out of sample predictions of non-responses to see whether they are
# # consistent with the model.
load("data/exp_data/okdats.E1.NR.RData")

# load("data/exp_data/okdats.A1.NR.RData")

#remove TP from old okdats obj
# okdats <- okdats[, !(names(okdats) %in% "TP")]

names(okdats)[length(okdats)] <- "trial.pos"
levels(okdats$block)<- c("2", "3")

for (i in 1:length(samples)) {
  data <- okdats[okdats$s==levels(okdats$s)[i],]
  data <-data[,c(2,3,4,5,6,7)]
  attr(samples[[i]], "NRdata") <- data
}
# 

#simulate with order matched to the relevant data frames.
h.ordermatched.sims <- h.post.predict.dmc.MATCHORDER(samples)
save(h.ordermatched.sims,
     file="data/after_sampling/E1NRPP.RData")
load("data/after_sampling/E1NRPP.RData")


##get cumulative sum of RTs for sim and data
NRsim <- do.call(rbind, h.ordermatched.sims)
getNRdata <- lapply(h.ordermatched.sims, function(x) attr(x, "data"))
NRdata <- do.call(rbind, getNRdata)
NRdata <- add.trial.cumsum.data(NRdata)
#the nrdata object is fed below just as a convenient (crappy) way to index 
#the trials for the simulated data frame.
NRsim <- add.trial.cumsum.sim (NRsim, NRdata)
###

# A 12 b 8 c 20 d 10

#add to the data the subject names so that 
NRdata <- cbind(NRdata, rep(names(h.ordermatched.sims), each=1280))
NRdata <- NRdata[!is.na(NRdata$RT),]
NRdata <- NRdata[!(NRdata$R=="M"),]
NRdata[NRdata$cumsum>12 & NRdata$cond=="A",]
NRdata[NRdata$cumsum>8 & NRdata$cond=="B",]
NRdata[NRdata$cumsum>20 & NRdata$cond=="C",]
NRdata[NRdata$cumsum>10 & NRdata$cond=="D",]
names(NRdata)[length(NRdata)] <- "s"



colMeans(get.trials.missed.E1_A4(NRsim)
sim.NAs <- (NRsim$cumsum>12 & NRsim$cond=="A") | (NRsim$cumsum>8 & NRsim$cond=="B")|
  (NRsim$cumsum>20 & NRsim$cond=="C") | (NRsim$cumsum>10 & NRsim$cond=="D")

sim <- NRsim[!sim.NAs,]
colMeans(get.trials.missed.E1_A4(sim))
         
## It turns out this is just getting the original data frame
# so i didn't need to do it and we could just take it out from the samples object
#like so:
# check <- lapply(E1PP, function(x) attr(x, "data"))
# check <- do.call(rbind, check)
# 

#anyway..
clean <- function(df) {
  dfc <- df
  n=tapply(df$RT,list(df$s),length)
  ns=tapply(df$RT,list(df$s),length)
  mn=tapply(df$RT,list(df$s),mean)
  sd=tapply(df$RT,list(df$s),IQR)
  upper <- mn+3*(sd/1.349)
  lower <- 0.2
  bad <- logical(dim(df)[1])
  levs <- paste(df$s,sep=".")
  for (i in levels(df$s)){
    lev <- i
    bad[levs==lev] <- df[levs==lev,"RT"] > upper[i] | df[levs==lev,"RT"] < lower
  }
  df=df[!bad,]
  nok=tapply(df$RT,list(df$s),length)
  pbad=100-100*nok/n
  nok=tapply(df$RT,list(df$s),length)
  pbad=100-100*nok/ns
  print(sort(round(pbad,5)))
  print(mean(pbad,na.rm=T))
  df
}

data<- clean(NRdata)
data<-data[,-length(data)]
data$R <- factor(as.character(data$R))

sim <- sim[,!colnames(sim)  %in% c("trial.pos", "cumsum", "trial")]
data <- data[,!colnames(data)  %in% c("trial.pos", "cumsum", "trial")]

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

# #ggsave("E1.Fits.Acc.Ongoing.png", plot = last_plot())

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

# #ggsave("E1.Fits.Acc.PM.png", plot = last_plot())

# The grid arrange will depend on your design
E1.acc.plots <- grid.arrange(ongoing.acc.plots, PM.acc.plots,layout_matrix = cbind(
  c(1,1,2), c(1,1,2)))

#ggsave("E1.Fits.Acc.png", plot = E1.acc.plots, width = 9, height = 12)

## stolp ycaruccA

# RT plots
# Ongoing Task RT - dump PM, dump false alarms (so rare)
oRT.obj <- GGLIST[[2]][(GGLIST[[2]]$S=="nn"|
                          GGLIST[[2]]$S=="cc") & GGLIST[[2]]$R!="P",]

oRT.obj <- oRT.obj[,c(1,3,2,4,5,6,7,8,9)]
head(oRT.obj)

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
  ylim(0,7) +
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") +
  ggtitle("Ongoing Task Correct RTs by PM Block and Time Pressure")
corr.RTgraph

ggsave("E1.Fits.RT.Correct.OOS.png", plot = last_plot())


err.RTs <- oRT.obj[(oRT.obj$S=="Conflict" & oRT.obj$R=="NR")|
                     (oRT.obj$S=="Nonconflict" & oRT.obj$R=="CR"),  ]
err.RTs <- err.RTs[,-4]

err.RTgraph <- ggplot.RT.dmc(err.RTs, panels.ppage=4, xaxis="cond") +
  ylim(0,7) +
  xlab("Time Pressure/Traffic Load") +
  ylab ("RT (s)") +
  ggtitle("Ongoing Task Error RTs by PM Block and Time Pressure")
err.RTgraph

ggsave("E1.Fits.RT.Error.OOS.png", plot = last_plot())

# There is not that much RT data for PM trials.
# so aggregate differently
# and show fit to TOTAL RT rather than separate correct/error.

pRT.obj <-get.fitgglist.dmc(sim, data, noR=T)[[2]]
pRT.obj <- pRT.obj[is.finite(pRT.obj$data),c(1,3,4,5,6,7,8)]
head(pRT.obj)

levels(pRT.obj$S) <- lev.S
levels(pRT.obj$cond) <- lev.cond

pRT.obj <- pRT.obj[pRT.obj$S=="PM (Nonconflict)"|pRT.obj$S=="PM (Conflict)",]

PM.RTgraph <- ggplot.RT.dmc(pRT.obj, xaxis="cond") +
  ylim(0,5) +
  xlab("Time Pressure/Traffic Load") +
  ylab("RT (s)") +
  ggtitle("PM RT by Time Pressure")
PM.RTgraph

# #ggsave("E1.Fits.RT.PM.png", plot = last_plot())



# fix xaxis names
# corr.RTgraph <- corr.RTgraph + xlab("")
# err.RTgraph <- err.RTgraph + xlab("")
# PM.RTgraph <- PM.RTgraph + xlab("")

E1.RT.plots <- grid.arrange(corr.RTgraph, err.RTgraph, PM.RTgraph, layout_matrix = cbind(
  c(1,1,2,2,3), c(1,1,2,2,3)))

#ggsave("E1.Fits.RT.png", plot = E1.RT.plots, width = 9, height = 12)
#





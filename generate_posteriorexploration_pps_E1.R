
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples

cat("Running a lot of sims, go do something else for a while")

E1PP <- h.post.predict.dmc(samples, save.simulation=
                             TRUE)
save(E1PP, file="data/after_sampling/E1PP.RData")

### First generate posterior predictives
#with certain parameters averaged across conditions.
#avps.post.predict.dmc & avps.h.post.predict.dmc
# grep match all parameters to the av.posts vector,
#and ones that match each other are averaged.


#Aggregate over BLOCK
av.posts.thresblock <- c(
  "^B.A*C", "^B.B*C", "^B.C*C", "^B.D*C", "^B.A*N" ,"^B.B*N",
  "^B.C*N" ,"^B.D*N"
)

av.posts.thresblock <- glob2rx(av.posts.thresblock)

#note this function converts your av.posts to regex using glob2rx
fixed.thresholds.block <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.thresblock)

av.posts.ratesblock<- c(
  "mean_v.ccA*C", "mean_v.nnA*C", "mean_v.ccB*C",
  "mean_v.nnB*C", "mean_v.ccC*C", "mean_v.nnC*C",
  "mean_v.ccD*C", "mean_v.nnD*C", "mean_v.ccA*N",
  "mean_v.nnA*N", "mean_v.ccB*N", "mean_v.nnB*N",
  "mean_v.ccC*N" ,"mean_v.nnC*N", "mean_v.ccD*N",
  "mean_v.nnD*N"
  
)

av.posts.ratesblock <- glob2rx(av.posts.ratesblock)
fixed.rates.block <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.ratesblock)

save(fixed.thresholds.block, fixed.rates.block, file="data/after_sampling/av.posts.block.PPs.RData")



#Aggregate over COND
av.posts.thresholdscond  <- c(
  "B.*2C", "B.*3C", "B.*2N" ,"B.*3N", "B.*3P"
) 

av.posts.thresholdscond <- unique(av.posts.thresholdscond)
av.posts.thresholdscond <- glob2rx(av.posts.thresholdscond)

fixed.thresholds.cond <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.thresholdscond )


av.posts.ratescond<- c("mean_v.cc*2C", "mean_v.nn*2C", "mean_v.cc*3C", "mean_v.nn*3C",
                       "mean_v.pc*3C", "mean_v.pn*3C", "mean_v.cc*2N", "mean_v.nn*2N",
                       "mean_v.cc*3N", "mean_v.nn*3N", "mean_v.pc*3N", "mean_v.pn*3N",
                       "mean_v.pp*3P" )

av.posts.ratescond <- glob2rx(av.posts.ratescond)

fixed.rates.cond <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.ratescond)
save(fixed.thresholds.cond, fixed.rates.cond, file="data/after_sampling/av.posts.cond.PPs.RData")


### Slightly different thing- don't average but remove PM related mechanisms.

##No proactive:
## Pick out the PM block ongoing task thresholds and replace them with the control
#thresholds.

thres2 <- c("B.A2C",        "B.B2C",        "B.C2C" ,      
            "B.D2C"  ,             "B.A2N"  ,      "B.B2N"  ,      "B.C2N"  ,     
            "B.D2N")

thres3  <- c("B.A3C"   ,     "B.B3C",        "B.C3C" ,      
             "B.D3C",        "B.A3N"   ,     "B.B3N" ,       "B.C3N" ,      
             "B.D3N" )

no_proactive <- pickps.h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples, 
                                          pickps_others=thres3,
                                          pickps_set=thres2, save.simulation=T)


parnames<- names(attr(attr(samples[[1]]$data, "model"), "p.vector"))
rates<- parnames[grep("mean_v.", parnames)]
rates<- rates[!(grepl("2", rates))]

nonPM_ongoingrates <- c(rates[grep("nn", rates)], rates[grep("cc", rates)])
PM_ongoingrates <-c(rates[grep("pn", rates)], rates[grep("pc", rates)])


#No reactive:
## Pick out the ongoing task accumulatoin rates on PM trials 
#and replace with ongoing task accumulation rates from non-PM trials (in PM blokcs)
no_reactive <- pickps.h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples, 
                                         pickps_others=PM_ongoingrates,
                                         pickps_set=nonPM_ongoingrates, save.simulation=T)

defaultparams <- c(nonPM_ongoingrates, thres2)
controlparams <- c(PM_ongoingrates, thres3)

no_control <- pickps.h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples, 
                                        pickps_others=controlparams,
                                        pickps_set=defaultparams, save.simulation=T)

save(no_proactive , 
     no_reactive,
     no_control,
file="data/after_sampling/PM_mechanism_piecetest_PPs.RData")


cat("Your sims are done, now calculating effects for each.")

full_effects.cond <- get.effects.dmc(E1PP, fun=cond.effects)
thresholdfix_effects.cond <- get.effects.dmc(fixed.thresholds.cond, fun=cond.effects)
ratesfix_effects.cond <- get.effects.dmc(fixed.rates.cond, fun=cond.effects)
save(full_effects.cond, thresholdfix_effects.cond,ratesfix_effects.cond, file="data/after_sampling/av.posts.cond.effects.RData")


full_effects.block <-get.effects.dmc(E1PP, fun=block.effects.E1A4)
control_thres_effects_all <-get.effects.dmc(no_proactive, fun=block.effects.E1A4)
no_reactive_effects <- get.effects.dmc(no_reactive, fun=block.effects.E1A4)
no_control_effects <- get.effects.dmc(no_control, fun=block.effects.E1A4)

save(full_effects.block , 
     control_thres_effects_all,
     no_reactive_effects,
     no_control_effects, file="data/after_sampling/PM_mechanism_piecetest_effects.RData")

full_effects.block <- get.effects.dmc(E1PP, fun=block.effects.E1A4)
thresholdfix_effects.block <- get.effects.dmc(fixed.thresholds.block, fun=block.effects.E1A4)
ratesfix_effects.block <- get.effects.dmc(fixed.rates.block, fun=block.effects.E1A4)

save(full_effects.block, thresholdfix_effects.block,ratesfix_effects.block,
     file="data/after_sampling/av.posts.block.effects.RData")


samples <- E1.block.B.V_cond.B.V.PMV.samples
names(attr(attr(samples[[1]]$data, "model"), "p.vector"))
load("data/after_sampling/E1PP.RData")

av.posts.thresblock <- unique(c(
  "^B.A*C"  ,      "^B.B*C"  ,      "^B.C*C"    ,    "^B.D*C"    ,   
 "^B.A*C"  ,      "^B.B*C"  ,      "^B.C*C"    ,    "^B.D*C"  ,      "^B.A*N" ,      
"^B.B*N"    ,    "^B.C*N"    ,    "^B.D*N"   ,     "^B.A*N"   ,     "^B.B*N"  ,     
 "^B.C*N"  ,      "^B.D*N"
  
))
av.posts.thresblock <- glob2rx(av.posts.thresblock)

#note this function converts your av.posts to regex using glob2rx
fixed.thresholds.block <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                                  av.posts.thresblock)



rates<- c("mean_v.ccA*C", "mean_v.nnA*C", "mean_v.ccB*C",
 "mean_v.nnB*C", "mean_v.ccC*C", "mean_v.nnC*C", "mean_v.ccD*C", "mean_v.nnD*C",
 "mean_v.ccA*C", "mean_v.nnA*C", "mean_v.ccA*C", "mean_v.nnA*C", "mean_v.ccB*C",
 "mean_v.nnB*C", "mean_v.ccB*C", "mean_v.nnB*C", "mean_v.ccC*C", "mean_v.nnC*C",
"mean_v.ccC*C", "mean_v.nnC*C", "mean_v.ccD*C", "mean_v.nnD*C", "mean_v.ccD*C",
"mean_v.nnD*C", "mean_v.ccA*N", "mean_v.nnA*N", "mean_v.ccB*N", "mean_v.nnB*N",
"mean_v.ccC*N", "mean_v.nnC*N", "mean_v.ccD*N", "mean_v.nnD*N", "mean_v.ccA*N",
 "mean_v.nnA*N", "mean_v.ccA*N", "mean_v.nnA*N", "mean_v.ccB*N", "mean_v.nnB*N",
"mean_v.ccB*N", "mean_v.nnB*N", "mean_v.ccC*N", "mean_v.nnC*N", "mean_v.ccC*N",
"mean_v.nnC*N", "mean_v.ccD*N", "mean_v.nnD*N", "mean_v.ccD*N", "mean_v.nnD*N")

av.posts.ratesblock<- unique(rates)

av.posts.ratesblock <- glob2rx(av.posts.ratesblock)
fixed.rates.block <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.ratesblock)


# get.accfits.E1.A4(fixed.thresholds.block)
# get.RTfits.E1.A4(fixed.thresholds.block)


full_effects <- get.effects.dmc(E1PP, fun=block.effects.E1A4)
thresholdfix_effects <- get.effects.dmc(fixed.thresholds.block, fun=block.effects.E1A4)
ratesfix_effects <- get.effects.dmc(fixed.rates.block, fun=block.effects.E1A4)

effects<-full_effects
finish.blockdf.E1_A4 <- function(effects) {
  effects$S <- NA
  effects$DV <- NA
  effects$S[grep ("Conflict", rownames(effects))] <- "Conflict"
  effects$S[grep ("Nonconflict", rownames(effects))] <- "Nonconflict"
  effects$DV <- "Accuracy"
  effects$DV[grep ("RT", rownames(effects))] <- "Correct RT"
  effects$DV[grep ("(FA)", rownames(effects))] <- "Error RT"
  effects$DV[rownames(effects)=="PM Accuracy"] <- "PM Accuracy"
  effects$S[rownames(effects)=="PM Accuracy"] <- "PM"
  effects$S[rownames(effects)=="PM RT"] <- "PM"
  effects$DV[rownames(effects)=="PM RT"] <- "RT"
  effects
}

full_noPM <- finish.blockdf.E1_A4(full_effects)[-(1:3),]; full_noPM$model <- "Full"
threshold_noPM <- finish.blockdf.E1_A4(thresholdfix_effects)[-(1:2),]; threshold_noPM$model <- 
  "Averaged Thresholds"
rates_noPM <- finish.blockdf.E1_A4(ratesfix_effects)[-(1:3),]; rates_noPM$model <- 
  "Averaged Accumulation Rates"

all_effects_predictives<- rbind(full_noPM, threshold_noPM, rates_noPM)
all_effects_predictives$model<- factor(all_effects_predictives$model
                                       , levels= c("Full", "Averaged Thresholds",
                                                   "Averaged Accumulation Rates"))




plot <- ggplot(all_effects_predictives, aes(S, mean)) 
plot<- plot+ facet_grid(DV  ~  model, scales = "free", space = "fixed") + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + 
  geom_point(aes(S, data), pch=21, size=4, colour="black") +xlab("Stimulus Type")+
  theme(text = element_text(size=18)) +ylab("Cost (PM - Control)")
plot

rates<- c("mean_v.cc*2C","mean_v.nn*2C","mean_v.cc*2C",
"mean_v.nn*2C","mean_v.cc*2C","mean_v.nn*2C","mean_v.cc*2C","mean_v.nn*2C",
"mean_v.cc*3C","mean_v.nn*3C","mean_v.pc*3C","mean_v.pn*3C","mean_v.cc*3C",
"mean_v.nn*3C","mean_v.pc*3C","mean_v.pn*3C","mean_v.cc*3C","mean_v.nn*3C",
 "mean_v.pc*3C","mean_v.pn*3C","mean_v.cc*3C","mean_v.nn*3C","mean_v.pc*3C",
 "mean_v.pn*3C","mean_v.cc*2N","mean_v.nn*2N","mean_v.cc*2N","mean_v.nn*2N",
 "mean_v.cc*2N","mean_v.nn*2N","mean_v.cc*2N","mean_v.nn*2N","mean_v.cc*3N",
 "mean_v.nn*3N","mean_v.pc*3N","mean_v.pn*3N","mean_v.cc*3N","mean_v.nn*3N",
 "mean_v.pc*3N","mean_v.pn*3N","mean_v.cc*3N","mean_v.nn*3N","mean_v.pc*3N",
 "mean_v.pn*3N","mean_v.cc*3N","mean_v.nn*3N","mean_v.pc*3N","mean_v.pn*3N",
 "mean_v.pp*3P","mean_v.pp*3P","mean_v.pp*3P","mean_v.pp*3P")


av.posts.ratescond<- unique(rates)
av.posts.ratescond <- glob2rx(av.posts.ratescond)



fixed.rates.cond <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.ratescond)

thresholds <- c("B.*2C",        "B.*2C",        "B.*2C" ,      
 "B.*2C"  ,      "B.*3C"   ,     "B.*3C",        "B.*3C" ,      
"B.*3C" ,       "B.*2N"  ,      "B.*2N"  ,      "B.*2N"  ,     
 "B.*2N",        "B.*3N"   ,     "B.*3N" ,       "B.*3N" ,      
 "B.*3N"  ,      "B.*3P"   ,     "B.*3P"   ,     "B.*3P" ,      
 "B.*3P") 

av.posts.thresholdscond <- unique(thresholds)
av.posts.thresholdscond <- glob2rx(av.posts.thresholdscond)

fixed.thresholds.cond <-
  avps.h.post.predict.dmc(samples, save.simulation=T, av.posts=
                            av.posts.thresholdscond )
 
 
full_effects_cond <- get.effects.dmc(E1PP, fun=cond.effects)
thresholdfix_effects_cond <- get.effects.dmc(fixed.thresholds.cond, fun=cond.effects)
ratesfix_effects_cond <- get.effects.dmc(fixed.rates.cond, fun=cond.effects)


effects <- full_effects_cond
finish.conddf.E1_A4 <- function(effects) {
  effects$S <- NA
  effects$DV <- NA
  effects$contrast <- NA
  effects$S[grep ("Conflict", rownames(effects))] <- "Conflict"
  effects$S[grep ("Nonconflict", rownames(effects))] <- "Nonconflict"
  effects$DV <- "Accuracy"
  effects$DV[grep ("RT", rownames(effects))] <- "Correct RT"
  effects$DV[grep ("(FA)", rownames(effects))] <- "Error RT"
  # effects$DV[grep ("PM Acc", rownames(effects))] <- "PM Accuracy"
  effects$S[grep ("PM", rownames(effects))] <- "PM"
  # effects$S[rownames(effects)=="PM Accuracy"]
  effects<- effects[grepl("Diff", rownames(effects)),]
  effects$contrast[grep ("A-B", rownames(effects))] <- "A-B"
  effects$contrast[grep ("B-C", rownames(effects))] <- "B-C"
  effects$contrast[grep ("C-D", rownames(effects))] <- "C-D"
  effects
}

full_noPM <- finish.conddf.E1_A4(full_effects_cond); full_noPM$model <- "Full"
threshold_noPM <- finish.conddf.E1_A4(thresholdfix_effects_cond); threshold_noPM$model <- 
  "Averaged Thresholds"
rates_noPM <- finish.conddf.E1_A4(ratesfix_effects_cond); rates_noPM$model <- 
  "Averaged Accumulation Rates"

all_effects_predictives<- rbind(full_noPM, threshold_noPM, rates_noPM)
all_effects_predictives$model<- factor(all_effects_predictives$model
                                       , levels= c("Full", "Averaged Thresholds",
                                                   "Averaged Accumulation Rates"))

# ongoing_predictives <- all_effects_predictives[!(grepl("PM", 
#                                   rownames(all_effects_predictives))),]
# 
# PM_predictives <- all_effects_predictives[(grepl("PM", 
#                                                   rownames(all_effects_predictives))),]
# 


plotC <- ggplot(ongoing_predictives[ongoing_predictives$S=="Conflict",], aes(contrast, mean)) 
plotC <- plotC+
  
  # facet_wrap(c("S", "DV", "model"), scales = 'free_x', 
  #                   nrow=9, ncol=3) +
  # 
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  
  
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + 
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +xlab("Contrast")+
  theme(text = element_text(size=18)) +ylab("Difference")

plotC

plotN <- ggplot(ongoing_predictives[ongoing_predictives$S=="Nonconflict",], aes(contrast, mean)) 

plotN <- plotN+
  
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  
  
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + 
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +xlab("Contrast")+
  theme(text = element_text(size=18)) +ylab("Difference")

plotN

plotP <- ggplot(all_effects_predictives[all_effects_predictives$S=="PM",], aes(contrast, mean)) 

plotP <- plotP+
  
  facet_grid(DV  ~ model, scales = "free", space = "fixed") +
  
  
  geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + 
  geom_point(aes(contrast, data), pch=21, size=4, colour="black") +xlab("Contrast")+
  theme(text = element_text(size=18)) +ylab("Difference")

plotP





pickps_set <- c("B.A2C",        "B.B2C",        "B.C2C" ,      
                "B.D2C"  ,             "B.A2N"  ,      "B.B2N"  ,      "B.C2N"  ,     
                "B.D2N")

pickps_other <- c("B.A3C"   ,     "B.B3C",        "B.C3C" ,      
                  "B.D3C",        "B.A3N"   ,     "B.B3N" ,       "B.C3N" ,      
                  "B.D3N" )

no_proactive <- pickps.h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples, 
                                               pickps_others=pickps_other,
                          pickps_set=pickps_set, save.simulation=T)

control_thres_effects_all <-get.effects.dmc(no_proactive, fun=block.effects.E1A4)
full_effects <- get.effects.dmc(E1PP, fun=block.effects.E1A4)

full_noPM_PMacc <- finish.blockdf.E1_A4(full_effects)[1:3,]; full_noPM_PMacc$model <- "Full"
no_proactive_PMacc <- finish.blockdf.E1_A4(control_thres_effects_all)[1:3,]; no_proactive_PMacc$model <- "proactive"
rbind(full_noPM_PMacc, no_proactive_PMacc)




c(
  
  "mean_v.ccA3C","mean_v.nnA3C","mean_v.pcA3C","mean_v.pnA3C","mean_v.ccB3C",
  "mean_v.nnB3C","mean_v.pcB3C","mean_v.pnB3C","mean_v.ccC3C","mean_v.nnC3C",
  "mean_v.pcC3C","mean_v.pnC3C","mean_v.ccD3C","mean_v.nnD3C","mean_v.pcD3C",
  "mean_v.pnD3C",
  
  "mean_v.ccA3N",
  "mean_v.nnA3N","mean_v.pcA3N","mean_v.pnA3N","mean_v.ccB3N","mean_v.nnB3N",
  "mean_v.pcB3N","mean_v.pnB3N","mean_v.ccC3N","mean_v.nnC3N","mean_v.pcC3N",
  "mean_v.pnC3N","mean_v.ccD3N","mean_v.nnD3N","mean_v.pcD3N","mean_v.pnD3N")


c(
  
  
)

parnames<- names(attr(attr(samples[[1]]$data, "model"), "p.vector"))
rates<- parnames[grep("mean_v.", parnames)]
rates<- rates[!(grepl("2", rates))]

nonPM_ongoingrates <- c(rates[grep("nn", rates)], rates[grep("cc", rates)])
PM_ongoingrates <-c(rates[grep("pn", rates)], rates[grep("pc", rates)])

no_reactive <- pickps.h.post.predict.dmc(E1.block.B.V_cond.B.V.PMV.samples, 
                                          pickps_others=PM_ongoingrates,
                                          pickps_set=nonPM_ongoingrates, save.simulation=T)

no_reactive_effects <- get.effects.dmc(no_reactive, fun=block.effects.E1A4)


full_noPM_PMacc <- finish.blockdf.E1_A4(full_effects)[1:3,]; full_noPM_PMacc$model <- "Full"
no_reactive_PMacc <- finish.blockdf.E1_A4(no_reactive_effects)[1:3,]; no_reactive_PMacc$model <- "reactive"
rbind(full_noPM_PMacc, no_proactive_PMacc, no_reactive_PMacc)




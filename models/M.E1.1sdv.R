# # # PMFA Model with 1 PM Drift Rate and 1 fixed sdv # # #
#
#     B ~ cond * block * R
#
#     mean_v ~ S * cond * block * R
#
#     1 PM False Alarm Rate
#     
#     1 sdv fixed
#
#     1 PM Drift
#

rm(list=ls()) 
setwd("~/Modelling")
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")
setwd("~/Modelling/x1/samples")
load("okdats.E1.RData")

# # # Create match maps # # #

# # # MAPMEANV # # # 
#
#     S * cond * block * R
#
#     1 PM False Alarm Rate
#
#     N Rates = 53

mapmeanv <- empty.map(list(S=c("cc","nn","pc","pn"),
                           cond=c("A","B","C","D"),
                           block=c("2","3"),
                           R=c("C","N","P")),
                      
                      levels=c("ccA2C","nnA2C",
                               "ccB2C","nnB2C",
                               "ccC2C","nnC2C",
                               "ccD2C","nnD2C",
                               
                               "ccA3C","nnA3C","pcA3C","pnA3C",
                               "ccB3C","nnB3C","pcB3C","pnB3C",
                               "ccC3C","nnC3C","pcC3C","pnC3C",
                               "ccD3C","nnD3C","pcD3C","pnD3C",
                               
                               "ccA2N","nnA2N",
                               "ccB2N","nnB2N",
                               "ccC2N","nnC2N",
                               "ccD2N","nnD2N",
                               
                               "ccA3N","nnA3N","pcA3N","pnA3N",
                               "ccB3N","nnB3N","pcB3N","pnB3N",
                               "ccC3N","nnC3N","pcC3N","pnC3N",
                               "ccD3N","nnD3N","pcD3N","pnD3N",
                               
                               "ppA3P",
                               "ppB3P",
                               "ppC3P",
                               "ppD3P",
                               
                               "PMFA",
                               
                               "FAKERATE"))

mapmeanv

mapmeanv[1:96] <- c("ccA2C","nnA2C","FAKERATE","FAKERATE",
                    "ccB2C","nnB2C","FAKERATE","FAKERATE",
                    "ccC2C","nnC2C","FAKERATE","FAKERATE",
                    "ccD2C","nnD2C","FAKERATE","FAKERATE",
                    
                    "ccA3C","nnA3C","pcA3C","pnA3C",
                    "ccB3C","nnB3C","pcB3C","pnB3C",
                    "ccC3C","nnC3C","pcC3C","pnC3C",
                    "ccD3C","nnD3C","pcD3C","pnD3C",
                    
                    "ccA2N","nnA2N","FAKERATE","FAKERATE",
                    "ccB2N","nnB2N","FAKERATE","FAKERATE",
                    "ccC2N","nnC2N","FAKERATE","FAKERATE",
                    "ccD2N","nnD2N","FAKERATE","FAKERATE",
                    
                    "ccA3N","nnA3N","pcA3N","pnA3N",
                    "ccB3N","nnB3N","pcB3N","pnB3N",
                    "ccC3N","nnC3N","pcC3N","pnC3N",
                    "ccD3N","nnD3N","pcD3N","pnD3N",
                    
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    
                    "PMFA","PMFA","ppA3P","ppA3P",
                    "PMFA","PMFA","ppB3P","ppB3P",
                    "PMFA","PMFA","ppC3P","ppC3P",
                    "PMFA","PMFA","ppD3P","ppD3P")

mapmeanv


# # # MAPB # # # 
#
#     cond * block * R
#
#     N Thresholds = 20

mapB <- empty.map(list(S=c("cc","nn","pc","pn"),
                       cond=c("A","B","C","D"),
                       block=c("2","3"),
                       R=c("C","N","P")),
                  
                  levels=c("A2C",
                           "B2C",
                           "C2C",
                           "D2C",
                           
                           "A3C",
                           "B3C",
                           "C3C",
                           "D3C",
                           
                           "A2N",
                           "B2N",
                           "C2N",
                           "D2N",
                           
                           "A3N",
                           "B3N",
                           "C3N",
                           "D3N",
                           
                           "A3P",
                           "B3P",
                           "C3P",
                           "D3P",
                           
                           "FAKEB"))

# mapB

mapB[1:96] <- c("A2C","A2C","FAKEB","FAKEB",
                "B2C","B2C","FAKEB","FAKEB",
                "C2C","C2C","FAKEB","FAKEB",
                "D2C","D2C","FAKEB","FAKEB",
                
                "A3C","A3C","A3C","A3C",
                "B3C","B3C","B3C","B3C",
                "C3C","C3C","C3C","C3C",
                "D3C","D3C","D3C","D3C",
                
                "A2N","A2N","FAKEB","FAKEB",
                "B2N","B2N","FAKEB","FAKEB",
                "C2N","C2N","FAKEB","FAKEB",
                "D2N","D2N","FAKEB","FAKEB",
                
                "A3N","A3N","A3N","A3N",
                "B3N","B3N","B3N","B3N",
                "C3N","C3N","C3N","C3N",
                "D3N","D3N","D3N","D3N",
                
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                
                "A3P","A3P","A3P","A3P",
                "B3P","B3P","B3P","B3P",
                "C3P","C3P","C3P","C3P",
                "D3P","D3P","D3P","D3P")

# mapB



# # # Build model # # # 

model <- model.dmc(p.map = list(A="1",
                                B=c("MAPB"),
                                t0="1",
                                mean_v=c("MAPMV"),
                                sd_v="1",
                                st0="1",
                                N="block"), 
                   
                   match.map = list(M=list(cc="C",nn="N",pc="P",pn="P"),
                                    
                                    MAPB=mapB,
                                    MAPMV=mapmeanv
                                    
                   ),
                   
                   factors = list(S=c("cc","nn","pc","pn"),
                                  cond=c("A","B","C","D"),
                                  block=c("2","3")),
                   
                   constants = c(N.2=2,N.3=3,sd_v=0.5,st0=0,
                                 
                                 B.FAKEB=Inf,
                                 mean_v.FAKERATE=1
                                 
                   ), 
                   
                   responses = c("C","N","P"),
                   
                   type = "normN")


# # # Create parameter vector # # # 

p.vector <- c(t0=0.3, A=3,
              
              B.A2C=2,  B.B2C=2,  B.C2C=2,  B.D2C=2,  
              B.A3C=2,  B.B3C=2,  B.C3C=2,  B.D3C=2,  
              
              B.A2N=2,  B.B2N=2,  B.C2N=2,  B.D2N=2,  
              B.A3N=2,  B.B3N=2,  B.C3N=2,  B.D3N=2,  
              
              B.A3P=2,  B.B3P=2,  B.C3P=2,  B.D3P=2,
              
              mean_v.ccA2C=1, mean_v.nnA2C=0,  
              mean_v.ccA2N=0, mean_v.nnA2N=1, 
              
              mean_v.ccB2C=1, mean_v.nnB2C=0,  
              mean_v.ccB2N=0, mean_v.nnB2N=1, 
              
              mean_v.ccC2C=1, mean_v.nnC2C=0,  
              mean_v.ccC2N=0, mean_v.nnC2N=1, 
              
              mean_v.ccD2C=1, mean_v.nnD2C=0,  
              mean_v.ccD2N=0, mean_v.nnD2N=1, 
              
              mean_v.ccA3C=1, mean_v.nnA3C=0, mean_v.pnA3C=0, mean_v.pcA3C=1, 
              mean_v.ccA3N=0, mean_v.nnA3N=1, mean_v.pnA3N=1, mean_v.pcA3N=0,
                                              mean_v.ppA3P=1, 
              
              mean_v.ccB3C=1, mean_v.nnB3C=0, mean_v.pnB3C=0, mean_v.pcB3C=1, 
              mean_v.ccB3N=0, mean_v.nnB3N=1, mean_v.pnB3N=1, mean_v.pcB3N=0,
                                              mean_v.ppB3P=1, 
              
              mean_v.ccC3C=1, mean_v.nnC3C=0, mean_v.pnC3C=0, mean_v.pcC3C=1, 
              mean_v.ccC3N=0, mean_v.nnC3N=1, mean_v.pnC3N=1, mean_v.pcC3N=0,
                                              mean_v.ppC3P=1, 
              
              mean_v.ccD3C=1, mean_v.nnD3C=0, mean_v.pnD3C=0, mean_v.pcD3C=1, 
              mean_v.ccD3N=0, mean_v.nnD3N=1, mean_v.pnD3N=1, mean_v.pcD3N=0,
                                              mean_v.ppD3P=1, 
              
              mean_v.PMFA=0
)


# # # Check parameter vector matches model # # # 

check.p.vector(p.vector, model)

# simulate.dmc(p.vector, model)


# # # Set priors # # #  

p.prior <- prior.p.dmc(
  dists = c("beta",rep("tnorm", length(p.vector)-1)),
  p1=c(t0=1, p.vector[-1]),                           
  p2=c(1, 1, 
       rep(1,20), 
       rep(2,53)), 
  lower=c(0.1, 0, rep(0,20), rep(NA,53)),
  upper=c(1, 10, rep(Inf, length(p.vector)-2))
)

E1.1sdv <- data.model.dmc(okdats, model)
save(E1.1sdv, file="E1.1sdv.data.model.RData")

# # # Sampling # # #
  
E1.1sdv.samples <- h.samples.dmc(nmc=120, p.prior, E1.1sdv, thin=20)

save(E1.1sdv.samples, file="E1.1sdv.samples.RData")


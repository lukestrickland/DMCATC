rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


load("okdats.E1.RData")

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
                               
                               "pcA3P","pnA3P",
                               "pcB3P","pnB3P",
                               "pcC3P","pnC3P",
                               "pcD3P","pnD3P",
                               
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
                    
                    "PMFA","PMFA","pcA3P","pnA3P",
                    "PMFA","PMFA","pcB3P","pnB3P",
                    "PMFA","PMFA","pcC3P","pnC3P",
                    "PMFA","PMFA","pcD3P","pnD3P")

mapmeanv


# # # MAPB # # # 
#
#     block * R
#
#     N Thresholds = 5

mapB <- empty.map(list(S=c("cc","nn","pc","pn"),
                       cond=c("A","B","C","D"),
                       block=c("2","3"),
                       R=c("C","N","P")),
                  
                  levels=c("2C","3C",
                           
                           "2N","3N",
                           
                           "3P",
                           
                           "FAKEB"))

# mapB

mapB[1:96] <- c("2C","2C","FAKEB","FAKEB",
                "2C","2C","FAKEB","FAKEB",
                "2C","2C","FAKEB","FAKEB",
                "2C","2C","FAKEB","FAKEB",
                
                "3C","3C","3C","3C",
                "3C","3C","3C","3C",
                "3C","3C","3C","3C",
                "3C","3C","3C","3C",
                
                "2N","2N","FAKEB","FAKEB",
                "2N","2N","FAKEB","FAKEB",
                "2N","2N","FAKEB","FAKEB",
                "2N","2N","FAKEB","FAKEB",
                
                "3N","3N","3N","3N",
                "3N","3N","3N","3N",
                "3N","3N","3N","3N",
                "3N","3N","3N","3N",
                
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                "FAKEB","FAKEB","FAKEB","FAKEB",
                
                "3P","3P","3P","3P",
                "3P","3P","3P","3P",
                "3P","3P","3P","3P",
                "3P","3P","3P","3P")

mapB


# # # MAPSDV # # # 
#
#     S * R
#
#     1 PM False Alarm Rate
#
#     N SDVs = 11

mapsdv <- empty.map(list(S=c("cc","nn","pc","pn"),
                         cond=c("A","B","C","D"),
                         block=c("2","3"),
                         R=c("C","N","P")),
                    
                    levels=c("ccC","nnC","pcC","pnC",
                             
                             "ccN","nnN","pcN","pnN",
                             
                             "pcP","pnP",
                             
                             "PMFAsdv",
                             
                             "FAKESDV"))

mapsdv

mapsdv[1:96] <- c("ccC","nnC","FAKESDV","FAKESDV",
                  "ccC","nnC","FAKESDV","FAKESDV",
                  "ccC","nnC","FAKESDV","FAKESDV",
                  "ccC","nnC","FAKESDV","FAKESDV",
                  
                  "ccC","nnC","pcC","pnC",
                  "ccC","nnC","pcC","pnC",
                  "ccC","nnC","pcC","pnC",
                  "ccC","nnC","pcC","pnC",
                  
                  "ccN","nnN","FAKESDV","FAKESDV",
                  "ccN","nnN","FAKESDV","FAKESDV",
                  "ccN","nnN","FAKESDV","FAKESDV",
                  "ccN","nnN","FAKESDV","FAKESDV",
                  
                  "ccN","nnN","pcN","pnN",
                  "ccN","nnN","pcN","pnN",
                  "ccN","nnN","pcN","pnN",
                  "ccN","nnN","pcN","pnN",
                  
                  "FAKESDV","FAKESDV","FAKESDV","FAKESDV",
                  "FAKESDV","FAKESDV","FAKESDV","FAKESDV",
                  "FAKESDV","FAKESDV","FAKESDV","FAKESDV",
                  "FAKESDV","FAKESDV","FAKESDV","FAKESDV",
                  
                  "PMFAsdv","PMFAsdv","pcP","pnP",
                  "PMFAsdv","PMFAsdv","pcP","pnP",
                  "PMFAsdv","PMFAsdv","pcP","pnP",
                  "PMFAsdv","PMFAsdv","pcP","pnP")

mapsdv


# # # Build model # # # 

model <- model.dmc(p.map = list(A="1",
                                B=c("MAPB"),
                                t0="1",
                                mean_v=c("MAPMV"),
                                sd_v=c("MAPSDV"),
                                st0="1",
                                N="block"), 
                   
                   match.map = list(M=list(cc="C",nn="N",pc="P",pn="P"),
                                    
                                    MAPB=mapB,
                                    MAPMV=mapmeanv,
                                    MAPSDV=mapsdv
                                    
                   ),
                   
                   factors = list(S=c("cc","nn","pc","pn"),
                                  cond=c("A","B","C","D"),
                                  block=c("2","3")),
                   
                   constants = c(N.2=2,N.3=3,sd_v.PMFAsdv=0.5,st0=0,
                                 
                                 B.FAKEB=Inf,
                                 mean_v.FAKERATE=1,
                                 sd_v.FAKESDV=1
                                 
                   ), 
                   
                   responses = c("C","N","P"),
                   
                   type = "normN")


# # # Create parameter vector # # # 

p.vector <- c(t0=0.3, A=3,
              sd_v.ccC=0.5,  sd_v.nnC=0.5,  sd_v.pnC=0.5,  sd_v.pcC=0.5,  
              sd_v.ccN=0.5,  sd_v.nnN=0.5,  sd_v.pnN=0.5,  sd_v.pcN=0.5,       
              sd_v.pnP=0.5,  sd_v.pcP=0.5,   
              
              B.2C=2,
              B.3C=2,
              
              B.2N=2,
              B.3N=2,
              
              B.3P=2,
              
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
                                              mean_v.pnA3P=1, mean_v.pcA3P=1,
              
              mean_v.ccB3C=1, mean_v.nnB3C=0, mean_v.pnB3C=0, mean_v.pcB3C=1, 
              mean_v.ccB3N=0, mean_v.nnB3N=1, mean_v.pnB3N=1, mean_v.pcB3N=0,
                                              mean_v.pnB3P=1, mean_v.pcB3P=1,
              
              mean_v.ccC3C=1, mean_v.nnC3C=0, mean_v.pnC3C=0, mean_v.pcC3C=1, 
              mean_v.ccC3N=0, mean_v.nnC3N=1, mean_v.pnC3N=1, mean_v.pcC3N=0,
                                              mean_v.pnC3P=1, mean_v.pcC3P=1,
              
              mean_v.ccD3C=1, mean_v.nnD3C=0, mean_v.pnD3C=0, mean_v.pcD3C=1, 
              mean_v.ccD3N=0, mean_v.nnD3N=1, mean_v.pnD3N=1, mean_v.pcD3N=0,
                                              mean_v.pnD3P=1, mean_v.pcD3P=1,
              
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
       rep(1,10), 
       rep(1,5), 
       rep(2,57)), 
  lower=c(0.1, 0, rep(0,10), rep(0,5), rep(NA,57)),
  upper=c(1, 10, rep(Inf, length(p.vector)-2))
)

E1.block.B.V_cond.V <- data.model.dmc(okdats, model)

E1.condVonly.samples <- h.samples.dmc(nmc=120, p.prior, E1.block.B.V_cond.V , thin=20)
save(E1.condVonly.samples, file="E1.condVonly.samples.RData")


rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


run.grid.dmc("E1.condVonly.samples",model.dir="LBA",model.file="lbaN_B.R",user="ljs392",n.add=40,wall.hours=30,GB=5)  


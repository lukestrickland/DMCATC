rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


load("okdats.E1.RData")

mapmeanv <- empty.map(list(S=c("cc","nn","pc","pn"),
                           cond=c("A","B","C","D"),
                           block=c("2","3"),
                           R=c("C","N","P")),
                      
                      levels=c("cc2C","nn2C",
                               
                               "cc3C","nn3C","pc3C","pn3C",
                               
                               "cc2N","nn2N",
                               
                               "cc3N","nn3N","pc3N","pn3N",
                               
                               "pc3P","pn3P",
                               
                               "PMFA",
                               
                               "FAKERATE"))

mapmeanv

mapmeanv[1:96] <- c("cc2C","nn2C","FAKERATE","FAKERATE",
                    "cc2C","nn2C","FAKERATE","FAKERATE",
                    "cc2C","nn2C","FAKERATE","FAKERATE",
                    "cc2C","nn2C","FAKERATE","FAKERATE",
                    
                    "cc3C","nn3C","pc3C","pn3C",
                    "cc3C","nn3C","pc3C","pn3C",
                    "cc3C","nn3C","pc3C","pn3C",
                    "cc3C","nn3C","pc3C","pn3C",
                    
                    "cc2N","nn2N","FAKERATE","FAKERATE",
                    "cc2N","nn2N","FAKERATE","FAKERATE",
                    "cc2N","nn2N","FAKERATE","FAKERATE",
                    "cc2N","nn2N","FAKERATE","FAKERATE",
                    
                    "cc3N","nn3N","pc3N","pn3N",
                    "cc3N","nn3N","pc3N","pn3N",
                    "cc3N","nn3N","pc3N","pn3N",
                    "cc3N","nn3N","pc3N","pn3N",
                    
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    
                    "PMFA","PMFA","pc3P","pn3P",
                    "PMFA","PMFA","pc3P","pn3P",
                    "PMFA","PMFA","pc3P","pn3P",
                    "PMFA","PMFA","pc3P","pn3P")

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
                   
                   constants = c(N.2=2,N.3=3,sd_v.nnC=0.5,st0=0,
                                 
                                 B.FAKEB=Inf,
                                 mean_v.FAKERATE=1,
                                 sd_v.FAKESDV=1
                                 
                   ), 
                   
                   responses = c("C","N","P"),
                   
                   type = "normN")


# # # Create parameter vector # # # 

p.vector <- c(t0=0.3, A=3,
              sd_v.ccC=0.5,                 sd_v.pnC=0.5,  sd_v.pcC=0.5,  
              sd_v.ccN=0.5,  sd_v.nnN=0.5,  sd_v.pnN=0.5,  sd_v.pcN=0.5,       
              sd_v.pnP=0.5,  sd_v.pcP=0.5,  sd_v.PMFAsdv=0.5,   
              
              B.A2C=2,  B.B2C=2,  B.C2C=2,  B.D2C=2,  
              B.A3C=2,  B.B3C=2,  B.C3C=2,  B.D3C=2,  
              
              B.A2N=2,  B.B2N=2,  B.C2N=2,  B.D2N=2,  
              B.A3N=2,  B.B3N=2,  B.C3N=2,  B.D3N=2,  
              
              B.A3P=2,  B.B3P=2,  B.C3P=2,  B.D3P=2,
              
              mean_v.cc2C=1, mean_v.nn2C=0,  
              mean_v.cc2N=0, mean_v.nn2N=1, 
              
              mean_v.cc3C=1, mean_v.nn3C=0, mean_v.pn3C=0, mean_v.pc3C=1, 
              mean_v.cc3N=0, mean_v.nn3N=1, mean_v.pn3N=1, mean_v.pc3N=0,
                                            mean_v.pn3P=1, mean_v.pc3P=1,
              
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
       rep(1,20), 
       rep(2,15)), 
  lower=c(0.1, 0, rep(0,10), rep(0,20), rep(NA,15)),
  upper=c(1, 10, rep(Inf, length(p.vector)-2))
)

E1.block.B.V_cond.B <- data.model.dmc(okdats, model)

E1.condBonly.samples <- h.samples.dmc(nmc=120, p.prior, E1.block.B.V_cond.B , thin=20)

save(E1.condBonly.samples, file="E1.condBonly.samples.RData")


rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


run.grid.dmc("E1.condBonly.samples",model.dir="LBA",model.file="lbaN_B.R",user="ljs392",n.add=40,wall.hours=30,GB=5)  


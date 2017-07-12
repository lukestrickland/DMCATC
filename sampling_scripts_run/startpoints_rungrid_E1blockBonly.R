rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


load("okdats.E1.RData")


#     S * cond * R
#
#     1 PM False Alarm Rate
#
#     N Rates = 41

mapmeanv <- empty.map(list(S=c("cc","nn","pc","pn"),
                           cond=c("A","B","C","D"),
                           block=c("2","3"),
                           R=c("C","N","P")),
                      
                      levels=c("ccAC","nnAC","pcAC","pnAC",
                               "ccBC","nnBC","pcBC","pnBC",
                               "ccCC","nnCC","pcCC","pnCC",
                               "ccDC","nnDC","pcDC","pnDC",
                               
                               "ccAN","nnAN","pcAN","pnAN",
                               "ccBN","nnBN","pcBN","pnBN",
                               "ccCN","nnCN","pcCN","pnCN",
                               "ccDN","nnDN","pcDN","pnDN",
                               
                               "pcAP","pnAP",
                               "pcBP","pnBP",
                               "pcCP","pnCP",
                               "pcDP","pnDP",
                               
                               "PMFA",
                               
                               "FAKERATE"))

mapmeanv

mapmeanv[1:96] <- c("ccAC","nnAC","FAKERATE","FAKERATE",
                    "ccBC","nnBC","FAKERATE","FAKERATE",
                    "ccCC","nnCC","FAKERATE","FAKERATE",
                    "ccDC","nnDC","FAKERATE","FAKERATE",
                    
                    "ccAC","nnAC","pcAC","pnAC",
                    "ccBC","nnBC","pcBC","pnBC",
                    "ccCC","nnCC","pcCC","pnCC",
                    "ccDC","nnDC","pcDC","pnDC",
                    
                    "ccAN","nnAN","FAKERATE","FAKERATE",
                    "ccBN","nnBN","FAKERATE","FAKERATE",
                    "ccCN","nnCN","FAKERATE","FAKERATE",
                    "ccDN","nnDN","FAKERATE","FAKERATE",
                    
                    "ccAN","nnAN","pcAN","pnAN",
                    "ccBN","nnBN","pcBN","pnBN",
                    "ccCN","nnCN","pcCN","pnCN",
                    "ccDN","nnDN","pcDN","pnDN",
                    
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    "FAKERATE","FAKERATE","FAKERATE","FAKERATE",
                    
                    "PMFA","PMFA","pcAP","pnAP",
                    "PMFA","PMFA","pcBP","pnBP",
                    "PMFA","PMFA","pcCP","pnCP",
                    "PMFA","PMFA","pcDP","pnDP")

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
              
              B.A2C=2,  B.B2C=2,  B.C2C=2,  B.D2C=2,  
              B.A3C=2,  B.B3C=2,  B.C3C=2,  B.D3C=2,  
              
              B.A2N=2,  B.B2N=2,  B.C2N=2,  B.D2N=2,  
              B.A3N=2,  B.B3N=2,  B.C3N=2,  B.D3N=2,  
              
              B.A3P=2,  B.B3P=2,  B.C3P=2,  B.D3P=2,
              
              mean_v.ccAC=1, mean_v.nnAC=0, mean_v.pnAC=0, mean_v.pcAC=1, 
              mean_v.ccAN=0, mean_v.nnAN=1, mean_v.pnAN=1, mean_v.pcAN=0,
                                            mean_v.pnAP=1, mean_v.pcAP=1,
              
              mean_v.ccBC=1, mean_v.nnBC=0, mean_v.pnBC=0, mean_v.pcBC=1, 
              mean_v.ccBN=0, mean_v.nnBN=1, mean_v.pnBN=1, mean_v.pcBN=0,
                                            mean_v.pnBP=1, mean_v.pcBP=1,
              
              mean_v.ccCC=1, mean_v.nnCC=0, mean_v.pnCC=0, mean_v.pcCC=1, 
              mean_v.ccCN=0, mean_v.nnCN=1, mean_v.pnCN=1, mean_v.pcCN=0,
                                            mean_v.pnCP=1, mean_v.pcCP=1,
              
              mean_v.ccDC=1, mean_v.nnDC=0, mean_v.pnDC=0, mean_v.pcDC=1, 
              mean_v.ccDN=0, mean_v.nnDN=1, mean_v.pnDN=1, mean_v.pcDN=0,
                                            mean_v.pnDP=1, mean_v.pcDP=1,
              
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
       rep(2,41)), 
  lower=c(0.1, 0, rep(0,10), rep(0,20), rep(NA,41)),
  upper=c(1, 10, rep(Inf, length(p.vector)-2))
)

E1.block.B_cond.B.V <- data.model.dmc(okdats, model)


E1.blockBonly.samples <- h.samples.dmc(nmc=120, p.prior, E1.block.B_cond.B.V  , thin=20)
save(E1.blockBonly.samples, file="E1.blockBonly.samples.RData")


rm(list=ls()) 
setwd("~/DMCATC_grid")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")


run.grid.dmc("E1.blockBonly.samples",model.dir="LBA",model.file="lbaN_B.R",user="ljs392",n.add=40,wall.hours=30,GB=5)  


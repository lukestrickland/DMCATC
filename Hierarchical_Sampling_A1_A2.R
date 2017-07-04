rm(list=ls()) 
setwd("D:/Software/DMC_ATCPMDC")
source ("dmc/dmc.R")
load_model ("LBA","lbaN_B.R")
load("data/samples/A1.block.B.V_cond.B.V.PMV.samples.RData")
cores=64

load("data/exp_data/A1.block.B.V_cond.B.V.PMV.data.model.RData")

dm <- A1.block.B.V_cond.B.V.PM

A1 <- A1.block.B.V_cond.B.V.PMV.samples

p.prior <- A1[[1]]$p.prior
p1 <- get.p.vector(A1[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
  dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)

hstart <- make.hstart(A1)
theta1 <- make.theta1(A1)
hA1 <- 
  h.samples.dmc(nmc=120,p.prior,dm,pp.prior,
  hstart.prior=hstart,theta1=theta1,thin=20)

# STAGE 2
hA1  <- h.run.unstuck.dmc(hA1, p.migrate = .05, cores = cores)
save(A1,hA1,file="rsA1.RData")
hA1 <- h.run.converge.dmc(h.samples.dmc(nmc=120, samples=hA1), 
  thorough=TRUE,cores=cores,finalrun=TRUE,finalI=120)
save(A1,hA1,file="rsA1.RData")

# POST PROCESS

pdf("rshA1_chains.pdf",height=6,width = 8)
plot.dmc(hA1,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hA1,hyper=TRUE,layout=c(2,5))
plot.dmc(hA1,hyper=TRUE,layout=c(2,5),p.prior=pp.prior)
dev.off()

ppA1 <- h.post.predict.dmc(hA1,cores=32)
save(A1,hA1,ppA1,file="rsA1.RData")
pdf("rshA1_fit.pdf",height=6, width = 8)
plot.pp.dmc(ppA1,model.legend = FALSE,layout=c(2,3))
dev.off()

A1.pll <- group_trial_log_likes(hA1,thin_pointwise = 10,cores=16)
A1.waic <- waic.dmc(A1.pll,digits=2,save=TRUE)
save(A1,hA1,ppA1,A1.waic,file="rsA1.RData")
h.IC.dmc(hA1,DIC=TRUE)



gelman.diag.dmc(hA1,hyper=TRUE)
gelman.diag.dmc(hA1)

tmp <- summary.dmc(hA1,hyper=TRUE)$quantiles
round(tmp[,c(1,3,5)],3)


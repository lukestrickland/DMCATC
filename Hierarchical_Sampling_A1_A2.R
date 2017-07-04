rm(list=ls()) 
source ("dmc/dmc.R")
load_model ("LBA","lba_BcsuPPT.R")
load("rscsuPPT.RData")
cores=30

# STAGE 1
csuPPT <- h.RUN.dmc(csuPPT,cores=32,verbose=TRUE)
save(dm,csuPPT,file="rscsuPPT.RData")

# Make hierarchical
p.prior <- csuPPT[[1]]$p.prior
p1 <- get.p.vector(csuPPT[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
  dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)

hstart <- make.hstart(csuPPT)
theta1 <- make.theta1(csuPPT)
hcsuPPT <- h.samples.dmc(nmc=100,p.prior,dm,pp.prior,
  hstart.prior=hstart,theta1=theta1,thin=10)

# STAGE 2
hcsuPPT  <- h.run.unstuck.dmc(hcsuPPT, p.migrate = .05, cores = cores)
save(csuPPT,hcsuPPT,file="rscsuPPT.RData")
hcsuPPT <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hcsuPPT), 
  thorough=TRUE,nmc=50,cores=cores,finalrun=TRUE,finalI=500)
save(csuPPT,hcsuPPT,file="rscsuPPT.RData")
# hcsuPPT <- h.run.dmc(h.samples.dmc(samples=hcsuPPT,nmc=2000),cores=cores,report=100)
# save(csuPPT,hcsuPPT,file="rscsuPPT.RData")
# hcsuPPT <- h.run.dmc(h.samples.dmc(samples=hcsuPPT,nmc=2000),cores=cores,report=100)
# save(csuPPT,hcsuPPT,file="rscsuPPT.RData")
# hcsuPPT <- h.run.dmc(h.samples.dmc(samples=hcsuPPT,nmc=500),cores=cores,report=100)
# save(csuPPT,hcsuPPT,file="rscsuPPT.RData")

# POST PROCESS

pdf("rshcsuPPT_chains.pdf",height=6,width = 8)
plot.dmc(hcsuPPT,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hcsuPPT,hyper=TRUE,layout=c(2,5))
plot.dmc(hcsuPPT,hyper=TRUE,layout=c(2,5),p.prior=pp.prior)
dev.off()

ppcsuPPT <- h.post.predict.dmc(hcsuPPT,cores=32)
save(csuPPT,hcsuPPT,ppcsuPPT,file="rscsuPPT.RData")
pdf("rshcsuPPT_fit.pdf",height=6, width = 8)
plot.pp.dmc(ppcsuPPT,model.legend = FALSE,layout=c(2,3))
dev.off()

csuPPT.pll <- group_trial_log_likes(hcsuPPT,thin_pointwise = 10,cores=16)
csuPPT.waic <- waic.dmc(csuPPT.pll,digits=2,save=TRUE)
save(csuPPT,hcsuPPT,ppcsuPPT,csuPPT.waic,file="rscsuPPT.RData")
h.IC.dmc(hcsuPPT,DIC=TRUE)



gelman.diag.dmc(hcsuPPT,hyper=TRUE)
gelman.diag.dmc(hcsuPPT)

tmp <- summary.dmc(hcsuPPT,hyper=TRUE)$quantiles
round(tmp[,c(1,3,5)],3)


# pdf("rscsuPPT_fits.pdf",height=10,width = 12)
# for (i in 1:length(hcsuPPT))
#   plot.pp.dmc(ppcsuPPT[[i]],"cdf",model.legend = FALSE,layout=c(5,6))
# dev.off()

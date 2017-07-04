

samples=samples[[1]];n.post=100;probs=c(1:99)/100;random=TRUE
bw="nrd0";report=10;save.simulation=FALSE;factors=NA
save.simulation.as.attribute=FALSE;ignore.R2=FALSE
gglist=FALSE; probs.gglist=c(0.1, 0.5, 0.9);CI.gglist=c(0.025, 0.975)


post.predict.dmc.MATCHORDER <- function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
         bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
         save.simulation.as.attribute=FALSE,ignore.R2=FALSE,
         gglist=FALSE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975))
  # make list of posterior preditive density, quantiles and response p(robability)
  # NB: quantiles only calcualted for 2 or more RTs
{ 
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  cvs <- samples$data[,attr(model,"cvs")]
  attr(cvs,"row.facs") <- apply(apply(
    samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  if ( ignore.R2 & any(names(samples$data)=="R2") )
    samples$data <- samples$data[,names(samples$data)[names(samples$data)!="R2"]]
  if (!is.null(factors) ) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  resp <- names(attr(model,"responses"))
  ##LUKE: plug in data size from okdats.noNR
  data <- attr(samples, "NRdata")
  ns <- table(data[,facs])
  # ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(data)[2]))
  names(sim) <- names(data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity)) 
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,n=ns,SSD=SSD,cvs=NULL)

    if ( (i %% report) == 0) cat(".")
    ###Luke plug in order matching
    
    callargs.data <- list()
    callargs.sim <- list()
    for (p in 1:length(facs)) {callargs.data[[p]] <- data[,facs[p]]
    callargs.sim[[p]] <- tmp[,facs[p]]
    }
    
    data.ind <- factor(do.call(paste, callargs.data))
    sim.ind <- factor(do.call(paste, callargs.sim))
    swappedsim <- data
    for(q in levels(data.ind)){swappedsim$R[data.ind==q] <- tmp$R[sim.ind==q]
    swappedsim$RT[data.ind==q] <- tmp$RT[sim.ind==q]
    }
    tmp <- swappedsim
    
    if (ignore.R2) tmp <- tmp[,names(tmp)[names(tmp)!="R2"]]
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    
    ########
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="") 
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="") 
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  reps <- rep(1:n.post,each=dim(data)[1])
  if ( save.simulation ) {
    sim <- cbind(reps,sim)
    attr(sim,"data") <- data
    # attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[reps==i,]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(out,"dpqs") <- dpqs
    if (save.simulation.as.attribute) 
      attr(out,"sim") <- cbind(reps,sim)
    if (gglist) attr(out, "gglist") <- 
      get.fitgglist.dmc(cbind(reps,sim),samples$data,factors=factors, noR=FALSE, 
                        quantiles.to.get= probs.gglist, CI = CI.gglist)
    out
  }
}

###non-response RATES
psim.match.data.order <- function (samples,n.post=100,report=10) {
  data <- attr(samples, "NRdata")
  model<- attr(samples$data, "model")
  facs <- names(attr(model,"factors")); nfacs <- length(facs)
  ns <- table(data[,facs])
#get the average p.vector
  
  n.rep <- sum(ns)     ## total number of data points
  n.par <- dim(samples$theta)[2]

  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:nrow(thetas)), n.post, replace=F),]
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  report <- 1
  sim <- data.frame()
  for (i in 1:n.post) {
    currentsim <- simulate.dmc(posts[i,],model,ns)
    if ( (i %% report) == 0) message(".", appendLF = FALSE)
    
    
    callargs.data <- list()
    callargs.sim <- list()
    for (p in 1:length(facs)) {callargs.data[[p]] <- data[,facs[p]]
      callargs.sim[[p]] <- currentsim[,facs[p]]
    }

    data.ind <- factor(do.call(paste, callargs.data))
     sim.ind <- factor(do.call(paste, callargs.sim))
    swappedsim <- data
    for(q in levels(data.ind)){swappedsim$R[data.ind==q] <- currentsim$R[sim.ind==q]
        swappedsim$RT[data.ind==q] <- currentsim$RT[sim.ind==q]
       }
      currentsim <- swappedsim
      currentsim$rep <- i
      if (i==1) sim <- currentsim else sim <- rbind(sim, currentsim)
  }
  
  
  
    sim
}


h.psim.match.data.order <- function(hsamples, n.post=100, cores=1) {
  
  if (cores==1) {lapply (hsamples,
                         psim.match.data.order,
                  n.post=n.post)} else {mclapply (hsamples,
                         psim.match.data.order,
                  n.post=n.post, mc.cores = getOption("mc.cores", cores))}

}


get.trials.missed.E1_A4 <- function (sim) {

  Amissed <- c(); Bmissed <- c(); Cmissed <- c(); Dmissed <- c()
  for (i in 1:length(unique(sim$rep))) {
    sim1 <- sim[sim$rep==i,]
    TrialRTsA <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="A"] + 
      sim1$RT[sim1$trial.pos==2 & sim1$cond=="A"]
    TrialRTsB <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="B"] +
      sim1$RT[sim1$trial.pos==2 & sim1$cond=="B"]
    TrialRTsC <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="C"] + 
      sim1$RT[sim1$trial.pos==2 & sim1$cond=="C"]+ 
      sim1$RT[sim1$trial.pos==3 & sim1$cond=="C"]+ 
      sim1$RT[sim1$trial.pos==4 & sim1$cond=="C"]+ 
      sim1$RT[sim1$trial.pos==5 & sim1$cond=="C"]
    TrialRTsD <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="D"] + 
      sim1$RT[sim1$trial.pos==2 & sim1$cond=="D"]+ 
      sim1$RT[sim1$trial.pos==3 & sim1$cond=="D"]+ 
      sim1$RT[sim1$trial.pos==4 & sim1$cond=="D"]+
      sim1$RT[sim1$trial.pos==5 & sim1$cond=="D"]
    Amissed[i] <- sum(TrialRTsA>12) / length(TrialRTsA)
    Bmissed [i] <- sum(TrialRTsB>8)/ length(TrialRTsB)
    Cmissed [i] <- sum(TrialRTsC>20)/ length(TrialRTsC)
    Dmissed [i] <- sum(TrialRTsD>10) / length(TrialRTsD)

  }

  missed <- cbind(Amissed, Bmissed, Cmissed, Dmissed)
  colnames (missed) <- c("A", "B", "C", "D")
  missed
}


get.missavs <- function(ordered_sims, cores=cores, hsamples, fun=NA) {
  
  # if (is.na(fun)) stop("You need to specify a NR function")
  
## First apply the non-response function to get the simulated missrates
#for all participants. 
missrates <- mclapply(FUN=fun, ordered_sims,
                    mc.cores = getOption("mc.cores", cores))

#Then use this to get the sim df for the entire experiment
#we are averaging over participants here and then averaging the result
ggplot.misses.averaged <- function (missrates, lower=0.025, upper= 0.975) {
  stats <- lapply (missrates, function (y) {
      df <- data.frame(t(apply(y, 2, FUN= function(x) c(mean(x), 
                          quantile(x, probs= lower), quantile(x, probs= upper)))))
     df$condition <- rownames(df); colnames(df)[1:3] <- c("mean", "upper",
                                                          "lower")
     df})
  
  
  simdf<- do.call("rbind", stats)
  averaged <- data.frame (cbind(tapply(simdf$mean, list(simdf$condition), mean), 
                                tapply(simdf$upper, list(simdf$condition), mean), 
                                tapply(simdf$lower, list(simdf$condition), mean)))
  averaged$cond <- rownames(averaged)
  colnames(averaged)[1:3] <-  c("mean", "upper", "lower") 
  averaged
}


get.data.misses<- function(data) {
  datatrialmiss <- ddply(data, .(cond, trial), summarize, M=any(R=="M"))
  length(datatrialmiss$M[datatrialmiss$M]) / length(datatrialmiss$M)
  tapply (datatrialmiss$M, list(datatrialmiss$cond), mean)
}


names(missrates) <- names( hsamples)
sim_miss <- ggplot.misses.averaged(missrates)

for (i in names(hsamples)){
   
  data <- attr(hsamples[[i]], "NRdata")
  data$trial <- NA
    g=1
   for (t in 1:length(data$RT))  {
    if (t==1) data$trial[t] <- 1 else if (data$trial.pos[t] ==  data$trial.pos[t-1] +1) data$trial[t] <- g else {
       g <- g+1 
        data$trial[t] <- g} 
   }
  
    if (i==names(hsamples)[1]) outmat <- get.data.misses(data) else outmat <- cbind(outmat, get.data.misses(data))
    }
data_miss <- rowMeans(outmat)
df <- cbind(sim_miss, data_miss)
colnames(df)[5] <- "data"
df

}

#####Non-response RTs




psim.with.cut.RTs.E1.A4 <- function (samples,n.post=100,report=10) {
get.trialsums <- function(test) {
sumtest<-c()
for (i in 1:nrow(test)) {
  if (i==1) sumtest[i] <- test$RT[i] else if (test$trial[i]==test$trial[i-1]) sumtest [i] <- test$RT[i] + sumtest[i-1] else sumtest[i] <- test$RT[i]
    
}
cbind(test, sumtest)
}


cut.NRs <- function (test2){
test3 <- test2[!(test2$cond=="A" & test2$sumtest>12) & !(test2$cond=="B" & test2$sumtest>8)& !(test2$cond=="C" & test2$sumtest>20)& !(test2$cond=="D" & test2$sumtest>10),]
}

  data <- attr(samples, "NRdata")
  model<- attr(samples$data, "model")
  facs <- names(attr(model,"factors")); nfacs <- length(facs)
  ns <- table(data[,facs])
#get the average p.vector
  
  n.rep <- sum(ns)     ## total number of data points
  n.par <- dim(samples$theta)[2]

  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:nrow(thetas)), n.post, replace=F),]
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  report <- 1
  sim <- data.frame()
  for (i in 1:n.post) {
    currentsim <- simulate.dmc(posts[i,],model,ns)
    if ( (i %% report) == 0) message(".", appendLF = FALSE)
    
    
    callargs.data <- list()
    callargs.sim <- list()
    for (p in 1:length(facs)) {callargs.data[[p]] <- data[,facs[p]]
      callargs.sim[[p]] <- currentsim[,facs[p]]
    }

    data.ind <- factor(do.call(paste, callargs.data))
     sim.ind <- factor(do.call(paste, callargs.sim))
    swappedsim <- data
    for(q in levels(data.ind)){swappedsim$R[data.ind==q] <- currentsim$R[sim.ind==q]
        swappedsim$RT[data.ind==q] <- currentsim$RT[sim.ind==q]
    }
      swappedsim <- get.trialsums(swappedsim)
      swappedsim <- cut.NRs(swappedsim)
      currentsim <- swappedsim
      currentsim$rep <- i
      if (i==1) sim <- currentsim else sim <- rbind(sim, currentsim)
  }
   attr(sim, "data") <-samples$data
   sim
}

#lapply(hsamples,psim.with.cut.RTs.E1.A4, n.post=5)



    
    
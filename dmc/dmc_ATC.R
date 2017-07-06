
load_model ("LBA","lbaN_B.R")
require("gridExtra")
require("lme4")
require("plyr")
require("dplyr")
require("data.table")
theme_set(theme_simple())

# 
# samples=samples[[1]];n.post=100;probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=FALSE;factors=NA
# save.simulation.as.attribute=FALSE;ignore.R2=FALSE
# gglist=FALSE; probs.gglist=c(0.1, 0.5, 0.9);CI.gglist=c(0.025, 0.975)


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

h.post.predict.dmc.MATCHORDER <- function(hsamples) {
  lapply(hsamples, post.predict.dmc.MATCHORDER, save.simulation=TRUE)
}

get.trials.missed.E1_A4 <- function (sim) {

  Amissed <- c(); Bmissed <- c(); Cmissed <- c(); Dmissed <- c()
  for (i in 1:length(unique(sim$reps))) {
    sim1 <- sim[sim$reps==i,]
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


# 
# get.data.misses<- function(data) {
#   datatrialmiss <- ddply(data, .(cond, trial), summarize, M=any(R=="M"))
#   length(datatrialmiss$M[datatrialmiss$M]) / length(datatrialmiss$M)
#   tapply (datatrialmiss$M, list(datatrialmiss$cond), mean)
# }
# 
# get.groupNRs.ATCDMC <- function(sim, data, fun=NA, lower=.025, upper=.975) {
#   
#   
#   df <- data.frame(t(apply(fun(sim), 2, FUN= function(x) c(mean(x), 
#                                                            quantile(x, probs= lower), 
#                                                            quantile(x, probs= upper)))))
#   data$trial <- NA
#   data$trial.pos <-as.numeric(data$trial.pos)
#   g=1
#   for (t in 1:length(data$RT))  {
#     if (t==1) data$trial[t] <- 1 else if (data$trial.pos[t] ==  data$trial.pos[t-1] +1) data$trial[t] <- g 
#     else {
#       g <- g+1 
#       data$trial[t] <- g} 
#   }
#   cbind(df,get.data.misses(data))
#   condition <- rownames(df)
#   df<-cbind(condition, cbind(df,get.data.misses(data)))
#   colnames(df)[1:5] <- c("cond", "mean", "lower", "upper", "data")
#   df
# }

add.trial.cumsum.data <- function(df) {
  df$trial <- NA
  df$trial.pos <-as.numeric(df$trial.pos)
  
  df$trial <- NA
  df$trial.pos <-as.numeric(df$trial.pos)
  g=1
  for (t in 1:length(df$RT))  {
    if (t==1) df$trial[t] <- 1 else if (df$trial.pos[t] ==  df$trial.pos[t-1] +1) df$trial[t] <- g 
    else {
      g <- g+1 
      df$trial[t] <- g} 
  }
  df<-cbind(df,unlist(by(df$RT,df$trial,cumsum)))
  names(df)[length(df)]<- "cumsum" 
  df
}



add.trial.cumsum.sim <- function (sim, data) {
  if (!all.equal(rep(data$trial.pos, 100), sim$trial.pos)){
    stop("Data and Sim Trial Positions do not match")
  }
  sim$trial <- rep(data$trial, max(sim$reps))
  
  for(i in 1:max(sim$reps)) {
    sim$trial[sim$reps==i] <- data$trial + (i-1)*max(data$trial)
  }
  sim<-data.table(sim) 
  sim<-sim[,list(
    reps,      cond ,     block,     S   ,   
    R   ,      RT ,       trial.pos,

    
    cumsum=cumsum(RT)),list(trial)] 
  data.frame(sim)
}


# 
# 
# get.grouptrials.missed.E1_A4 <- function (sim) {
#   
#   Amissed <- c(); Bmissed <- c(); Cmissed <- c(); Dmissed <- c()
#   for (i in 1:length(unique(sim$reps))) {
#     sim1 <- sim[sim$reps==i,]
#     TrialRTsA <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="A"] + 
#       sim1$RT[sim1$trial.pos==2 & sim1$cond=="A"]
#     TrialRTsB <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="B"] +
#       sim1$RT[sim1$trial.pos==2 & sim1$cond=="B"]
#     TrialRTsC <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="C"] + 
#       sim1$RT[sim1$trial.pos==2 & sim1$cond=="C"]+ 
#       sim1$RT[sim1$trial.pos==3 & sim1$cond=="C"]+ 
#       sim1$RT[sim1$trial.pos==4 & sim1$cond=="C"]+ 
#       sim1$RT[sim1$trial.pos==5 & sim1$cond=="C"]
#     TrialRTsD <- sim1$RT[sim1$trial.pos==1 & sim1$cond=="D"] + 
#       sim1$RT[sim1$trial.pos==2 & sim1$cond=="D"]+ 
#       sim1$RT[sim1$trial.pos==3 & sim1$cond=="D"]+ 
#       sim1$RT[sim1$trial.pos==4 & sim1$cond=="D"]+
#       sim1$RT[sim1$trial.pos==5 & sim1$cond=="D"]
#     Amissed[i] <- sum(TrialRTsA>12) / length(TrialRTsA)
#     Bmissed [i] <- sum(TrialRTsB>8)/ length(TrialRTsB)
#     Cmissed [i] <- sum(TrialRTsC>20)/ length(TrialRTsC)
#     Dmissed [i] <- sum(TrialRTsD>10) / length(TrialRTsD)
#     
#   }
#   
#   missed <- cbind(Amissed, Bmissed, Cmissed, Dmissed)
#   colnames (missed) <- c("A", "B", "C", "D")
#   missed
# }
# 
# 

get.trials.missed.E1_A4 <- function(sim) {
  Amissed <- c(); Bmissed <- c(); Cmissed <- c(); Dmissed <- c()
  for (i in 1:max(sim$reps)) {
    simi <- sim[sim$reps==i,]
    Amissed[i]= length(simi$RT[simi$cond=="A" & simi$cumsum >12])/ length(simi$RT[simi$cond=="A"])
    Bmissed [i]= length(simi$RT[simi$cond=="B" & simi$cumsum >8])/ length(simi$RT[simi$cond=="B"])
    Cmissed [i]= length(simi$RT[simi$cond=="C" & simi$cumsum >20])/ length(simi$RT[simi$cond=="C"])
    Dmissed [i]= length(simi$RT[simi$cond=="D" & simi$cumsum >10])/ length(simi$RT[simi$cond=="D"])
  }
  missed <- cbind(Amissed, Bmissed, Cmissed, Dmissed)
  colnames (missed) <- c("A", "B", "C", "D")
  missed
}


get.NRs.ATCDMC <- function(sim, data, miss_fun=NA, lower=.025, upper=.975) {
  missrates <- miss_fun(sim)
  df <- data.frame(t(apply(missrates, 2, FUN= function(x) c(mean(x), 
                                                           quantile(x, probs= lower), 
                                                           quantile(x, probs= upper)))))

  data_misses<- tapply(data$R=="M", data$cond, mean)
  df<- cbind(df,data_misses)
  condition <- rownames(df)
  df<-cbind(condition, df)
  colnames(df)[1:5] <- c("cond", "mean", "lower", "upper", "data")
  df
}


    
    
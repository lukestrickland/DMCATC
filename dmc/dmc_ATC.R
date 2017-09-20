
load_model ("LBA","lbaN_B.R")
require("gridExtra")
require("lme4")
require("plyr")
require("dplyr")
require("data.table")
theme_set(theme_simple())

# A few functions for posterior predicctive p values and z scores.
minp <- function (effect) min(ecdf(effect)(0), 1-ecdf(effect)(0))

zandp <- function(samples, fun){
  effect<- group.inference.dist(samples, fun)
  Z <- mean(effect)/sd(effect) }

  Z.p.acrossexp <- function(samples1,samples2, fun){
    effect1<- group.inference.dist(samples1, fun)
    effect2 <- group.inference.dist(samples2, fun)
    effect<- effect1 - effect2
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), " (", round(p,3), ")", sep="")
  }

##accepts a function and does it to the thetas for each subject then averages after
group.inference.dist <- function (hsamples, fun) {
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  apply(inf2, c(1,2,3), mean)
}

#The below function, post.predict.dmc.MATCHORDER, accepts a samples object and
#an okdats object. The latter must have the FULL original data frame for each
#participant, before cleaning and including non-responses.
# The function simulates the same amount of data as in the data frame with the same
#design. Then it arranges the orders of trials to correspond to the data.
# This will allow subsequent functions both to truncate simulations that would
#be non-responses, and to calculate the non-response rate.
#
# samples=samples[[1]];n.post=100;probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=FALSE;factors=NA
# save.simulation.as.attribute=FALSE;ignore.R2=FALSE
# gglist=FALSE; probs.gglist=c(0.1, 0.5, 0.9);CI.gglist=c(0.025, 0.975)
#assumes samples have an attribute called NRdata. This attribute is the unfiltered data for each
#sample i.e. before we removed non-responses, cleaned etc.
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

    #getting the factor structure of whatever design you feed in
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


#this function merely lapplys the function above.
h.post.predict.dmc.MATCHORDER <- function(hsamples, n.post=100) {
  lapply(hsamples, post.predict.dmc.MATCHORDER, save.simulation=TRUE, n.post=n.post)
}

# Gets the cumulative sum of RTs for data.
# assigns trials to the data
#with a loop.  Note this trial index is used for the sim as well.
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


#gets the cumulative sum of RTs for each trial for sim.
#for some reason I needed a data.table trick to make this work
# Relies on input from a data object that has been run through the above function add.trial
# .cumsum.data along with the
# sim object.
add.trial.cumsum.sim <- function (sim, datawithtrialcumsum) {
  data <- datawithtrialcumsum
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

#Calculates non-response rates for a data frame or simulated df for E1-A4
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


#This function takes missrates, calculatted with some miss_fun (
# our miss fun for E1-A4 is get.trials.missed.E1_A4), and puts
# in a nice ggdf containing all information needed for plotting -
#posterior mean, posterior credible intervals. Needs sim and data.

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
#

## Below funciton averages parameters across conditions by grepping out
# from av.posts. So if you set av.posts to match anything containing mean_v,
# it would average all rates together and replace all values with the avg before
# simming. We use it to parse the effects of rates/thresholds on the manifests
# in terms of block and cond.

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA
# av.posts<-av.posts.threscond
avps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                    bw="nrd0",report=10,save.simulation=TRUE,factors=NA, av.posts=c())
  # make list of posterior preditive density, quantiles and response p(robability)
{


  get.dqp <- function(sim,facs,probs,n.post=NA) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})

    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA

    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }

    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]


  cat("Below is how I'm averaging (each row is averaged). If this is wrong, adjust your
      av.posts to grep correctly.")
  for (q in 1:length(av.posts)) print(colnames(posts[, grep(av.posts[q], colnames(posts))]))
  ###tweak to average av.posts
  q=1

  if(length(av.posts)!= 0) {
    ### loop through all the averaged posts
    for (q in 1:length(av.posts)) {

      num.params <- dim(posts[, grep(av.posts[q], colnames(posts))])[2]
      average.params <- rowMeans(posts[, grep(av.posts[q], colnames(posts))])
      posts[, grep(av.posts[q], colnames(posts))] <- matrix(average.params,nrow=length(average.params),ncol=num.params,byrow=F)

    }
  }

  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}


#lapplys the above function on everybody
avps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                          bw="nrd0",
                                     save.simulation=FALSE, av.posts=c())
  # apply lost.predict to each subject
{
  lapply(samples,avps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, av.posts=av.posts)
}

# PPs = E1PP
# fun = block.effects.E1A4
# lower=.025
# upper=.975
get.effects.dmc <- function (PPs, fun = function (x) {mean (x)}, lower=.025, upper=.975) {

  simdata<- do.call(rbind, PPs)
  data <- lapply(PPs, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  nreps=max(PPs[[1]]$reps)


  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")
  ######

  ##calculate effects separately for each rep
  for (j in 1:nreps) {

    currentsim.effects <- simdata[simdata$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }

  ##Get a ggplot df with posterior mean, lower, and upper.
  effects.ggdf <-  t(apply(sim.effects, c(2), function(x) c(mean(x),
                                                quantile(x, probs=c(lower,upper)))))
  effects.ggdf <- data.frame(effects.ggdf)
  effects.ggdf <- effects.ggdf[(!rownames(effects.ggdf) %in% "n.rep"),]
  colnames(effects.ggdf) <- c("mean", "lower", "upper")
  contrast <- rownames(effects.ggdf)
  effects.ggdf$data<-as.vector(data.effects)
  attr(effects.ggdf, "post.effects.samples") <- sim.effects
  effects.ggdf
}


block.effects.E1A4 <- function (currentsim) {

  costccC = NA;costccN = NA
  costnnN = NA;costnnC = NA
  accC = NA; accN = NA
  # nonaccC = NA; nonaccN = NA
  pmacc <- NA
  pmrt <- NA

  pmrt <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn")])


  pmacc <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P"])/
    length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn")])

  costccC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])

  costccN <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="3"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="2"])

  costnnN <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])

  costnnC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="3"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="2"])

  accC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="3"]) -
    length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="2"])

  accN <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="3"]) -
    length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="2"])


  out <- c(pmacc,
           pmrt,
           costccC,costnnC,
           costnnN,costccN,
           # noncostccC,noncostccN,
           # noncostnnN,noncostnnC,
           accC,
           # nonaccC,
           accN
           # nonaccN
  )

  names(out) <- c("PM Accuracy",
                  "PM RT",
                  "RT Cost Conflict","RT Cost Conflict (FA)",
                  "RT Cost Nonconflict","RT Cost Nonconflict (FA)",
                  # "noncostccC","noncostccN",
                  # "noncostnnN","noncostnnC",
                  "Accuracy Cost Conflict",
                  # "nonaccC",
                  "Accuracy Cost Nonconflict"
                  # "nonaccN"
  )
  out

}



#This is your function. I added something to get the aggregated RT of
#all rsponses on PM TRIALS (i.e., conf, nonconf, PM all aggregated) which
#I confusingly called pmcrt rather than pmrt.
cond.effects <- function (currentsim) {

  RTccCA <- NA;RTccNA <- NA;RTnnNA <- NA;RTnnCA <- NA
  RTccCB <- NA;RTccNB <- NA;RTnnNB <- NA;RTnnCB <- NA
  RTccCC <- NA;RTccNC <- NA;RTnnNC <- NA;RTnnCC <- NA
  RTccCD <- NA;RTccND <- NA;RTnnND <- NA;RTnnCD <- NA

  RTdiffccCAB = NA;RTdiffccNAB = NA
  RTdiffnnNAB = NA;RTdiffnnCAB = NA

  RTdiffccCBC = NA;RTdiffccNBC = NA
  RTdiffnnNBC = NA;RTdiffnnCBC = NA

  RTdiffccCCD = NA;RTdiffccNCD = NA
  RTdiffnnNCD = NA;RTdiffnnCCD = NA

  accCA <- NA;accCB <- NA;accCC <- NA;accCD <- NA
  accNA <- NA;accNB <- NA;accNC <- NA;accND <- NA

  accdiffCAB = NA; accdiffNAB = NA
  accdiffCBC = NA; accdiffNBC = NA
  accdiffCCD = NA; accdiffNCD = NA
  pmcrtA=NA;pmcrtB=NA;pmcrtC=NA;pmcrtD=NA


  pmaccA <- NA; pmaccB <- NA; pmaccC <- NA; pmaccD <- NA
  pmaccdiffAB <- NA; pmaccdiffBC <- NA; pmaccdiffCD <- NA

  pmaccA <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="A"])/
    length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="A"])
  pmaccB <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="B"])/
    length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="B"])
  pmaccC <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="C"])/
    length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="C"])
  pmaccD <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P" & currentsim$cond=="D"])/
    length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$cond=="D"])

  pmaccdiffAB <- pmaccA - pmaccB
  pmaccdiffBC <- pmaccB - pmaccC
  pmaccdiffCD <- pmaccC - pmaccD

  #
  pmcrtA <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") &  currentsim$cond=="A"])
  pmcrtB <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") &  currentsim$cond=="B"])
  pmcrtC <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") &  currentsim$cond=="C"])
  pmcrtD <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") &  currentsim$cond=="D"])

  RTdiffPMcAB <- pmcrtA -  pmcrtB
  RTdiffPMcBC <- pmcrtB -  pmcrtC
  RTdiffPMcCD <- pmcrtC -  pmcrtD



  # RT for each stim by response by cond
  RTccCA <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])
  RTccNA <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="A"])
  RTnnNA <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])
  RTnnCA <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="A"])

  RTccCB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])
  RTccNB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"])
  RTnnNB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])
  RTnnCB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"])

  RTccCC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])
  RTccNC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"])
  RTnnNC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])
  RTnnCC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"])

  RTccCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])
  RTccND <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="D"])
  RTnnND <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])
  RTnnCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="D"])

  # RT differences between time pressure conditions for each stim by response
  RTdiffccCAB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])
  RTdiffccNAB <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="A"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"])
  RTdiffnnNAB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])
  RTdiffnnCAB <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="A"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"])

  RTdiffccCBC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])
  RTdiffccNBC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="B"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"])
  RTdiffnnNBC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])
  RTdiffnnCBC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="B"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"])

  RTdiffccCCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])
  RTdiffccNCD <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="C"]) -
    mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$cond=="D"])
  RTdiffnnNCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])
  RTdiffnnCCD <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="C"]) -
    mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$cond=="D"])

  # Accuracy by stimulus by condition
  accCA <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="A"])
  accCB <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"])
  accCC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"])
  accCD <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="D"])

  accNA <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="A"])
  accNB <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"])
  accNC <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"])
  accND <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="D"])

  # Accuracy differences between time pressure conditions for each stimulus
  accdiffCAB <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="A"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="A"]) -
    length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"])
  accdiffNAB <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="A"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="A"]) -
    length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"])

  accdiffCBC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="B"]) -
    length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"])
  accdiffNBC <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="B"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="B"]) -
    length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"])

  accdiffCCD <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="C"]) -
    length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$cond=="D"])/
    length(currentsim$RT[currentsim$S=="cc" & currentsim$cond=="D"])
  accdiffNCD <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="C"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="C"]) -
    length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$cond=="D"])/
    length(currentsim$RT[currentsim$S=="nn" & currentsim$cond=="D"])


  RTdiffPMcAB <- pmcrtA -  pmcrtB
  RTdiffPMcBC <- pmcrtB -  pmcrtC
  RTdiffPMcCD <- pmcrtC -  pmcrtD


  out <- c(pmaccA,pmaccB,pmaccC,pmaccD,

           pmaccdiffAB,pmaccdiffBC,pmaccdiffCD,

           RTdiffPMcAB,
           RTdiffPMcBC,
           RTdiffPMcCD ,


           RTccCA,RTccCB,RTccCC,RTccCD,
           RTnnNA,RTnnNB,RTnnNC,RTnnND,
           RTnnCA,RTnnCB,RTnnCC,RTnnCD,
           RTccNA,RTccNB,RTccNC,RTccND,

           RTdiffccCAB,RTdiffccCBC,RTdiffccCCD,
           RTdiffnnNAB,RTdiffnnNBC,RTdiffnnNCD,
           RTdiffnnCAB,RTdiffnnCBC,RTdiffnnCCD,
           RTdiffccNAB,RTdiffccNBC,RTdiffccNCD,

           accCA,accCB,accCC,accCD,
           accNA,accNB,accNC,accND,

           accdiffCAB,accdiffCBC,accdiffCCD,
           accdiffNAB,accdiffNBC,accdiffNCD

           # pmrtdiff,

  )

  names(out) <- c("PM Accuracy A","PM Accuracy B","PM Accuracy C","PM Accuracy D",
                  "PM Acc Diff A-B","PM Acc Diff B-C","PM Acc Diff C-D",
                  "RT Diff PM A-B",
                  "RT Diff PM B-C",
                  "RT Diff PM C-D",

                  "RT Conflict A","RT Conflict B","RT Conflict C","RT Conflict D",
                  "RT Nonconflict A","RT Nonconflict B","RT Nonconflict C","RT Nonconflict D",
                  "RT Conflict (FA) A","RT Conflict (FA) B","RT Conflict (FA) C","RT Conflict (FA) D",
                  "RT Nonconflict (FA) A","RT Nonconflict (FA) B","RT Nonconflict (FA) C","RT Nonconflict (FA) D",

                  "RT Diff Conflict A-B","RT Diff Conflict B-C","RT Diff Conflict C-D",
                  "RT Diff Nonconflict A-B","RT Diff Nonconflict B-C","RT Diff Nonconflict C-D",
                  "RT Diff Conflict (FA) A-B","RT Diff Conflict (FA) B-C","RT Diff Conflict (FA) C-D",
                  "RT Diff Nonconflict (FA) A-B","RT Diff Nonconflict (FA) B-C","RT Diff Nonconflict (FA) C-D",

                  "Accuracy Conflict A","Accuracy Conflict B","Accuracy Conflict C","Accuracy Conflict D",
                  "Accuracy Nonconflict A","Accuracy Nonconflict B","Accuracy Nonconflict C","Accuracy Nonconflict D",

                  "Acc Diff Conflict A-B","Acc Diff Conflict B-C","Acc Diff Conflict C-D",
                  "Acc Diff Nonconflict A-B","Acc Diff Nonconflict B-C","Acc Diff Nonconflict C-D"

  )
  out

}


#The below function picks certain parameters from a vector called pickps_other,
# and replaces them with pickps_set, before performing posterior prediciton.
# We use it to turn control mechanisms off in the model. To turn off proactive
# control, we set the ongoing task thresholds equal in the PM block to the control
#threshols. To turn off reactive, we set the ongiong rates on PM trials (PM block)
# to the ongoing rates on non-PM trials (PM block)

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA; n.post=100
# pickps_set <- c("B.A2C",        "B.B2C",        "B.C2C" ,
#                 "B.D2C"  ,             "B.A2N"  ,      "B.B2N"  ,      "B.C2N"  ,
#                 "B.D2N")
#
# pickps_others <- c("B.A3C"   ,     "B.B3C",        "B.C3C" ,
#                   "B.D3C",        "B.A3N"   ,     "B.B3N" ,       "B.C3N" ,
#                   "B.D3N" )
#
#


pickps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                   bw="nrd0",report=10,save.simulation=TRUE,factors=NA, pickps_others, pickps_set)
  # make list of posterior preditive density, quantiles and response p(robability)
{
  
  
  get.dqp <- function(sim,facs,probs,n.post=NA) {
    
    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    
    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA
    
    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }
    
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  
  #more robust
  
  
  ###Replace some parameter vlaues with others.
  posts[,colnames(posts) %in% pickps_others][,pickps_others] <- 
    posts[,colnames(posts) %in% pickps_set][,pickps_set] 
  
  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}

#lapply the above to the whole samples object.
pickps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                   bw="nrd0",
                                   save.simulation=FALSE, pickps_set, pickps_others)
  # apply lost.predict to each subject
{
  lapply(samples,pickps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, pickps_set=pickps_set, pickps_others=pickps_others)
}


## After you have collected up the effects from a posterior sim
# with get.effects.dmc, this function will process the data to get
# block relevant effects for E1_A4.
finish.blockdf.E1_A4 <- function(effects) {
  effects$S <- NA
  effects$DV <- NA
  effects$S[grep ("Conflict", rownames(effects))] <- "Conflict"
  effects$S[grep ("Nonconflict", rownames(effects))] <- "Nonconflict"
  effects$DV <- "Accuracy"
  effects$DV[grep ("RT", rownames(effects))] <- "Correct RT"
  effects$DV[grep ("PM RT", rownames(effects))] <- "RT"
  effects$DV[grep ("(FA)", rownames(effects))] <- "Error RT"
  effects$S[grep ("PM", rownames(effects))] <- "PM"
  effects
}

## After you have collected up the effects from a posterior sim
# with get.effects.dmc, this function will process the data to get
# cond relevant effects for E1_A4.
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
  effects$DV[grep ("RT Diff PM", rownames(effects))] <- "RT"
  effects<- effects[grepl("Diff", rownames(effects)),]
  effects$contrast[grep ("A-B", rownames(effects))] <- "A-B"
  effects$contrast[grep ("B-C", rownames(effects))] <- "B-C"
  effects$contrast[grep ("C-D", rownames(effects))] <- "C-D"
  effects
}




####  Unused functions ######
#redundant stuff that I'm not entirely sure we won't need again.

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


# finish.blockdf.E1_A4 <- function(effects) {
#   effects$S <- NA
#   effects$DV <- NA
#   effects$S[grep ("Conflict", rownames(effects))] <- "Conflict"
#   effects$S[grep ("Nonconflict", rownames(effects))] <- "Nonconflict"
#   effects$DV <- "Accuracy"
#   effects$DV[grep ("RT", rownames(effects))] <- "Correct RT"
#   effects$DV[grep ("eRT", rownames(effects))] <- "Error RT"
#   effects$DV[grep ("(FA)", rownames(effects))] <- "Error RT"
#   effects$S[grep ("PM", rownames(effects))] <- "PM"
#   effects
# }



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
#
# get.trials.missed.E1_A4 <- function (sim) {
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



# block.effects.E1A4 <- function (currentsim) {
#
#   costccC = NA;costccN = NA
#   costnnN = NA;costnnC = NA
#   accC = NA; accN = NA
#   # nonaccC = NA; nonaccN = NA
#   pmacc <- NA
#   pmcrt <- NA
#   pmert <- NA
#
#   pmcrt <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P"])
#   pmert <- mean(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & !currentsim$R=="P"])
#
#   pmacc <- length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn") & currentsim$R=="P"])/
#     length(currentsim$RT[(currentsim$S=="pc"|currentsim$S=="pn")])
#
#   costccC <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"]) -
#     mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])
#
#   costccN <- mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="3"]) -
#     mean(currentsim$RT[currentsim$S=="cc" & currentsim$R=="N" & currentsim$block=="2"])
#
#   costnnN <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"]) -
#     mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])
#
#   costnnC <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="3"]) -
#     mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="C" & currentsim$block=="2"])
#
#   accC <- length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="3"])/
#     length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="3"]) -
#     length(currentsim$RT[currentsim$S=="cc" & currentsim$R=="C" & currentsim$block=="2"])/
#     length(currentsim$RT[currentsim$S=="cc" & currentsim$block=="2"])
#
#   accN <- length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="3"])/
#     length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="3"]) -
#     length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$block=="2"])/
#     length(currentsim$RT[currentsim$S=="nn" & currentsim$block=="2"])
#
#
#   out <- c(pmacc,
#            pmcrt,
#            pmert,
#            costccC,costnnC,
#            costnnN,costccN,
#            # noncostccC,noncostccN,
#            # noncostnnN,noncostnnC,
#            accC,
#            # nonaccC,
#            accN
#            # nonaccN
#   )
#
#   names(out) <- c("PM Accuracy",
#                   "PM cRT",
#                   "PM eRT",
#                   "RT Cost Conflict","RT Cost Conflict (FA)",
#                   "RT Cost Nonconflict","RT Cost Nonconflict (FA)",
#                   # "noncostccC","noncostccN",
#                   # "noncostnnN","noncostnnC",
#                   "Accuracy Cost Conflict",
#                   # "nonaccC",
#                   "Accuracy Cost Nonconflict"
#                   # "nonaccN"
#   )
#   out
#
# }

##LUKE additions 04/08 at Russ request- so we can sim proportion of NRs with parameters
#avgd


post.predict.dmc.MATCHORDER.AVG <- function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                            bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
                                            save.simulation.as.attribute=FALSE,ignore.R2=FALSE,
                                            gglist=FALSE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975), av.posts=c())
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
  ###Plug in averaging 
  cat("Below is how I'm averaging (each row is averaged). If this is wrong, adjust your
      av.posts to grep correctly.")
  for (q in 1:length(av.posts)) print(colnames(posts[, grep(av.posts[q], colnames(posts))]))
  ###tweak to average av.posts
  q=1
  
  if(length(av.posts)!= 0) {
    ### loop through all the averaged posts
    for (q in 1:length(av.posts)) {
      
      num.params <- dim(posts[, grep(av.posts[q], colnames(posts))])[2]
      average.params <- rowMeans(posts[, grep(av.posts[q], colnames(posts))])
      posts[, grep(av.posts[q], colnames(posts))] <- matrix(average.params,nrow=length(average.params),ncol=num.params,byrow=F)
      
    }
  }
  ######  
  
  
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
    
    #getting the factor structure of whatever design you feed in
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






#this function merely lapplys the function above.
h.post.predict.dmc.MATCHORDER.AVG <- function(hsamples, n.post=100, av.posts=c()) {
  lapply(hsamples, post.predict.dmc.MATCHORDER.AVG, save.simulation=TRUE, n.post=n.post, av.posts=av.posts)
}







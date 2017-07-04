library(car)
library(heplots)

wsAnova=function(dat,SStype=3,spss=F) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  dat=dat[do.call(order,dat[,-dim(dat)[2]]),]
  snams=levels(dat[,1]); ns=length(snams)
  dvnam=names(dat)[dim(dat)[2]]
  facn=names(dat)[-c(1,dim(dat)[2])]
  nifac=length(facn)
  idata=data.frame(dat[dat[,1]==snams[1],facn])
  names(idata)=facn
  for (i in facn) 
    if (i==facn[1]) ifacn=as.character(idata[,1]) else 
      ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
  facnr=facn[nifac:1]
  e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
              ncol=length(ifacn),dimnames=list(snams,ifacn))
  print(summary(Anova(lm(e.mv ~ 1),
                idata=idata,type=SStype, 
                idesign=formula(paste("~",paste(facn,collapse="*")))),
          multivariate=FALSE))
#LUKE PARTIAL ETA SQUARED  
  SSH <- (summary(Anova(lm(e.mv ~ 1),
                        idata=idata,type=SStype, 
                        idesign=formula(paste("~",paste(facn,collapse="*")))),
                  multivariate=FALSE))$univariate.tests[,1] 
  
  SST <- SSH + (summary(Anova(lm(e.mv ~ 1),
                              idata=idata,type=SStype, 
                              idesign=formula(paste("~",paste(facn,collapse="*")))),
                        multivariate=FALSE))$univariate.tests[,3] 

  out <- append(summary(Anova(lm(e.mv ~ 1),
                idata=idata,type=SStype, 
                idesign=formula(paste("~",paste(facn,collapse="*")))),
          multivariate=FALSE), list((SSH/SST)[-1]))
 names(out)[length(out)] <- "PartialEtaSq"
 print(out)
}


  #print("Partial Eta Sq:")
#  print((SSH/SST)[-1])
###  
  #if (spss) {
 #   e.mv=cbind.data.frame(s=row.names(e.mv),e.mv)
 #   row.names(e.mv)=NULL
 #   e.mv
#  }

mneffects=function(df,elist,digits=3,err=F,vars=F,dvnam="y") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  Tmns <- list()
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  for (i in 1:length(elist)) {
    cat(paste(paste(elist[[i]],collapse=":"),"\n"))
    mns=tapply(df[,dvnam],df[,elist[[i]]],mean,na.rm=T)
    if (err) print(round(plogis(mns),digits)) else
      if (vars) print(round(sqrt(mns),digits))  else
        print(round(mns,digits)) 
Tmns[[i]]<- mns; names(Tmns)[i] <- paste(elist[[i]],collapse=":")
    cat("\n")
  }  
  Tmns
}


se=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  mn=tapply(df[,dvnam],df[,facs],mean)
  se=tapply(df[,dvnam],df[,facs],sd)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (is.na(ci)) mn else {
    if (ci=="SE") se/sqrt(ns) else
     qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
  }
}

se2=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean, na.rm=T)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  mn=tapply(df[,dvnam],df[,facs],mean, na.rm=T)
  se=tapply(df[,dvnam],df[,facs],sd, na.rm=T)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (is.na(ci)) mn else {
    if (ci=="SE") se/sqrt(ns) else
     qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
  }
}

add.bars=function(mn,se,xvals=NA,len=.1,antiprobit=FALSE,col="black") {
  
  plotbars <- function(x,m,l,h,len,col="black") {
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],l[j],length=len,angle=90,col=col)
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],h[j],length=len,angle=90,col=col)    
  }
  
  if (any(is.na(xvals))) if (is.matrix(mn)) 
    xvals <- as.numeric(dimnames(mn)[[2]]) else
    xvals <- as.numeric(factor(names(mn)))
  lo <- mn-se
  hi <- mn+se
  if (antiprobit) {
    mn=pnorm(mn)
    lo=pnorm(lo)
    hi=pnorm(hi)
  }
  if (!is.matrix(mn)) 
      plotbars(xvals,mn,lo,hi,col=col,len=len) else
    for (i in 1:dim(mn)[1]) 
      plotbars(x=xvals,m=mn[i,],l=lo[i,],h=hi[i,],len=len,col=col)
}    

arr2df=function(arr) {
  if (is.null(dim(arr))) out=data.frame(y=arr) else {
    dn=dimnames(arr)
    if (length(dn)==1) {
      out=cbind.data.frame(factor(dn[[1]],dn[[1]]),arr)
      names(out)=c(names(dn),"y")
      row.names(out)=NULL
    } else {
      tmp=vector(mode="list",length=length(dn))
      names(tmp)=names(dn)
      k=1
      for (j in names(dn)) {
        n=length(dn[[j]])
        tmp[[j]]=gl(n,k,length(arr),dn[[j]])
        k=k*n
      }
      out=cbind(data.frame(tmp),y=as.vector(arr))
      row.names(out)=NULL
    }
  }
  out
}

# spss=F; dvnam=NA; SStype=3; sfac="s"
mixedAnova=function(dat,bsfacn,wsfacn=NULL,sfac="s",SStype=3,spss=F,dvnam=NA) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  if (is.na(dvnam)) dvnam <- names(dat)[dim(dat)[2]]
  dat <- dat[,c(sfac,bsfacn,wsfacn,dvnam)]
  if (length(wsfacn)>0) dat <- dat[do.call(order,dat[,c(sfac,wsfacn)]),]
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  snams=levels(dat[,sfac]); ns=length(snams)
  nifac=length(wsfacn)
  lev1s=unlist(lapply(lapply(dat[,wsfacn],levels),function(x){x[1]}))
  bsfacs=dat[apply(dat[,wsfacn,drop=F],1,function(x){all(x==lev1s)}),bsfacn,drop=F]  
  if ( nifac>0 ) {
    idata=data.frame(dat[dat[,sfac]==snams[1],wsfacn])
    names(idata)=wsfacn
    for (i in wsfacn) 
      if (i==wsfacn[1]) ifacn=as.character(idata[,1]) else 
        ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
    e.mv=matrix(unlist(
      tapply(dat[,dvnam],dat[,wsfacn[length(wsfacn):1]],function(x){x})),
                ncol=length(ifacn),dimnames=list(snams,ifacn))
    summary(Anova(
      lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),bsfacs),
      idata=idata,type=SStype,
      idesign=formula(paste("~",paste(wsfacn,collapse="*")))),multivariate=F)
  } else {
    e.mv <- cbind(y=dat[,dvnam]) 
    print(Anova(lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),
                   bsfacs),type=3))
  }
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),bsfacs,e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}

# scol = subject column (numeric or name, same for following)
# dvcol = dependent variable column
# ws and bs name corresponding factors in dat
# mixedAnova=function(dat,ws,bs=NA,scol=1,dvcol=dim(dat)[2],
#   SStype=3,icontrasts=c("contr.sum", "contr.poly"),
#   savedf=F) {
#   dat=dat[do.call(order,dat[,ws,drop=F]),]
#   snams=levels(dat[,scol]); ns=length(snams)
#   dvnam=names(dat)[dvcol]
#   nifac=length(ws)
#   idata=data.frame(dat[dat[,scol]==snams[1],ws])
#   names(idata)=ws  
#   for (i in ws) 
#     if (i==ws[1]) ifacn=as.character(idata[,1]) else 
#                   ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
#   facnr=ws[nifac:1]
#   e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
#               ncol=length(ifacn),dimnames=list(snams,ifacn))
#   if (length(bs)==1) bsf=bs else bsf=paste(bs,collapse="*")
#   if (any(is.na(bs))) {
#     form=formulua("e.mv~1")
#     tmp <- e.mv
#   } else {
#     form=formula(paste("e.mv ~",bsf))
#     trans=snams
#     names(trans)=snams
#     for (i in bs) {
#       for (j in snams) trans[j] <- as.character(dat[dat[,scol]==j,i][1])
#       if (i==bs[1]) bsfac <- factor(trans[dimnames(e.mv)[[1]]]) else
#         bsfac <- cbind.data.frame(bsfac,factor(trans[dimnames(e.mv)[[1]]]))
#     }
#     tmp <- cbind.data.frame(bsfac,e.mv)
#     names(tmp)[1:length(bs)] <- bs
#   }
#   summary(Anova(lm(form,tmp),
#     idata=idata,type=SStype,icontrasts=icontrasts,
#     idesign=formula(paste("~",paste(ws,collapse="*"))) ),
#     multivariate=FALSE)
#   invisible(tmp)
# }

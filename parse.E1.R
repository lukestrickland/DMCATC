rm(list=ls())
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
# setwd("~/Software/DMCATC")
pkgs <- c("plyr", "dplyr", "tidyr", "broom", "pander", "xtable")
# install.packages(pkgs) #install
sapply(pkgs, require, character.only = T) #load

# Load data
data <- read.csv("data/exp_data/ATC.E1.csv")
str(data)

## Exclude p17 due to no PM responses
data <- data[ data$s!="p17", ]
data$s <- factor(data$s)

## PARSE E1

dat <- data
dat$s <- factor(as.character(dat$s))
# dat$TP <- dat$TP/1000
dat$RT <- dat$RT/1000

str(dat)
dat$pool <- NULL
dat$load <- NULL
dat$TP <- NULL
dat$presOrder <- NULL
dat$pairNumber <- NULL
dat$DOMS <- NULL
dat$TTMS <- NULL
dat$OOP <- NULL
dat$C <- NULL
dat$angle <- NULL
dat$ac1_type <- NULL
dat$ac2_type <- NULL
dat$ac1_cs <- NULL
dat$ac2_cs <- NULL
dat$ac1_fl <- NULL
dat$ac2_fl <- NULL
dat$ac1_speed <- NULL
dat$ac2_speed <- NULL
dat$pm_status <- NULL
dat$conflict_status <- NULL
str(dat)

# Rename stimulus and response factor levels
levels(dat$S)
levels(dat$S) <- c("cc", "nn", "pc", "pn")

levels(dat$R)
levels(dat$R) <- c("C","N","P","M")

# Revalue Control and PM blocks by number of accumulators ("2", "3")
levels(dat$block)
levels(dat$block) <- c("2", "3")

# Label nonresponses with 'M' for 'miss'
dat[ which(is.na(dat$R)==TRUE), "R" ] <- "M"
plot(dat$R)

# Calculate % nonresponses
percent.NR <- length(dat$R[dat$R=="M"])/length(dat$R) * 100
percent.NR

# Calculate % PM responses in control blocks
# dat[dat$block=="2" & dat$R=="P",]
percent.PM.in.Control <- length(dat[dat$block=="2" & dat$R=="P", 1])/length(dat$R) * 100
percent.PM.in.Control

# Clean data, remove outliers and nonresponses
okdats <- dat
okdats <- okdats[okdats$R!="M", ]  # Remove misses
okdats$R <- factor(okdats$R)
okdats <- okdats[!is.na(okdats$RT),]  # Remove NA RT values
str(okdats)

# Remove PM false alarms
okdats<-okdats[!(okdats$block=="2" & okdats$R=="P"),]
str(okdats)

clean <- function(df) {
    dfc <- df
    n=tapply(df$RT,list(df$s),length)
    ns=tapply(df$RT,list(df$s),length)
    mn=tapply(df$RT,list(df$s),mean)
    sd=tapply(df$RT,list(df$s),IQR)
    upper <- mn+3*(sd/1.349)
    lower <- 0.2
    bad <- logical(dim(df)[1])
    levs <- paste(df$s,sep=".")
    for (i in levels(df$s)){
        lev <- i
        bad[levs==lev] <- df[levs==lev,"RT"] > upper[i] | df[levs==lev,"RT"] < lower
    }
    df=df[!bad,]
    nok=tapply(df$RT,list(df$s),length)
    pbad=100-100*nok/n
    nok=tapply(df$RT,list(df$s),length)
    pbad=100-100*nok/ns
    print(sort(round(pbad,5)))
    print(mean(pbad,na.rm=T))
    df
}


okdats <- clean(okdats)
str(okdats)

# Calculate % of responses cleaned / outliers removed
percent.cleaned <- (1 - length(okdats$RT)/length(dat$RT)) * 100
percent.cleaned



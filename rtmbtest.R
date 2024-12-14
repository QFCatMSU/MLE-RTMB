# Simple test script for RTMB built to do many of the calculations 
# described at: https://kaskr.r-universe.dev/articles/RTMB/RTMB-introduction.html

parameters <- list(
  mua=0,          ## Mean slope
  sda=1,          ## Std of slopes
  mub=0,          ## Mean intercept
  sdb=1,          ## Std of intercepts
  sdeps=1,        ## Residual Std
  a=rep(0, 50),   ## Random slope by chick
  b=rep(0, 50)    ## Random intercept by chick
)

f <- function(parms) {
  getAll(ChickWeight, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  weight <- OBS(weight)
  ## Initialize joint negative log likelihood
  nll <- 0
  ## Random slopes
  nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
  ## Random intercepts
  nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
  ## Data
  predWeight <- a[Chick] * Time + b[Chick]
  nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
  ## Get predicted weight uncertainties
  ADREPORT(predWeight)
  ## Return
  nll
}

obj <- MakeADFun(f, parameters, random=c("a", "b"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr

as.list(sdr, "Est") ## parameter estimates
as.list(sdr, "Std") ## parameter uncertainties

as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
as.list(sdr, "Std", report=TRUE) ## ADREPORT uncertainties

set.seed(1)
chk <- checkConsistency(obj)
chk

osa <- oneStepPredict(obj, method="fullGaussian", discrete=FALSE)
qqnorm(osa$res); abline(0,1)


# Modified from Cahill RTMB GitHub
# Bence July 2023, cahill added a bit of styling
# 2024 Simulates data or could still read in data used in past
# R script for fitting vonB model to multiple ponds with TMB
# Allowing for random effects on Linf

library(RTMB)

#Simulate data *************************************
simpar = list(
  logLinfmn = 6.38, logLinfsd = 0.03,
  logK = -1.63, t0 = -1.04, Sig = 9.0
)
ages=2:12 
nponds=10
simdat =function(p,ages,nponds){
  logLinfs=rnorm(nponds,p$logLinfmn,p$logLinfsd)
  Lens=matrix(nrow=length(ages),ncol=nponds)
  K=exp(p$logK)
  Linfs=exp(logLinfs)
  for (j in (1:nponds)){
    ExL = Linfs[j] * (1 - exp(-K * (ages - p$t0)))
    Lens[,j]=rnorm(length(ages),mean=ExL,sd=p$Sig)
  }
  Lens
}

set.seed(1234)
mvonBRdat =as.data.frame(simdat(simpar,ages,nponds))
#End Simulate data **********************************

#mvonBRdat = read.table("lesson4/data/multiLinfTMB.dat", header = F)

# Set up the data and starting value of parameters for TMB
data = list(L = as.matrix(mvonBRdat), A = 2:12)

Linfs = sapply(mvonBRdat, max)
pars = list(
  logLinfmn = 7, loglogLinfsd = -2.8,
  logLinfs = log(Linfs),
  logK = -1.6, t0 = 0, logSig = 4
)

f = function(pars) {
  getAll(data, pars)
  Linfmn = exp(logLinfmn)
  logLinfsd = exp(loglogLinfsd)
  Linfs = exp(logLinfs)
  K = exp(logK)
  Sig = exp(logSig)
  nponds = length(Linfs)
  nages = length(A)
  predL = matrix(0, nrow = nages, ncol = nponds)
  # fill one column (pond) at a time:
  for (i in 1:nponds) {
    predL[, i] = Linfs[i] * (1 - exp(-K * (A - t0)))
  }
  nll = -sum(dnorm(x = L, mean = predL, sd = Sig, log = TRUE))
  nprand = -sum(dnorm(x = logLinfs, mean = logLinfmn, sd = logLinfsd, log = TRUE))
  jnll = nll + nprand
  jnll
}

obj = MakeADFun(f, pars, random = c("logLinfs"))
fit = nlminb(obj$par, obj$fn, obj$gr)
fit 

sdr = sdreport(obj)
sdr



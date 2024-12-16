# R script for fitting musky L@Age with RTMB
# modification for demo of passing data using closure
# 12/2024 JRB

library(RTMB)

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T)

#Set up the data and starting value of parameters for RTMB
dat = list(lenobs=gmRdat[,"Length"],age=gmRdat[,"Age"])
par = list(logLinf=7,logK=-1.6,t0=0,logsd=4)

f = function(par,dat){
  getAll(par,dat)
  Linf = exp(logLinf)
  K = exp(logK)
  sd = exp(logsd)
  lenobs=OBS(lenobs) # New line!
  lenpred = Linf * (1 - exp(-K * (age - t0)))
  nll = -sum(dnorm(lenobs, lenpred, sd, TRUE))
  atagepred = Linf * (1 - exp(-K * ((1:11) - t0)))
  REPORT(atagepred)
  nll
}

#cmb is a closure.  Takes two arguments but returns function with one arg
cmb <- function(f, d) function(p) f(p, d)

obj=RTMB::MakeADFun(cmb(f,gmdat),par)
fit=nlminb(obj$par, obj$fn, obj$gr)
fit$convergence
sdr=sdreport(obj)
sdr

simdat=gmdat;  #copy real data
#Simulate obs lengths, write to data copy,refit model
set.seed(123456)
simdat$lenobs=obj$simulate()$lenobs

objsim = RTMB::MakeADFun(cmb(f,simdat),par)
simfit = nlminb(objsim$par, objsim$fn, objsim$gr)
simfit$convergence
simsdr=sdreport(objsim)
simsdr



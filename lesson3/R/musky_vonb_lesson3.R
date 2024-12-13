# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
# various modifications for lesson 3, Dec 2023/2023 workshops
library(RTMB);

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T);

#Set up the data and starting value of parameters for RTMB
datlst = list(lenobs=gmRdat[,"Length"],age=gmRdat[,"Age"]);
parlst = list(loglinf=7,logvbk=-1.6,t0=0,logsd=4);

#My NLL function
f = function(parlst){
  getAll(datlst,parlst);
  linf = exp(loglinf);
  vbk = exp(logvbk);
  sd = exp(logsd);
  lenpred = linf * (1 - exp(-vbk * (age - t0)));
  nll = -sum(dnorm(lenobs, lenpred, sd, TRUE));
  atagepred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atagepred);
  nll
}

obj <- MakeADFun(f,parlst); #RTMB::MakeADFun(.) if TMB loaded
fit = nlminb(obj$par, obj$fn, obj$gr);
sdr = sdreport(obj);
sdr

#refit with bounds
lower = c(4,-5,-5,0);
upper = c(10,1,5,10);
objb <- MakeADFun(f,parlst);
fitb = nlminb(objb$par, objb$fn, objb$gr,lower=lower,upper=upper);
sdrb = sdreport(objb);
sdrb

#refit with log_linf fixed
mymap=list(loglinf=factor(NA));
objf <- MakeADFun(f,parlst,map=mymap);
fitf = nlminb(objf$par, objf$fn, objf$gr);
sdrf = sdreport(objf);
sdrf

#refit with both bounds and log_linf fixed
#note removal of bounds for log_linf
lower = c(-5,-5,0);
upper = c(1,5,10);
objfb <- MakeADFun(f,parlst,map=mymap);
fitfb = nlminb(objfb$par, objfb$fn, objfb$gr,lower=lower,upper=upper);
sdrfb = sdreport(objfb);
sdrfb

#Do a likelihood profile CI for logLinf
library(TMB); #load for CI
loglinfprof = tmbprofile(obj,"loglinf");
confint(loglinfprof);
plot(loglinfprof);



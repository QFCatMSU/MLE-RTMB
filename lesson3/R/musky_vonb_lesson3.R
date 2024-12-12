# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
# various modifications for lesson 3, Dec 2023 workshop
library(RTMB);

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T);

#Set up the data and starting value of parameters for RTMB
gmdat = list(len_obs=gmRdat[,"Length"],age=gmRdat[,"Age"]);
gmpar = list(log_linf=7,log_vbk=-1.6,t0=0,log_sd=4);

#My NLL function
NLL_fun = function(par_lst){
  getAll(gmdat,par_lst);
  linf = exp(log_linf);
  vbk = exp(log_vbk);
  sd = exp(log_sd);
  len_pred = linf * (1 - exp(-vbk * (age - t0)));
  nll = -sum(dnorm(len_obs, len_pred, sd, TRUE));
  atage_pred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atage_pred);
  nll
}

obj <- MakeADFun(NLL_fun,gmpar); #RTMB::MakeADFun(.) if TMB loaded
fit = nlminb(obj$par, obj$fn, obj$gr);
sdr = sdreport(obj);
sdr

#refit with bounds
lower = c(4,-5,-5,0);
upper = c(10,1,5,10);
objb <- MakeADFun(NLL_fun,gmpar);
fitb = nlminb(objb$par, objb$fn, objb$gr,lower=lower,upper=upper);
sdrb = sdreport(objb);
sdrb

#refit with log_linf fixed
mymap=list(log_linf=factor(NA));
objf <- MakeADFun(NLL_fun,gmpar,map=mymap);
fitf = nlminb(objf$par, objf$fn, objf$gr);
sdrf = sdreport(objf);
sdrf

#refit with both bounds and log_linf fixed
#note removal of bounds for log_linf
lower = c(-5,-5,0);
upper = c(1,5,10);
objfb <- MakeADFun(NLL_fun,gmpar,map=mymap);
fitfb = nlminb(objfb$par, objfb$fn, objfb$gr,lower=lower,upper=upper);
sdrfb = sdreport(objfb);
sdrfb

#Do a likelihood profile CI for logLinf
library(TMB); #load for CI
log_linf_prof = tmbprofile(obj,"log_linf");
confint(log_linf_prof);
plot(log_linf_prof);



# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
library(RTMB);

gmRdat = read.table("lesson2/data/musky_vonb.dat",head=T);

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
  atagepred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atagepred);
  -sum(dnorm(lenobs, lenpred, sd, TRUE));
}

obj <- MakeADFun(f,parlst);

#Print report variables before fitting the model
GMreport=obj$report();
GMreport

## Call function minimizer
fit = nlminb(obj$par, obj$fn, obj$gr);

## Get parameter uncertainties and convergence diagnostics
sdr = sdreport(obj);
summary(sdr)

#Save the report variable list and print it after model fitting 
GMreport = obj$report();
GMreport;


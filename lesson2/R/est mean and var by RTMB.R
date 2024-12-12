# script for RTMB search for the mean and log SD

#load RTMB
library(RTMB);

#data
xvec=c(10.72,7.23,10.07,8.62,8.55);
datlst=list(xvec=xvec)

#starting guesses for parameters
parlst=list(mu=8,logsd=log(sqrt(3)));

#My NLL function
f = function(parlst){
  getAll(datlst,parlst)
  -sum(dnorm(xvec,mean=mu,sd=exp(logsd),log=T))
}

obj = MakeADFun(f,parlst);
opt = nlminb(obj$par,obj$fn,obj$gr);
sdrep = sdreport(obj)
sumsd = summary(sdrep)
sumsd

#calculate var est from log SD estimate to compare
exp(sumsd[2,1])^2

#calculate the analytical param estimates
mean(xvec);
var(xvec);

#var gives the unbiased rather than MLE est of var
#convert to MLE (see bias example from lecture)
n=length(xvec);
var(xvec)*(n-1)/n;

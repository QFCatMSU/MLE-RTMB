# script for search for MLE of detection probability
#uses nlminb without and with RTMB

#load RTMB
library(RTMB);


#data
n=30
y=2

#starting guess for parameter
p=0.1;
tranp=log(p/(1-p))#logit (inverse of logistic function)

#This version does not use RTMB! Gets derivs by finite difference
f= function(tranp){
  prob=1/(1+exp(-tranp)) #back transform detection prob using logistic function
  -dbinom(x=y,size=n,prob=prob,log=TRUE)
}

simpest<- nlminb(start=tranp,f)

#RTMB version
parlst=list(tranp=tranp)
datlst=list(n=n,y=y)

#My RTMB NLL function
f = function(parlst){
  getAll(datlst,parlst)
  prob=1/(1+exp(-tranp)) #back transform detection prob using logistic function
    -dbinom(x=y,size=n,prob=prob,log=TRUE)
}

obj = MakeADFun(f,parlst);
opt = nlminb(obj$par,obj$fn,obj$gr);
sdrep = sdreport(obj)
summary(sdrep)
opt






 #RTMB program to fit hierarchical vonB for "realistic" data
 # Created Feb 2021 in TMB for MLE Software course by JRB
 # Model with params Linf, K, L2 (rather than t0)
 # {logLinf, logK, logL2} MV Normal
 # Modified July-Dec 2023 for RTMB
 # Some styling and incorporation of data generation Dec 2024

library(RTMB);


#******** Data simulation
set.seed(100)
#vonB parameterized with Linf, K, L2 rather than t0)
#"pond" specific vonB parameters logLinf, logK and LogL2 MVN

#simulated variation about pond specific vonB is normal
# with logCV = intercept + slope*(expected length at age)

library(MASS);

#Mean over ponds vonB params
logLinfBar = log(500)
logKBar = log(.35)
logL2Bar = log(175)

#var-cov matrix Sigma for vonB params (variation over ponds)
SD = diag(c(0.3,0.2,0.15))
corvb = matrix(c(1.0 ,-0.5, -0.4,
                 -0.5, 1.0, 0.4,
                 -0.4, 0.4, 1.0), nrow=3,ncol=3)
Sigma = SD%*%corvb%*%SD

#generate pond-specific sample sizes, Z mortality rate, vonB params
nponds = 20
nfish=round(runif(n = nponds, min = 100, max = 500))
Z = rlnorm(n=nponds,meanlog = log(.45), sdlog = .2);
vbPond=mvrnorm(n = 20, mu = c(logLinfBar,logKBar,logL2Bar), 
               Sigma=Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

#log CV for among fish variation in length is linear function of expect L given age
CVInt = log(0.4)
CVSlp = -(log(0.4)-log(0.05))/500
#Note alt version just uses SDL
#SDL=20;

res = list()
for (p in 1:20){
  #age specific n Poisson with expectation based on pond-specific expected 
  # sample size (nfish) and pond-specific Z
  Na = rpois(n = length(2:25), nfish[p]*exp(-Z[p]*(2:25))/sum(exp(-Z[p]*(2:25))))
  #crude emulation of capping N by size category to avoid oversampling ages
  Na = round(pmin(Na,runif(length(Na), min = 6.5, max = 13.5)))
  Linf = exp(vbPond[p,1])
  K = exp(vbPond[p,2])
  L2=exp(vbPond[p,3])
  tmplst = list()
  for (a in 2:25){
    mnL=L2+(Linf-L2)*(1-exp(-K*(a-2)));
    sd=exp(CVInt+CVSlp*mnL)*mnL; #log CV linear function of predicted length
    #  sd=SDL;
    if(Na[a-1]>0) {
      tmplst[[a-1]]=cbind(pid=rep(p),age=a,L=round(rnorm(n=Na[a-1],mean=mnL,sd=sd)));
    }
  }
  res=c(res,tmplst);
}   

res2=data.frame(do.call(rbind,res));
#write.csv(x=res2,file = "newmulti.csv", row.names = FALSE);
#******** End Data simulation


#Data visualization

plot (res2$L~res2$age)

for(i in levels(as.factor(res2$pid))) {
  plot(res2$L[res2$pid==i]~res2$age[res2$pid==i],main=paste0("pid = ",i)) 
}





dat=list(pid = res2$pid,obsL = res2$L, age = res2$age)
 #find max length for each pond (use as starting Linf)
 maxL<-aggregate(L ~ pid, data = res2, max)
 
 pars = list(logLinfBar = log(500), logLinf = log(maxL$L), 
               logKBar = log(.3), logK = rep(log(.3),20),
               logL2Bar = log(175), logL2 = rep(log(175),20), 
               logvbsd=log(c(0.1,0.1,0.1)),
               theta = c(0,0,0), logcvInt = log(0.1), logcvSlp=0);

us = unstructured(3); 

f = function(pars){
 getAll(dat, pars)
 LinfBar = exp(logLinfBar)
 KBar = exp(logKBar)
 L2bar = exp(logL2Bar)
 mnvonb = c(logLinfBar,logKBar,logL2Bar)
 
 Linf = exp(logLinf)
 K = exp(logK)
 L2 = exp(logL2)
 vbsd = exp(logvbsd)
 nobs = length(obsL)
 nponds = length(logL2)
 
 predL = sd = obsL*0; #create predL and sd same length as obsL and initialize with zeros

 for(i in 1:nobs){
   #predL = l2 +pred increase from l2 (vonB model starting at l2 at age=2)
   predL[i]=L2[pid[i]]+(Linf[pid[i]]-L2[pid[i]])*(1-exp(-K[pid[i]]*(age[i]-2)))
   #log CV that determines SD for obs length linear function of predicted length
   sd[i]=exp(logcvInt+logcvSlp*predL[i])*predL[i];
 }

 nll = -sum(dnorm(x=obsL,mean=predL,sd=sd,log=TRUE))
 
 cm = us$corr(theta)
  #calculate var-cov matrix
 Sigma = diag(vbsd) %*% cm %*% diag(vbsd)
 vonbpars = cbind(logLinf, logK, logL2)
 nprand = -sum(dmvnorm(x = vonbpars, 
                       mu=mnvonb,Sigma=Sigma,log=TRUE))
 nll+nprand;
 }
 
 re = c("logLinf", "logK", "logL2")
 obj = RTMB::MakeADFun(f, pars,random=re)
 obj$report()

 ## Call function minimizer
 fit <- nlminb(obj$par, obj$fn, obj$gr)
 fit
 sdr <- sdreport(obj)
 sdr
 summary(sdr)



 


 
 

#Data generating  and estimating for age comp data
# intended to demo flexibility of rtmb - not fully tested

library(RTMB);

### code used to generate data
set.seed(123456)
logrbar=1 #(arbitary)
Z=.15;
sd=0.1;
sel=c(.1,.15,.2,.3,.5,.68,.8,.85,.9,.94, .97, .99, .995, 1)
sampsize=c(100,100,10,50,50,50,10,200,200);
N=rep(0,14)
for(i in 1:14){
N[i]=exp(rnorm(n=1,mean=logrbar-Z*(i-1),sd=sd)) #process error in Z! 
}
prop=N*sel/sum(N*sel)

sampdat=matrix(nrow=14,ncol=9)
for (i in 1:9){
sampdat[,i]=rmultinom(1,size=sampsize[i],prob=prop)
}

##### End data generation simulation code

realdat=list(sampdat=sampdat,sampsize=sampsize,sel=sel)
pars=list(Z=0.2)
dat=realdat #just so I can tweak dat but get back to what I simulated

f=function(pars){
  getAll(dat,pars)
  p=exp(-Z*(0:13))*sel #Survival of recruit
  p=p/sum(p) #normalized to predict proportion at age
  NLL=0
  for(i in 1:9){
    NLL=NLL-dmultinom(sampdat[,i],prob=p,log=TRUE)
  }
  NLL
}

obj = RTMB::MakeADFun(f,pars)
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt
sdrep = RTMB::sdreport(obj)
summary(sdrep)

#Assume sel=1 for last four ages, monitonically increasing
logit=function(x){
  log(x/(1-x)) #x is on the real number line
}

invlogit=function(x){ 
  1/(1+exp(-x)) # x is in (0,1)
}

adjust=rep(0.9,10)
tradjust=logit(adjust)
pars=list(Z=Z,tradjust=tradjust)

f=function(pars){
  getAll(dat,pars)
  adj=invlogit(tradjust)
  seluse=rep(1,14)
  for(i in 10:1){seluse[i]=seluse[i+1]*adj[i]}
  ptmp=exp(-Z*(0:13))*seluse; #Sel adjusted Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  NLL=0
  for(i in 1:9){
    NLL=NLL-dmultinom(sampdat[,i],prob=p,log=TRUE)
    REPORT(seluse)
  }
  NLL
}

obj = RTMB::MakeADFun(f,pars)
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt$convergence
sdrep <- RTMB::sdreport(obj)
sdrep
obj$report()
opt$convergence

#logistic selectivity
logistic=function(x,inf,slp){
1/(1+exp(-slp*(x-inf)))
}
logistic(1:14,5,0.5)

loglogistinf=log(5)
loglogistslp=log(0.5)
pars=list(Z=Z,loglogistinf=loglogistinf,loglogistslp
             =loglogistslp)

f=function(pars){
  getAll(dat,pars)
  logistslp=exp(loglogistslp)
  logistinf=exp(loglogistinf)
  seluse=logistic(1:14,logistinf,logistslp)
  ptmp=exp(-Z*(0:13))*seluse; #Sel adj Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  NLL=0
  for(i in 1:9){
    NLL=NLL-dmultinom(sampdat[,i],prob=p,log=TRUE)
    REPORT(seluse)
  }
  NLL;
}

obj = RTMB::MakeADFun(f,pars)
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt$convergence
sdrep <- RTMB::sdreport(obj)
sdrep
obj$report()

# Now for fun, assume a linear form of Dirichlet-multinomial dsn
# Check out e.g., Fisch et al CJFAS 79:1745–1764 (2022)for
# log likelihood

# Dirichlet-multinomial is a compound distribution which assumes 
# probabilities come from a Dirichlet and given those probs, 
# n from multinomial.
  
ddirmult=function(x,alpha,log=FALSE){
  n=sum(x);
  salpha=sum(alpha);
  lden=lgamma(salpha)+lgamma(n+1)-lgamma(n+salpha) +
    sum(lgamma(x+alpha))-sum(lgamma(alpha))-sum(lgamma(x+1))
  if(log) ret=lden else ret=exp(lden)
  ret
}

# apply my dirichlet multinomial density to see how it works
 Nobs=dat$sampdat[,1];
 ddirmult(Nobs,Nobs*.8+0.01,log=T)
 
 pars=list(Z=Z,trscale=logit(0.8))
 
 f=function(pars){
   getAll(dat,pars)
   ptmp=exp(-Z*(0:13))*sel #Sel adj Survival of recruit
   p=ptmp/sum(ptmp) #normalized to predict proportion at age
   scale=invlogit(trscale)
   NLL=0;
   for(i in 1:9){
     sampsize=sum(sampdat[,i])
     alpha=sampsize*scale*p
     NLL=NLL-ddirmult(sampdat[,i],alpha=alpha,log=TRUE)
   }
   NLL
 }
 
 obj = RTMB::MakeADFun(f,pars)
 opt <- nlminb(obj$par,obj$fn,obj$gr)
 opt$convergence
 sdrep <- RTMB::sdreport(obj)
 sdrep

  #Make the sample sizes way big for info content
 dat=realdat
 dat$sampdat=dat$sampdat*10
 dat$sampsize=dat$sampsize*10
 f=function(pars){
   getAll(dat,pars)
   ptmp=exp(-Z*(0:13))*sel #Sel adj Survival of recruit
   p=ptmp/sum(ptmp); #normalized to predict proportion at age
   scale=invlogit(trscale)
   NLL=0
   for(i in 1:9){
     sampsize=sum(sampdat[,i])
     alpha=sampsize*scale*p
     NLL=NLL-ddirmult(sampdat[,i],alpha=alpha,log=TRUE)
   }
   NLL
 }
 
 obj = RTMB::MakeADFun(f,pars)
 opt <- nlminb(obj$par,obj$fn,obj$gr)
 opt$convergence
 sdrep <- RTMB::sdreport(obj)
 sdrep
 
 #In case someone wants to simulate
 # would have to do this manually, i.e., not using obj$simulate()
 rdirmult=function(n,size,alpha){
   #not vectorized!
   res=list();
   k=length(alpha);
   for (i in 1:n){
     y=stats::rgamma(k,alpha,1);
     x=y/sum(y); #x~dirichlet
     res[[i]]=as.vector(stats::rmultinom(1,size=size,prob=x));
   }
   do.call(rbind,res);
 }
 rdirmult(1,sum(Nobs),alpha=Nobs*.5);
 
#"Manual simulation
# Use bloolean to determine if f generates sim data or calcs NLL
f=function(pars){
  getAll(dat,pars)
  if (simcode) simlst=list()
  # note that sel is assumed to exist (calc outside function) not estimated
  ptmp=exp(-Z*(0:13))*sel #Sel adj Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  scale=invlogit(trscale)
  NLL=0
  for(i in 1:9){
    sampsize=sum(sampdat[,i])
    alpha=sampsize*scale*p
    NLL=NLL-ddirmult(sampdat[,i],alpha=alpha,log=TRUE)
    if(simcode) simlst[[i]]=as.vector(rdirmult(1,sampsize,alpha))
  }
  if(simcode) {
    simdat=do.call(cbind,simlst)
    ret=simdat
  }
  else ret=NLL
  ret;
}

#Showing this still works to estimate when simcode is FALSE
simcode=FALSE
dat=realdat
dat$sampdat=dat$sampdat*10
dat$sampsize=dat$sampsize*10
obj = RTMB::MakeADFun(f,pars)
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt$convergence
sdrep <- RTMB::sdreport(obj)
sdrep;
invlogit(opt$par[2]) 

# Now using f to simulate

#get ests as list
ests=as.list(sdrep,"Est")
attr(ests,"what")=NULL
ests$trscale=-2.2
simcode=TRUE
simsampdat=f(ests)
simsampdat

#Fit model to simulated data
dat=realdat
dat$sampdat=simsampdat
dat$sampsize=colSums(simsampdat)
simcode=FALSE
simobj = RTMB::MakeADFun(f,pars)
simopt <- nlminb(simobj$par,simobj$fn,simobj$gr)
simopt$convergence
simsdrep <- RTMB::sdreport(simobj)
simsdrep
invlogit(simopt$par[2]) 


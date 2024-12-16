# R script for fitting musky L@Age with RTMB
# profile CI and simulation add ons
# Adapted from admb and TMB classes by Bence and Brenden
set.seed(123456)
library(RTMB)

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T)

#Set up the data and starting value of parameters for RTMB
gmdat = list(lenobs=gmRdat[,"Length"],age=gmRdat[,"Age"])
par = list(logLinf=7,logK=-1.6,t0=0,logsd=4)


f = function(par){
  getAll(dat,par)
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

#TIP - assign data you want to use to list used in your NLL function
#   immediately before MakeADFun (when using different datasets)
#   or use closure function (covered later)
#Create object and fit model to original data
dat=gmdat
obj=RTMB::MakeADFun(f,par)
fit=nlminb(obj$par, obj$fn, obj$gr)
fit$convergence
sdr=sdreport(obj)
sdr

simdat=gmdat;  #copy real data
#Simulate obs lengths, write to data copy,refit model
simdat$lenobs=obj$simulate()$lenobs
dat=simdat
objsim <- RTMB::MakeADFun(f,par)
simfit = nlminb(objsim$par, objsim$fn, objsim$gr)
simconv=simfit$convergence
simsdr=sdreport(objsim)

# Pulling out results we want
simest=as.list(simsdr, "Est")
attr(simest,"what")=NULL
simse=as.list(simsdr,"Std")
attr(simse,"what")=NULL

#Lists and do.call very useful for sims
#List sizes do not need to be specified in advance
#Writing to lists is numerically more efficient
lstest = list()
lstse = list()
lstconv = list()

#Here just to illustrate we write our set of simulation
# results to the first element of list then rerun our sim lines
# and write the second set to the second element;
lstest[[1]]=unlist(simest)
lstse[[1]]=unlist(simse)
lstconv[[1]]=unlist(simconv)

# REMEMBER to rerun the simulation before writing these lines
# In a full simulation you would run your simulation in a loop
# and write to the "ith" element each time
lstest[[2]]=unlist(simest)
lstse[[2]]=unlist(simse)
lstconv[[2]]=unlist(simconv)

# After running all the iterations we can use do.call
#  to put results of a particular type together in matrix
# e.g., each row a simulation run, each col the est for a different param

#First put all results in one list (not essential);
reslst<-list(est=lstest,se=lstse,conv=lstconv)

#Then create the matrix of estimates
matrixest=do.call(rbind,reslst$est)
matrixSEs=do.call(rbind,reslst$se)
vecconv=as.vector(do.call(rbind,reslst$conv))

#Full simulation code

# Say you wanted to simulate with starting pars (or others)
# rather than estimated
#  Recreate obj using desired pars and don't fit!
#datlst=gmdat; 
#obj=RTMB::MakeADFun(f,otherpar);


set.seed(54321)
datlst=gmdat  #copy real data
lstest = list()
lstse = list()
lstconv = list()


for(i in 1:1000){
#Simulate obs lengths, write to data copy,refit model
dat$lenobs=obj$simulate()$lenobs
objsim <- RTMB::MakeADFun(f,par,silent=TRUE)
simfit = nlminb(objsim$par, objsim$fn, objsim$gr)
simconv=simfit$convergence
simsdr=sdreport(objsim)

# Pulling out results we want
simest=as.list(simsdr, "Est")
attr(simest,"what")=NULL
simse=as.list(simsdr,"Std")
attr(simse,"what")=NULL

lstest[[i]]=unlist(simest)
lstse[[i]]=unlist(simse)
lstconv[[i]]=unlist(simconv)
}

#First put all results in one list (not essential)
reslst<-list(est=lstest,se=lstse,conv=lstconv)

#Then create the matrix of estimates
mest=do.call(rbind,reslst$est)
mSEs=do.call(rbind,reslst$se)
mconv=do.call(rbind,reslst$conv)

#Plot results and calculate 95% par bootstrap CI
hist(mest[,"logLinf"])
quantile(mest[,"logLinf"],probs=c(0.025,0.975))

#coverage
#note if you are simulating from other than estimated pars
# this needs minor changes.

ci=mest[,"logLinf"]+cbind(-mSEs[,"logLinf"],mSEs[,"logLinf"])*qnorm(0.975)
ci=as.data.frame(ci)
names(ci)=c("lb","ub")
head(ci)
tstci=with(ci, lb <= fit$par["logLinf"] & ub >= fit$par["logLinf"])
#count up the cases its in the CI (if procedure right should be ~950)
sum(tstci)

#Illustrate how to create a list of estimates, and change one
#Pull out parameter estimates and organize as parlst
estparlst=as.list(sdr,"Est");
attr(estparlst,"what")=NULL;
estparlst #same form as the starting value list
#change one parameter
estparlst$loglinf=7
estparlst

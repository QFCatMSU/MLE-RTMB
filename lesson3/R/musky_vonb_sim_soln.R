# R script for fitting musky L@Age with RTMB
# profile CI and simulation add ons
# Adapted from admb and TMB classes by Bence and Brenden
set.seed(123456)
library(RTMB);

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T);

#Set up the data and starting value of parameters for RTMB
gmdat_real = list(len_obs=gmRdat[,"Length"],age=gmRdat[,"Age"]);
gmpar = list(log_linf=7,log_vbk=-1.6,t0=0,log_sd=4);

NLL_fun = function(par_lst){
  getAll(gmdat,par_lst);
  linf = exp(log_linf);
  vbk = exp(log_vbk);
  sd = exp(log_sd);
  len_obs=OBS(len_obs); # New line!
  len_pred = linf * (1 - exp(-vbk * (age - t0)));
  nll = -sum(dnorm(len_obs, len_pred, sd, TRUE));
  atage_pred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atage_pred);
  nll
}

#TIP - assign data you want to use to list used in your NLL function
#   immediately before MakeADFun (when using different datasets)
#Create object and fit model to original data
gmdat=gmdat_real; 
obj=RTMB::MakeADFun(NLL_fun,gmpar);
fit=nlminb(obj$par, obj$fn, obj$gr);
fit$convergence;
sdr=sdreport(obj);
sdr;

simdat=gmdat_real;  #copy real data
#Simulate obs lengths, write to data copy,refit model
simdat$len_obs=obj$simulate()$len_obs;
gmdat=simdat;
objsim <- RTMB::MakeADFun(NLL_fun,gmpar);
simfit = nlminb(objsim$par, objsim$fn, objsim$gr);
sim_conv=simfit$convergence;
sim_sdr=sdreport(objsim);

# Pulling out results we want
sim_est=as.list(sim_sdr, "Est");
attr(sim_est,"what")=NULL;
sim_se=as.list(sim_sdr,"Std");
attr(sim_se,"what")=NULL;

#Lists and do.call very useful for sims
#List sizes do not need to be specified in advance
#Writing to lists is numerically more efficient
lst_est = list();
lst_se = list();
lst_conv = list();

#Here just to illustrate we write our set of simulation
# results to the first element of list then rerun our sim lines
# and write the second set to the second element;
lst_est[[1]]=unlist(sim_est);
lst_se[[1]]=unlist(sim_se);
lst_conv[[1]]=unlist(sim_conv);

# REMEMBER to rerun the simulation before writing these lines
# In a full simulation you would run your simulation in a loop
# and write to the "ith" element each time
lst_est[[2]]=unlist(sim_est);
lst_se[[2]]=unlist(sim_se);
lst_conv[[2]]=unlist(sim_conv);

# After running all the iterations we can use do.call
#  to put results of a particular type together in matrix
# e.g., each row a simulation run, each col the est for a different param

#First put all results in one list (not essential);
res_lst<-list(est=lst_est,se=lst_se,conv=lst_conv);

#Then create the matrix of estimates
matrix_est=do.call(rbind,res_lst$est);
matrix_SEs=do.call(rbind,res_lst$se);
vec_conv=as.vector(do.call(rbind,res_lst$conv));

#Full simulation code

set.seed(54321);
gmdat=gmdat_real;  #copy real data
lst_est = list();
lst_se = list();
lst_conv = list();

for(i in 1:1000){
#Simulate obs lengths, write to data copy,refit model
gmdat$len_obs=obj$simulate()$len_obs;
objsim <- RTMB::MakeADFun(NLL_fun,gmpar,silent=TRUE);
simfit = nlminb(objsim$par, objsim$fn, objsim$gr);
sim_conv=simfit$convergence;
sim_sdr=sdreport(objsim);

# Pulling out results we want
sim_est=as.list(sim_sdr, "Est");
attr(sim_est,"what")=NULL;
sim_se=as.list(sim_sdr,"Std");
attr(sim_se,"what")=NULL;

lst_est[[i]]=unlist(sim_est);
lst_se[[i]]=unlist(sim_se);
lst_conv[[i]]=unlist(sim_conv);
}

#First put all results in one list (not essential);
res_lst<-list(est=lst_est,se=lst_se,conv=lst_conv);

#Then create the matrix of estimates
m_est=do.call(rbind,res_lst$est);
m_SEs=do.call(rbind,res_lst$se);
m_conv=do.call(rbind,res_lst$conv);

#Plot results and calculate 95% par bootstrap CI
hist(m_est[,"log_linf"])
quantile(m_est[,"log_linf"],probs=c(0.025,0.975));

#coverage

ci=m_est[,"log_linf"]+cbind(-m_SEs[,"log_linf"],m_SEs[,"log_linf"])*qnorm(0.975);
ci=as.data.frame(ci)
names(ci)=c("lb","ub");
head(ci);
tst_ci=with(ci, lb<= fit$par["log_linf"] & ub >= fit$par["log_linf"])
#count up the cases its in the CI (if procedure right should be ~950)
sum(tst_ci);

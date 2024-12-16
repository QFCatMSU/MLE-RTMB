# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
# Here we demo fitting by REML or ML and compare
# Show that you get same estimates except sd becomes approx unbiased estimator

#Also added for 2024 (required adding observation indicator OBS to function)
#* -demo of checkConsistency
#* -demo of OneStepPredict

library(RTMB)
gmdat = read.table("lesson2/data/musky_vonb.dat", head = T)
# Set up the data and starting value of parameters for RTMB
dat = list(lenobs = gmdat[, "Length"], age = gmdat[, "Age"])
par = list(logLinf = 7, logK = -1.6, t0 = 0, logsd = 4)

# My NLL function
f = function(par) {
  getAll(dat, par)
  lenobs = OBS(lenobs)
  Linf = exp(logLinf)
  K = exp(logK)
  sd = exp(logsd)
  lenpred = Linf * (1 - exp(-K * (age - t0)))
  nll = -sum(dnorm(lenobs, lenpred, sd, TRUE))
  atagepred = Linf * (1 - exp(-K * ((1:11) - t0)))
  REPORT(atagepred)
  nll
}

# Create two versions of "obj" one set up for ML estimation and one for REML
obj = MakeADFun(f, par)
re = c("logLinf", "logK", "t0")
objreml = MakeADFun(f, par, random = re)

# Call function minimizer
fit = nlminb(obj$par, obj$fn, obj$gr)
fitreml = nlminb(objreml$par, objreml$fn, objreml$gr)

# Get parameter uncertainties and convergence diagnostics
sdr = sdreport(obj)
summary(sdr)
sdrreml = sdreport(objreml)
summary(sdrreml)

# For reml logsd is now listed first in summary output
# - From RTMB perspective it is the only fixed effect par!
# You won't see the other params if you just print sdrreml
# sd is larger with reml (adjusting for negative bias of ML est)
# Notice that the other "parameters" have nearly identical values

# Numerical comparison of sds
sd = exp(fit$par[4]) # ml sd
sdreml = exp(fitreml$par[1]) # reml sd
n = length(dat$lenobs)
# exclude sd (want pars that determine E(L$age))
p = length(obj$par) - 1
# By known theory this should be same as REML
sdadj = sqrt((sd^2) * (n / (n - p)))
#print sd ests to compare
sd
sdreml
sdadj

#Recursive quantile (osa) residuals
osa = oneStepPredict(obj)
#osa <- oneStepPredict(obj, method="fullGaussian")
osa <- oneStepPredict(obj, method="oneStepGeneric") #most accurate and slow
hist(osa$residual)
plot(osa$reesidual~dat$age)
qqnorm(osa$res); abline(0,1)

#Check consistency of Laplace approximation
#tests that expected joint score function =0
#Also estimates parameter bias caused by the Laplace approximation
set.seed(111)
chk = checkConsistency(obj,estimate=TRUE)
chk

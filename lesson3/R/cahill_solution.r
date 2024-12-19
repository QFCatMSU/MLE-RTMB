# R script for fitting musky L@Age with RTMB
# profile CI and simulation add ons
# Adapted from admb and TMB classes by Bence and Brenden
set.seed(123456)
library(RTMB)

gmRdat <- read.table("lesson3/data/musky_vonb.dat", head = T)

# Set up the data and starting value of parameters for RTMB
dat <- list(lenobs = gmRdat[, "Length"], age = gmRdat[, "Age"])
par <- list(logLinf = 7, logK = -1.6, t0 = 0, logsd = 4)

f <- function(par) {
    getAll(dat, par)
    Linf <- exp(logLinf)
    K <- exp(logK)
    sd <- exp(logsd)
    lenobs <- OBS(lenobs) # New line!
    lenpred <- Linf * (1 - exp(-K * (age - t0)))
    nll <- -sum(dnorm(lenobs, lenpred, sd, TRUE))
    atagepred <- Linf * (1 - exp(-K * ((1:11) - t0)))
    REPORT(atagepred)
    nll
}

obj <- MakeADFun(f, par)
opt <- nlminb(obj$par, obj$fn, obj$gr) # estimate

# sim study full
odat <- dat
doone <- function() {
    dat <<- list()
    dat$lenobs <<- obj$simulate()$lenobs
    objsim <- MakeADFun(f, par, silent = TRUE)
    fitsim <- nlminb(objsim$par, objsim$fn, objsim$gr)
    fitsim$par
}

set.seed(1)
sim <- replicate(100, doone())
dat <- odat

# plot it
boxplot(t(sim))
points(1:length(opt$par), opt$par,
    cex = 5,
    pch = 4, lwd = 3, col = "darkgreen"
)

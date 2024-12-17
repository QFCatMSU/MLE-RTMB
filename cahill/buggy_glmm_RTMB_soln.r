# see thorson and kristensen 2024 chp 2

library(RTMB)

data <- list(y_i = c(
    0, 3, 4, 0, 12, 57, 14, 4,
    14, 101, 1, 76, 39, 19, 2, 1
))

par <- list(
    "ln_mu" = 0, "ln_sd" = 0,
    "eps_i" = rep(0, length(data$y_i) + 1)
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    y_i <- OBS(y_i)
    jnll <- 0
    for (i in 1:length(y_i)) {
        jnll <- jnll - dnorm(eps_i[i], ln_mu, exp(ln_sd), TRUE)
    }
    yhat_i <- exp(eps_i)
    for (i in 1:length(y_i)) {
        jnll <- jnll - dpois(y_i[i], yhat_i[i], TRUE)
    }
    yhat_sum <- sum(yhat_i)
    REPORT(jnll)
    REPORT(yhat_i)
    REPORT(yhat_sum)
    ADREPORT(yhat_i)
    ADREPORT(yhat_sum)
    jnll
}

obj <- MakeADFun(f, par, random = c("eps_i"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

#---------------------------------------------------------------------
# tracking down bugs in the inner optimization
# pluck out the inner hessian

inner_hessian <- obj$env$spHess(random = TRUE)
Matrix::image(inner_hessian,
    main = "Inner Hessian in TMB\nwith non-invertible inner Hessian"
)

# map off fixed effects and treat random as fixed
FE <- setdiff(names(par), "eps_i")
map <- lapply(FE, function(x) factor(par[[x]] * NA))
names(map) <- FE

# rebuild and run optimizer holding outer (fixed effects) constant via map
test <- MakeADFun(f, par, map = map)
opt_random <- nlminb(test$par, test$fn, test$gr)

# optimize inner Hessian
inner_hessian <- optimHess(
    par = opt_random$par,
    fn = test$fn, gr = test$gr
)

# check eigendecomposition, i.e.,
# searching for local maxima or flat areas in solution space
Eigen <- eigen(inner_hessian)
bad_vectors <- Eigen$vectors[, which(Eigen$values <= 0)]

data.frame(
    "parameter" = names(test$par), "MLE" = test$par,
    "curvature" = ifelse(bad_vectors > 0.1, ":(", "Good")
)

# locate the misery and deal with it

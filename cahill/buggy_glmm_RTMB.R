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

# Define the function
f <- function(x) 2 + 0.5 * x + sin(x)
# ----------------------------------------------------------
# Symbolic differentiation
library(Deriv)
f_prime <- Deriv(f, "x") # First derivative
f_double_prime <- Deriv(f_prime, "x") # Second derivative
f_prime(1)
f_double_prime(1)
#-----------------------------------------------------------
# Automatic differentiation
library(RTMB)
F <- MakeTape(f, numeric(1)) # one valued input
f_prime <- F$jacfun()
f_prime(1)
f_double_prime <- f_prime$jacfun()
f_double_prime(1)

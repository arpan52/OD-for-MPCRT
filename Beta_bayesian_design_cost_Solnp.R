rm(list = ls())  # Clear all variables

library(Rsolnp)
library(pracma)  # For integral2 function

# Parameters
T_star <- 120
c_values <- c(rep(1,3),rep(1.1,3),rep(1.2,3),rep(0.7,3))

m <- length(c_values)  # Number of groups

sigma <- rep(1.8, m)



L <- 7  # Lower bound
U <- 25  # Upper bound

# Beta Bayesian parameters
y1 <- c(rep(1,3), rep(1,3), rep(1,3), rep(1,3))
z1 <- c(rep(999,3), rep(99,3), rep(49,3), rep(39,3))
y2 <- c(rep(18,3), rep(3,3), rep(17,3), rep(22,3))
z2 <- c(rep(7,3), rep(1,3), rep(3,3), rep(3,3))

# Feasible initial values satisfying sum(c * n) = T_star
x0 <- rep(T_star / sum(c_values), m)
x0 <- pmax(L, pmin(x0, U))  # Ensure within bounds

# Objective function: variance of beta design
var_beta_design <- function(n) {
  Int <- numeric(m)
  eff_beta1 <- numeric(m)
  
  for (j in 1:m) { 
    fun <- function(r1, r2) {
      ((r1^(y1[j] - 1)) * (r2^(y2[j] - 1)) * ((1 - r1)^(z1[j] - 1)) * ((1 - r2)^(z2[j] - 1))) /
        ((sigma[j]^2) * (1 + ((n[j] - 1) * r1 - r2)))
    }
    
    Int[j] <- integral2(fun, 0, 1, 0, 1, reltol = 1e-10)$Q
    eff_beta1[j] <- (n[j] / (beta(y1[j], z1[j]) * beta(y2[j], z2[j]))) * Int[j]
  }
  
  eff_beta <- -(1 / 2) * sum(eff_beta1)
  return(eff_beta)  # Return a numeric value
}

# Constraint function: total budget restriction (must return zero when satisfied)
const1 <- function(n) {
  return(sum(c_values * n) - T_star)
}

# Run optimization using solnp()
res <- solnp(
  pars = x0,               # Initial values
  fun = var_beta_design,   # Objective function
  eqfun = const1,          # Equality constraint function
  eqB = 0,                 # Ensures sum(c * n) = T_star
  LB = rep(L, m),          # Lower bounds
  UB = rep(U, m)           # Upper bounds
)

n_opt <- round(res$pars)  # Rounded final values

# Print results and check constraint satisfaction
print(n_opt)
cat("Check Constraint:", sum(c_values * n_opt), "vs", T_star, "\n")

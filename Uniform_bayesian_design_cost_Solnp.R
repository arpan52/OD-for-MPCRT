library(Rsolnp)

#set.seed(123)  # Ensures reproducibility

# Define input parameters
rho<- c(rep(0.001,3),rep(0.01,3),rep(0.02,3),rep(0.025,3))
rho1<- c(rep(0.72,3),rep(0.75,3),rep(0.85,3),rep(0.88,3))
sigma<- rep(1.8,12)
c<- c(rep(1,3),rep(1.1,3),rep(1.2,3),rep(0.7,3))
T_star <- 120  # Example value
m<-length(c)

# Define parameter bounds
a1 <- rep(0,m)
b1 <- rho + 0.02
a2 <- rho1 - 0.1
b2 <- rho1 + 0.1

lb <- rep(7, m)  # Lower bounds
ub <- rep(25, m)  # Upper bounds

# Compute a feasible initial value that satisfies sum(c * n) = T_star
x0 <- rep(T_star / sum(c), m)  # Equal allocation

# Ensure x0 is within bounds
x0 <- pmax(lb, pmin(x0, ub))

# Define the function var_design (objective function)
var_design <- function(n) {
  var <- numeric(m)
  
  for (j in 1:m) {
    k <- n[j] / (sigma[j]^2 * (b1[j] - a1[j]) * (b2[j] - a2[j]) * (n[j] - 1))  # Avoid div by zero
    k1 <- 1 - a2[j] + (b1[j] * (n[j] - 1))
    k2 <- 1 - b2[j] + (b1[j] * (n[j] - 1))
    k3 <- 1 - a2[j] + (a1[j] * (n[j] - 1))
    k4 <- 1 - b2[j] + (a1[j] * (n[j] - 1))
    var[j] <- k * ((k1 * (log(k1) - 1)) - (k2 * (log(k2) - 1)) -
                     (k3 * (log(k3) - 1)) + (k4 * (log(k4) - 1)))
    
  }
  
  return(-0.5 * sum(var))  # Return negative because we minimize
}

# Define equality constraint function (must return zero when satisfied)
eq_constraint <- function(n) {
  return(sum(c * n) - T_star)  # Must be zero
}

# Solve using `solnp()` from the `Rsolnp` package
res <- solnp(
  pars = x0,              # Initial values
  fun = var_design,       # Objective function
  eqfun = eq_constraint,  # Equality constraint
  eqB = 0,                # Right-hand side (must be exactly zero)
  LB = lb,                # Lower bounds
  UB = ub                 # Upper bounds
)

n_uni <- round(res$pars)  # Rounded final values

# Print the result and check constraint satisfaction
print(n_uni)
cat("Check Constraint:", sum(c * n_uni), "vs", T_star, "\n")

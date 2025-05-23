library(MASS)  # For unifrnd equivalent
#library(ggplot2)  # For boxplot visualization

al<-0.05
delta <- 0.4  # Placeholder for effect size
rho<- c(rep(0.001,3),rep(0.01,3),rep(0.02,3),rep(0.025,3))
rho1<- c(rep(0.72,3),rep(0.75,3),rep(0.85,3),rep(0.88,3))

# Initialize parameters
n_opt<- round(c(19,19,19,7,7,7,7,7,7,7,7,7))
n_uni <- round(c(13,13,13,8,8,8,7,7,7,14,14,14))
n_beta <- round(c(7,7,7,13,13,13,7,7,7,14,14,14))
n_minmax <- round(c(15,15,15,10,10,10,7,7,7,9,9,9))

m <- length(n_opt)
sigma <- rep(1.8, m)  # Assuming sigma is defined elsewhere

t1 = qt(al/2,m-1);
t2 = qt(1-(al/2),m-1);

# Define parameter bounds
a1 <- rep(0,m)
b1 <- rho + 0.02
a2 <-  rho1 - 0.1
b2 <- rho1 + 0.1

# Number of simulations
n_sim <- 1e4

# Initialize power vectors
pwr_opt <- numeric(n_sim)
pwr_uni <- numeric(n_sim)
pwr_beta <- numeric(n_sim)
pwr_minmax <- numeric(n_sim)
V_opt <- numeric(m)
V_uni <- numeric(m)
V_beta <- numeric(m)
V_minmax <- numeric(m)


# Loop over simulations
for (i in 1:n_sim) {
  
  r <- mapply(runif, n = 1, min = a1, max = b1) 
  
  r1 <- mapply(runif, n = 1, min = a2, max = b2) 
  

  for (j in 1:m) {
    V_opt[j] <- n_opt[j] / ((sigma[j]^2) * (1 + ((n_opt[j] - 1) * r[j] - r1[j])))
  }
    V_opt_sum<- 2/sum(V_opt)
    
 for (j in 1:m) {
      V_uni[j] <- n_uni[j] / ((sigma[j]^2) * (1 + ((n_uni[j] - 1) * r[j] - r1[j])))
    }
    V_uni_sum<- 2/sum(V_uni)
    
  for (j in 1:m) {
      V_beta[j] <- n_beta[j] / ((sigma[j]^2) * (1 + ((n_beta[j] - 1) * r[j] - r1[j])))
    }
    V_beta_sum<- 2/sum(V_beta)
    
    for (j in 1:m) {
      V_minmax[j] <- n_minmax[j] / ((sigma[j]^2) * (1 + ((n_minmax[j] - 1) * r[j] - r1[j])))
    }
    V_minmax_sum<- 2/sum(V_minmax)
    
   pwr_opt[i] <- 1 - pt(t2, df = m - 1, ncp = delta / sqrt(V_opt_sum)) +
    pt(t1, df = m - 1, ncp = delta / sqrt(V_opt_sum))
  
  pwr_uni[i] <- 1 - pt(t2, df = m - 1, ncp = delta / sqrt(V_uni_sum)) +
    pt(t1, df = m - 1, ncp = delta / sqrt(V_uni_sum))
  
  pwr_beta[i] <- 1 - pt(t2, df = m - 1, ncp = delta / sqrt(V_beta_sum)) +
    pt(t1, df = m - 1, ncp = delta / sqrt(V_beta_sum))
  
  pwr_minmax[i] <- 1 - pt(t2, df = m - 1, ncp = delta / sqrt(V_minmax_sum)) +
    pt(t1, df = m - 1, ncp = delta / sqrt(V_minmax_sum))
}

# Boxplot visualization

boxplot(pwr_opt, pwr_uni, pwr_beta, pwr_minmax,
        names = c("LOD", "Uniform", "Beta", "Min-max"),
        col = c("lightblue", "pink", "lightgreen", "lightcoral"),
        ylab = "Power",
        xlab = "Design",
        border = "darkblue"
       )

# Load necessary library
rm(list = ls())
library(MASS)  # For mvrnorm

# Define parameters
r <- c(rep(0.014,3),rep(0.014,3),rep(0.014,3),rep(0.014,3))#rep(0.014,12)
sig <- rep(1.8, 12)  # Standard deviation
a <- c(rep(0.17, 3), rep(0.17, 3), rep(0.17, 3), rep(0.17, 3))  # Lower bounds
b <- c(rep(0.23, 3), rep(0.23, 3), rep(0.23, 3), rep(0.23, 3))  # Upper bounds

n_des <- c(25, 25, 25, 25, 25, 25, 7, 7, 7, 7, 7, 7)  # Sizes for the optimal design
n_bal <- rep(10, 12)  # Sizes for the balanced design
alpha <- 0.05  # Significance level

n_samples <- 1000  # Number of samples
n_sim <- 1000 # Number of simulations

Pwr_opt <- rep(0, n_sim)  # Power for optimal design
Pwr_bal <- rep(0, n_sim)  # Power for balanced design

mu <- 5.6  # Mean for control group
beta <- 0.4  # Effect size

for (v in 1:n_sim) {
  
  # Simulate r1 for this specific simulation v (same r1 for both designs)
  r1 <- mapply(runif, n = 1, min = a, max = b)
  
  # ------------------- Simulations for Optimal Design -------------------
  v_des <- rep(0, length(n_des))
  for (k in 1:length(n_des)) {
    v_des[k] <- (n_des[k]) / ((sig[k]^2) * (1 + ((n_des[k] - 1) * (r[k]) - r1[k])))
  }
  var_n_des <- 2 / sum(v_des)
  
  # Simulate the optimal design samples
  mean_vec <- cbind(t(rep(mu, sum(n_des))), t(rep((mu + beta), sum(n_des))))
  correlation_matrix <- matrix(0, nrow = 2 * sum(n_des), ncol = 2 * sum(n_des))
  s <- 1
  u <- sum(n_des) + 1
  
  for (i in 1:length(n_des)) {
    matrix_ICC <- diag((1 - r[i]), n_des[i], n_des[i]) + matrix(r[i], nrow = n_des[i], ncol = n_des[i])
    matrix_match <- diag(r1[i], n_des[i], n_des[i])
    
    correlation_matrix[s:(s + n_des[i] - 1), s:(s + n_des[i] - 1)] <- matrix_ICC
    correlation_matrix[s:(s + n_des[i] - 1), u:(u + n_des[i] - 1)] <- matrix_match
    correlation_matrix[u:(u + n_des[i] - 1), s:(s + n_des[i] - 1)] <- matrix_match
    correlation_matrix[u:(u + n_des[i] - 1), u:(u + n_des[i] - 1)] <- matrix_ICC
    
    s <- s + n_des[i]
    u <- u + n_des[i]
  }
  
  samples <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = correlation_matrix)
  
  # Calculate power for optimal design
  ctr_opt <- 0
  for (j in 1:n_samples) {
    s <- 1  # Reset s for each sample
    beta_hat <- rep(0, n_samples)
    for (i in 1:length(n_des)) {
      Vj <- diag((1 - r[i]), n_des[i], n_des[i]) + matrix(r[i], nrow = n_des[i], ncol = n_des[i])
      Dj <- diag(r1[i], n_des[i], n_des[i])
      V_inv <- solve(Vj)
      Cj <- (1 / (sig[i])^2) * solve(Vj - ((r1[i])^2 * V_inv))
      Ej <- -(1 / (sig[i])^2) * V_inv %*% Dj %*% Cj
      
      # Calculate b1, b2, and b3 for each group
      b1 <- (1 / (sig[i])^2) * t(rep(1, n_des[i])) %*% (Cj - Ej) %*% (samples[j, (sum(n_des) + s):(sum(n_des) + s + n_des[i] - 1)] - samples[j, s:(s + n_des[i] - 1)])
      b2 <- (1 / (sig[i])^2) * t(rep(1, n_des[i])) %*% (Cj - Ej) %*% matrix(1, nrow = n_des[i], ncol = 1)
      
      s <- s + n_des[i]
    }
    beta_hat <- sum(b1)  / sum(b2)

    t_stat <- (beta_hat) / sqrt(var_n_des)
    p_value <- 2 * (1 - pt(abs(t_stat), df = (length(n_des)) - 1))  # Two-tailed p-value
    if (p_value < alpha) {
      ctr_opt <- ctr_opt + 1
    }
  }
  Pwr_opt[v] <- ctr_opt / n_samples
  
  # ------------------- Simulations for Balanced Design -------------------
  v_bal <- rep(0, length(n_bal))
  for (k in 1:length(n_bal)) {
    v_bal[k] <- (n_bal[k]) / ((sig[k]^2) * (1 + ((n_bal[k] - 1) * (r[k]) - r1[k])))
  }
  var_n_bal <- 2 / sum(v_bal)
  
  # Simulate the balanced design samples
  mean_vec <- cbind(t(rep(mu, sum(n_bal))), t(rep((mu + beta), sum(n_bal))))
  correlation_matrix <- matrix(0, nrow = 2 * sum(n_bal), ncol = 2 * sum(n_bal))
  s <- 1
  u <- sum(n_bal) + 1
  
  for (i in 1:length(n_bal)) {
    matrix_ICC <- diag((1 - r[i]), n_bal[i], n_bal[i]) + matrix(r[i], nrow = n_bal[i], ncol = n_bal[i])
    matrix_match <- diag(r1[i], n_bal[i], n_bal[i])
    
    correlation_matrix[s:(s + n_bal[i] - 1), s:(s + n_bal[i] - 1)] <- matrix_ICC
    correlation_matrix[s:(s + n_bal[i] - 1), u:(u + n_bal[i] - 1)] <- matrix_match
    correlation_matrix[u:(u + n_bal[i] - 1), s:(s + n_bal[i] - 1)] <- matrix_match
    correlation_matrix[u:(u + n_bal[i] - 1), u:(u + n_bal[i] - 1)] <- matrix_ICC
    
    s <- s + n_bal[i]
    u <- u + n_bal[i]
  }
  
  samples <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = correlation_matrix)
  
  # Calculate power for balanced design
  ctr_bal <- 0
  for (j in 1:n_samples) {
    s <- 1  # Reset s for each sample
    beta_hat <- rep(0, n_samples)
    for (i in 1:length(n_bal)) {
      Vj <- diag((1 - r[i]), n_bal[i], n_bal[i]) + matrix(r[i], nrow = n_bal[i], ncol = n_bal[i])
      Dj <- diag(r1[i], n_bal[i], n_bal[i])
      V_inv <- solve(Vj)
      Cj <- (1 / (sig[i])^2) * solve(Vj - ((r1[i])^2 * V_inv))
      Ej <- -(1 / (sig[i])^2) * V_inv %*% Dj %*% Cj
      
      # Calculate b1, b2, and b3 for each group
      b1 <- (1 / (sig[i])^2) * t(rep(1, n_bal[i])) %*% (Cj - Ej) %*% (samples[j, (sum(n_bal) + s):(sum(n_bal) + s + n_bal[i] - 1)] - samples[j, s:(s + n_bal[i] - 1)])
      b2 <- (1 / (sig[i])^2) * t(rep(1, n_bal[i])) %*% (Cj - Ej) %*% matrix(1, nrow = n_bal[i], ncol = 1)
      
      s <- s + n_bal[i]
    }
    beta_hat <- sum(b1) / sum(b2)
    

    t_stat <- (beta_hat) / sqrt(var_n_bal)
    p_value <- 2 * (1 - pt(abs(t_stat), df = (length(n_bal)) - 1))  # Two-tailed p-value
    if (p_value < alpha) {
      ctr_bal <- ctr_bal + 1
    }
  }
  Pwr_bal[v] <- ctr_bal / n_samples
}

# Boxplot for comparison
boxplot(Pwr_opt, Pwr_bal,
        names = c("LOD", "Balanced"),
        col = c("lightblue", "pink"),
        #main = "Power comparison of the designs",
        ylab = "Power",
        xlab = "Design",
        border = "darkblue")
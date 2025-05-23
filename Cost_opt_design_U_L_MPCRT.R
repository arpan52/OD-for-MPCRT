rm(list = ls())  # Remove all existing objects from the R environment

optimal_design_MPCRT_varying_cost<- function(T,L,U,rho,rho1,sigma,c){


# Initialization
T_star <- T / 2  # Budget per arm
a <- rep(0, length(c))  # Will store 1 - rho - rho1 for each cluster
x1 <- rep(0, length(c)) # Temporary storage for formula components
y1 <- rep(0, length(c))

# Calculate a[i] = 1 - rho[i] - rho1[i]
for (i in 1:length(c)) {
  a[i] <- 1 - rho[i] - rho1[i]
}

# Compute x1 and y1 used to find unconstrained optimal n_new
for (i in 1:length(c)) {
  x1[i] <- (sqrt(c[i] * a[i])) / (sigma[i] * rho[i])
  y1[i] <- (c[i] * a[i]) / rho[i]
}

x <- sum(x1)  # Aggregate x1
y <- sum(y1)  # Aggregate y1

numbers <- 1:length(c)

# Generate all combinations of indices to assign L or U to some cluster sizes
all_combinations <- list()
c_n <- rep(0, length(c))

for (k in 1:length(numbers)) {
  all_combinations[[k]] <- combn(numbers, k, simplify = FALSE)
}

# Compute length of design space
len <- 0
for(i in 1:(length(all_combinations) - 1)){
  len <- len + length(all_combinations[[i]])
}

n_new <- rep(0, length(c))     # Store current design
n_til <- rep(0, length(c))     # Will store all valid designs

# Compute unconstrained optimal design
for (i in 1:(length(c))) {
  n_new[i] <- (T_star - ((x * sigma[i] * sqrt(c[i] * a[i])) - y)) / (x * sigma[i] * rho[i] * sqrt(c[i] / a[i]))
}

# Check if unconstrained design satisfies constraints
for (l in 1:length(c)) {
  c_n[l] <- c[l] * n_new[l]
}

if (all(n_new >= L & n_new <= U)) {
  if (sum(c_n) >= T_star - 2 & sum(c_n) <= T_star + 2) {
    n_til <- rbind(n_til, n_new)
  }
}

# Consider combinations where some clusters have n_j = L
# and optimize for remaining clusters
for (i in 1:(length(all_combinations) - 1)) {
  for (j in 1:length(all_combinations[[i]])) {
    for (k1 in 1:length(all_combinations[[i]][[j]])) {
      n_new[all_combinations[[i]][[j]][k1]] <- L
    }
    set <- setdiff(1:length(c), all_combinations[[i]][[j]])
    
    c_minus <- sum(c[all_combinations[[i]][[j]]])
    T_star1 <- T_star - (L * c_minus)
    
    xl1 <- sapply(set, function(l1) (sqrt(c[l1] * a[l1])) / (sigma[l1] * rho[l1]))
    yl1 <- sapply(set, function(l1) (c[l1] * a[l1]) / rho[l1])
    x <- sum(xl1)
    y <- sum(yl1)
    
    for (k2 in 1:length(set)) {
      n_new[set[k2]] <- (T_star1 - ((x * sigma[set[k2]] * sqrt(c[set[k2]] * a[set[k2]])) - y)) / 
        (x * sigma[set[k2]] * rho[set[k2]] * sqrt(c[set[k2]] / a[set[k2]]))
    }
    
    for (l in 1:length(c)) {
      c_n[l] <- c[l] * n_new[l]
    }
    if (all(n_new >= L & n_new <= U)) {
      if (sum(c_n) >= T_star - 2 & sum(c_n) <= T_star + 2) {
        n_til <- rbind(n_til, n_new)
      }
    }
  }
}

# Repeat the same process where some clusters are fixed at U
for (i in 1:(length(all_combinations) - 1)) {
  for (j in 1:length(all_combinations[[i]])) {
    for (k1 in 1:length(all_combinations[[i]][[j]])) {
      n_new[all_combinations[[i]][[j]][k1]] <- U
    }
    set <- setdiff(1:length(c), all_combinations[[i]][[j]])
    
    c_minus <- sum(c[all_combinations[[i]][[j]]])
    T_star1 <- T_star - (U * c_minus)
    
    xl1 <- sapply(set, function(l1) (sqrt(c[l1] * a[l1])) / (sigma[l1] * rho[l1]))
    yl1 <- sapply(set, function(l1) (c[l1] * a[l1]) / rho[l1])
    x <- sum(xl1)
    y <- sum(yl1)
    
    for (k2 in 1:length(set)) {
      n_new[set[k2]] <- (T_star1 - ((x * sigma[set[k2]] * sqrt(c[set[k2]] * a[set[k2]])) - y)) / 
        (x * sigma[set[k2]] * rho[set[k2]] * sqrt(c[set[k2]] / a[set[k2]]))
    }
    
    for (l in 1:length(c)) {
      c_n[l] <- c[l] * n_new[l]
    }
    if (all(n_new >= L & n_new <= U)) {
      if (sum(c_n) >= T_star - 2 & sum(c_n) <= T_star + 2) {
        n_til <- rbind(n_til, n_new)
      }
    }
  }
}

# Combinations with clusters fixed at L and remaining at U
for (i in 1:(length(all_combinations) - 1)) {
  for (j in 1:length(all_combinations[[i]])) {
    for (k1 in 1:length(all_combinations[[i]][[j]])) {
      n_new[all_combinations[[i]][[j]][k1]] <- L
    }
    set <- setdiff(1:length(c), all_combinations[[i]][[j]])
    for (k2 in 1:length(set)) {
      n_new[set[k2]] <- U
    }
    for (l in 1:length(c)) {
      c_n[l] <- c[l] * n_new[l]
    }
    if (all(n_new >= L & n_new <= U)) {
      if (sum(c_n) >= T_star - 2 & sum(c_n) <= T_star + 2) {
        n_til <- rbind(n_til, n_new)
      }
    }
  }
}

# Generate mixed designs where some clusters are L, some U, and others optimized
for (i in 1:(length(all_combinations) - 2)) {
  for (j in 1:length(all_combinations[[i]])) {
    for (k1 in 1:length(all_combinations[[i]][[j]])) {
      n_new[all_combinations[[i]][[j]][k1]] <- L
      rem_set <- setdiff(1:length(c), all_combinations[[i]][[j]])
      
      all_combinations_rem <- list()
      for (u in 1:length(rem_set)) {
        all_combinations_rem[[u]] <- combn(rem_set, u, simplify = FALSE)
      }
      
      for (i1 in 1:(length(all_combinations_rem) - 1)) {
        for (j1 in 1:length(all_combinations_rem[[i1]])) {
          for (k3 in 1:length(all_combinations_rem[[i1]][[j1]])) {
            n_new[all_combinations_rem[[i1]][[j1]][k3]] <- U
          }
          
          rem_set_new <- setdiff(rem_set, all_combinations_rem[[i1]][[j1]])
          
          c_minus_new1 <- sum(c[all_combinations[[i]][[j]]])
          c_minus_new2 <- sum(c[all_combinations_rem[[i1]][[j1]]])
          T_star1 <- T_star - (L * c_minus_new1) - (U * c_minus_new2)
          
          xl1 <- sapply(rem_set_new, function(l1) (sqrt(c[l1] * a[l1])) / (sigma[l1] * rho[l1]))
          yl1 <- sapply(rem_set_new, function(l1) (c[l1] * a[l1]) / rho[l1])
          x <- sum(xl1)
          y <- sum(yl1)
          
          for (k4 in 1:length(rem_set_new)) {
            n_new[rem_set_new[k4]] <- (T_star1 - ((x * sigma[rem_set_new[k4]] * sqrt(c[rem_set_new[k4]] * a[rem_set_new[k4]])) - y)) /
              (x * sigma[rem_set_new[k4]] * rho[rem_set_new[k4]] * sqrt(c[rem_set_new[k4]] / a[rem_set_new[k4]]))
          }
          
          for (l in 1:length(c)) {
            c_n[l] <- c[l] * n_new[l]
          }
          if (all(n_new >= L & n_new <= U)) {
            if (sum(c_n) >= T_star - 2 & sum(c_n) <= T_star + 2) {
              n_til <- rbind(n_til, n_new)
            }
          }
        }
      }
    }
  }
}

# Stop if no valid designs were found
if (all(n_til == 0)){
  stop("Lower and Upper bounds are extremely narrow")
}

# Remove initial placeholder row
n_til <- n_til[-1, ]


v<-rep(0,length(c))
var_n_til<-rep(0,nrow(n_til))
for (i in 1: nrow(n_til)){
  for (j in 1:(length(c))){
    v[j] <- (n_til[i,j])/((sigma[j]^2)*(1+((n_til[i,j]-1)*(rho[j])-rho1[j])))
  }
  var_n_til[i]<- 2/sum(v)
}
min_var<- which.min(var_n_til)

opt_design<- n_til[min_var,]

n_bal<-rep(T_star/sum(c),length(c))
v_bal<-rep(0,length(c))
for (j in 1:(length(c))){
  v_bal[j] <- (n_bal[j])/((sigma[j]^2)*(1+((n_bal[j]-1)*(rho[j])-rho1[j])))
}
var_n_bal<- 2/sum(v_bal)


print('Optimal Design containg the vector of cluster sizes of the paired clusters')
print(round(opt_design))

print('Total Subjects recruited')
print(2*sum(round(opt_design)))

print('Variance of the optimal design')
print(min(var_n_til))

print('Efficiency of the Optimal design comapred to the Balanced design')
print(var_n_bal/min(var_n_til))

}


###########################################################################################################

# Set up model parameters
T <- 240      # Total budget constraint for two arms
L <- 7        # Lower bound for cluster sizes
U <- 25       # Upper bound for cluster sizes

# Vectors defining characteristics of each of the m = 12 clusters
rho <- c(rep(0.014, 12))          # Intra-class correlations (ICCs)
rho1 <- c(rep(0.4, 12))           # Matching correlations
sigma <- rep(1.8, 12)             # Standard deviations
c <- c(rep(0.5,3), rep(1,6), rep(1.5,3))  # Sampling costs

optimal_design_MPCRT_varying_cost(T,L,U,rho,rho1,sigma,c)   # compute the output of Optimal Design
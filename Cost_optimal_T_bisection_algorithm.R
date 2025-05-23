rm(list = ls())


T_u<-200
T_l<-100
T_star<- (T_u+T_l)/2
p0<-0.8
p<-0.1
Iteration<-1
print('Iteration')
print(Iteration)


rho<- c(rep(0.014,3),rep(0.014,3),rep(0.14,3),rep(0.014,3))
rho1<- c(rep(0.7,3),rep(0.75,3),rep(0.85,3),rep(0.9,3))
sigma<- rep(1.8,12)
c<-c(rep(0.8,3),rep(1.1,3),rep(1.1,3),rep(1,3))
m<-length(c)

# lower and upper bounds

L<-7
U<-25
al<-0.05

t1 = qt(al/2,m-1)
t2 = qt(1-(al/2),m-1)
delta<-0.4

while(abs(p-p0)>0.005){
  
  T_star<- (T_u+T_l)/2
  
  print('Interval T')
  print(c(T_l,T_u))
  
  print('Budget')
  print(T_star)
  
  
  a <- rep(0,length(c))
  x1<- rep(0,length(c))
  y1<- rep(0,length(c))
  
  for (i in 1:length(c)) {
    a[i]<- 1-rho[i]-rho1[i]
  }
  
  for (i in 1:length(c)){
    x1[i]<- (sqrt(c[i]*a[i]))/(sigma[i]*rho[i])
  }
  
  x<-sum(x1)
  
  for (i in 1:length(c)) {
    y1[i]<-(c[i]*a[i])/rho[i]
  }
  y<-sum(y1)
  
  
  numbers <- 1:length(c)
  # Generate all the combinations where the cluster size can be L or U combinations
  all_combinations <- list()
  c_n<-rep(0,length(c))
  
  for (k in 1:length(numbers)) {
    all_combinations[[k]] <- combn(numbers, k, simplify = FALSE)
  }
  
  # Find the total number of design vectors n
  len<-0
  for(i in 1:(length(all_combinations)-1)){
    len<-len+length(all_combinations[[i]])
  }
  
  n_new<-rep(0,length(c))
  
  n_til<-rep(0,length(c))        # initialize the optimal design without bounds
  
  for (i in 1:(length(c))){
    n_new[i]<- (T_star-((x*sigma[i]*sqrt(c[i]*a[i]))-y))/(x*sigma[i]*rho[i]*sqrt(c[i]/a[i]))
    
  }
  #n_new<-round(n_new)
  # check for KKT conditions for the vector n
  
  for (l in 1:(length(c))){
    c_n[l]<- c[l]* n_new[l]
  }
  if (all(n_new >= L & n_new <= U)) {
    
    if(sum(c_n)>=T_star-2 & sum(c_n)<=T_star+2){
      n_til <- rbind(n_til, n_new)
    }
  }
  
  for (i in 1:(length(all_combinations)-1)){
    
    for (j in 1:length(all_combinations[[i]])){
      
      for (k1 in 1:length(all_combinations[[i]][[j]])){
        
        n_new[all_combinations[[i]][[j]][k1]]<- L                   # Those n_j's which are equal to L
        
      }
      set<- setdiff(1:length(c), all_combinations[[i]][[j]])        # remaining n_j's
      
      c_minus<- 0
      for (p in 1:length(all_combinations[[i]][[j]])){
        c_minus<- c_minus+c[all_combinations[[i]][[j]][p]]
      }
      
      
      T_star1<- T_star-(L*c_minus)                                 # new T_star for remaining n_j's
      xl1<-rep(0,length(set))
      yl1<-rep(0,length(set))
      for (l1 in 1:length(set)){
        xl1[l1]<- (sqrt(c[set[l1]]*a[set[l1]]))/(sigma[set[l1]]*rho[set[l1]])
      }
      
      x<-sum(xl1)
      for (l1 in 1:length(set)) {
        yl1[l1]<-(c[set[l1]]*a[set[l1]])/rho[set[l1]]
      }
      y<-sum(yl1)
      
      for (k2 in 1:length(set)){                                   # formula for remaining n_j's
        
        n_new[set[k2]]<- (T_star1-((x*sigma[set[k2]]*sqrt(c[set[k2]]*a[set[k2]]))-y))/(x*sigma[set[k2]]*rho[set[k2]]*sqrt(c[set[k2]]/a[set[k2]]))
        
      }
      #n_new<-round(n_new)
      # check for KKT conditions for the vector n
      
      for (l in 1:(length(c))){
        c_n[l]<- c[l]* n_new[l]
      }
      if (all(n_new >= L & n_new <= U)) {
        
        if(sum(c_n)>=T_star-2 & sum(c_n)<=T_star+2){
          n_til <- rbind(n_til, n_new)
        }
      }
      
      
    }
  }
  
  for (i in 1:(length(all_combinations)-1)){
    
    for (j in 1:length(all_combinations[[i]])){
      
      for (k1 in 1:length(all_combinations[[i]][[j]])){
        
        n_new[all_combinations[[i]][[j]][k1]]<- U                   # Those n_j's which are equal to U
        
      }
      set<- setdiff(1:length(c), all_combinations[[i]][[j]])        # remaining n_j's
      
      c_minus<- 0
      for (p in 1:length(all_combinations[[i]][[j]])){
        c_minus<- c_minus+c[all_combinations[[i]][[j]][p]]
      }
      T_star1<- T_star-(U*c_minus)                                 # new T_star for remaining n_j's
      
      xl1<-rep(0,length(set))
      yl1<-rep(0,length(set))
      for (l1 in 1:length(set)){
        xl1[l1]<- (sqrt(c[set[l1]]*a[set[l1]]))/(sigma[set[l1]]*rho[set[l1]])
      }
      
      x<-sum(xl1)
      for (l1 in 1:length(set)) {
        yl1[l1]<-(c[set[l1]]*a[set[l1]])/rho[set[l1]]
      }
      y<-sum(yl1)
      
      
      
      for (k2 in 1:length(set)){                                   # formula for remaining n_j's
        
        n_new[set[k2]]<- (T_star1-((x*sigma[set[k2]]*sqrt(c[set[k2]]*a[set[k2]]))-y))/(x*sigma[set[k2]]*rho[set[k2]]*sqrt(c[set[k2]]/a[set[k2]]))
        
      }
      #n_new<-round(n_new)
      
      # check for KKT conditions for the vector n
      
      for (l in 1:(length(c))){
        c_n[l]<- c[l]* n_new[l]
      }
      if (all(n_new >= L & n_new <= U)) {
        
        if(sum(c_n)>=T_star-2 & sum(c_n)<=T_star+2){
          n_til <- rbind(n_til, n_new)
        }
      }
      
      
    }
  }
  
  for (i in 1:(length(all_combinations)-1)){
    
    for (j in 1:length(all_combinations[[i]])){
      
      for (k1 in 1:length(all_combinations[[i]][[j]])){
        
        n_new[all_combinations[[i]][[j]][k1]]<- L                   # Those n_j's which are equal to L
        
      }
      set<- setdiff(1:length(c), all_combinations[[i]][[j]])        # remaining n_j's
      
      for (k2 in 1:length(set)){                                   # formula for remaining n_j's
        
        n_new[set[k2]]<- U
        
      }
      
      #n_new<-round(n_new)
      
      # check for KKT conditions for the vector n
      for (l in 1:(length(c))){
        c_n[l]<- c[l]* n_new[l]
      }
      if (all(n_new >= L & n_new <= U)) {
        
        if(sum(c_n)>=T_star-2 & sum(c_n)<=T_star+2){
          n_til <- rbind(n_til, n_new)
        }
      }
      
      
    }
  }
  
  # Those n's which contains all three type of n_j's
  
  for (i in 1:(length(all_combinations)-2)){
    for (j in 1:length(all_combinations[[i]])){
      for (k1 in 1:length(all_combinations[[i]][[j]])){
        
        # n[i,j,all_combinations[[i]][[j]][k1]]<- L                   # Those n_j's which are equal to L
        n_new[all_combinations[[i]][[j]][k1]]<- L                   # Those n_j's which are equal to L
        
        rem_set<- setdiff(1:length(c), all_combinations[[i]][[j]])        # remaining n_j's
        
        
        
        # Generate all the combinations where the cluster size can be L or U combinations
        all_combinations_rem <- list()
        
        for (u in 1:length(rem_set)) {
          all_combinations_rem[[u]] <- combn(rem_set, u, simplify = FALSE)
        }  
        
        # assign U and n_tilde to the remaining n_j's
        for (i1 in 1:(length(all_combinations_rem)-1)){
          
          for (j1 in 1:length(all_combinations_rem[[i1]])){
            
            for (k3 in 1:length(all_combinations_rem[[i1]][[j1]])){
              
              n_new[all_combinations_rem[[i1]][[j1]][k3]]<- U                   # Those n_j's which are equal to L
              
            }
            
            
            rem_set_new<- setdiff(rem_set, all_combinations_rem[[i1]][[j1]])        # remaining n_j's
            
            c_minus_new1<- 0
            c_minus_new2<- 0
            
            for (p in 1:length(all_combinations[[i]][[j]])){
              c_minus_new1<- c_minus_new1+c[all_combinations[[i]][[j]][p]]
            }
            for (p in 1:length(all_combinations_rem[[i1]][[j1]])){
              c_minus_new2<- c_minus_new2+c[all_combinations_rem[[i1]][[j1]][p]]
            }
            
            T_star1<- T_star-(L*c_minus_new1)-(U*c_minus_new2) # new T_star for remaining n_j's
            
            xl1<-rep(0,length(rem_set_new))
            yl1<-rep(0,length(rem_set_new))
            for (l1 in 1:length(rem_set_new)){
              xl1[l1]<- (sqrt(c[rem_set_new[l1]]*a[rem_set_new[l1]]))/(sigma[rem_set_new[l1]]*rho[rem_set_new[l1]])
            }
            
            x<-sum(xl1)
            for (l1 in 1:length(rem_set_new)) {
              yl1[l1]<-(c[rem_set_new[l1]]*a[rem_set_new[l1]])/rho[rem_set_new[l1]]
            }
            y<-sum(yl1)
            
            for (k4 in 1:length(rem_set_new)){                                   # formula for remaining n_j's
              
              n_new[rem_set_new[k4]]<- (T_star1-((x*sigma[rem_set_new[k4]]*sqrt(c[rem_set_new[k4]]*a[rem_set_new[k4]]))-y))/(x*sigma[rem_set_new[k4]]*rho[rem_set_new[k4]]*sqrt(c[rem_set_new[k4]]/a[rem_set_new[k4]]))
              
            }
            
            #n_new<-round(n_new)
            
            # check for KKT conditions for the vector n
            for (l in 1:(length(c))){
              c_n[l]<- c[l]* n_new[l]
            }
            
            if (all(n_new >= L & n_new <= U)) {
              
              if(sum(c_n)>=T_star-2 & sum(c_n)<=T_star+2){
                n_til <- rbind(n_til, n_new)
              }
            }
            
            
          }
        }  
        
      }
      
      
      
    }
  }
  
  if (all(n_til == 0)){
    stop("Lower and Upper bounds are extremely narrow")
  }
  
  n_til <- n_til[-1, ]
  
  v<-rep(0,length(c))
  var_n_til<-rep(0,nrow(n_til))
  for (i in 1: nrow(n_til)){
    for (j in 1:(length(c))){
      v[j] <- (n_til[i,j])/((sigma[j]^2)*(1+((n_til[i,j]-1)*(rho[j])-rho1[j])))
    }
    var_n_til[i]<- 2/sum(v)
  }
  
  min_des<- which.min(var_n_til)
  
  opt_design<- n_til[min_des,]
  print('Optimal Design')
  print(opt_design)
  print('Total Subjects')
  print(sum(round(opt_design)))
  
  min_var<- min(var_n_til)
  
  pwr_opt <- 1 - pt(t2, df = m - 1, ncp = delta / sqrt(min_var)) +
    pt(t1, df = m - 1, ncp = delta / sqrt(min_var))
  
  if(pwr_opt>=0.8){T_u<- (T_l+T_u)/2}
  if(pwr_opt<0.8){T_l<- (T_l+T_u)/2}
  
  p<-pwr_opt
  
  print('Iteration')
  Iteration<-Iteration+1
  print(Iteration)
  
  print('Power')
  print(p)
  
  
}
print('Optimal Design')
min_v<- which.min(var_n_til)
opt_design<- n_til[min_v,]
print(round(opt_design))


print('Optimal Budget')
print(T_star)
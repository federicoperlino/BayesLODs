######### Bayesian MI Models for data suffering from Limit of Detections ######

# Packages ----------------------------------------------------------------

library(plyr)
library(tmvtnorm)

# Preparatory Functions --------------------------------------------------
# As in Hoff, Peter D. A first course in Bayesian statistical methods (2009) 

rmvnorm <- function(n,mu,Sigma) {
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  res
}

rwish <- function(n,nu0,S0){
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}



# MVN Model ---------------------------------------------------------------

gibbs_mvn <- function(Z, mu_0, Omega, v0, S_0, cens, R, ag = "gibbs"){
  # Fix starting values for truncated sampler based on censoring thresholds
  fix_start_value <- function(cens_val) {
    if (cens_val >= 0) {
      return(cens_val / sqrt(2))
    } else {
      return(cens_val - abs(cens_val) / sqrt(2))
    }
  }
  
  # Impute missing entries of z using a truncated multivariate normal
  tnorm.impute <- function(z, mu, Sigma, cens){
    id.na <- which(is.na(z))
    if(length(id.na) > 0){
      # Conditional mean and covariance for the missing components
      mu.cond <- mu[id.na] + Sigma[id.na,-id.na] %*% solve(Sigma[-id.na, -id.na]) %*% (z[-id.na] - mu[-id.na])
      Sigma.cond <- Sigma[id.na, id.na] - Sigma[id.na, -id.na] %*% 
        solve(Sigma[-id.na, -id.na]) %*% Sigma[-id.na, id.na]
      
      # Compute start values for sampler using fix_start_value
      start.val <- sapply(cens[id.na], fix_start_value)
      
      # Draw from truncated multivariate normal up to the censoring limits
      z[id.na] <- rtmvnorm(
        n = 1,
        mean = as.numeric(mu.cond),
        sigma = Sigma.cond,
        lower = rep(-Inf, length(mu.cond)),
        upper = cens[id.na],
        algorithm = ag,
        start.value = start.val,
        thinning = 2
      )
    }
    return(z)
  }
  
  # Dimensions: q variables, n observations
  q <- ncol(Z)
  n <- nrow(Z)
  
  # Array Y to store imputed values across R iterations
  Y <- array(dim = c(n, q, R))
  
  #--------------------------------------------------
  # Initialization: draw initial mu and Sigma
  mu <- rnorm(q)                             # Initial mean vector
  Sigma <- diag(runif(1, 0, 2), q)          # Initial covariance matrix
  
  if (is.matrix(cens)){
    # Case: individual-specific censoring limits
    # First imputation based on initial parameter values
    Y[,,1] <- sapply(1:nrow(Z), function(x)
      tnorm.impute(Z[x, ], mu, Sigma, cens[x, ])
    ) |> t()
    
    # Gibbs sampling iterations
    for(s in 2:R){
      # Update mu
      ybar <- apply(Y[,,s-1], 2, mean)
      A_n <- solve(Omega) + n * solve(Sigma)
      b_n <- solve(Omega) %*% mu_0 + n * solve(Sigma) %*% ybar
      mu <- rmvnorm(1, solve(A_n) %*% b_n, solve(A_n))
      # Note: rmvnorm interprets the third argument as covariance
      
      # Update Sigma
      S_mu <- (t(Y[,,s-1]) - c(mu)) %*% t(t(Y[,,s-1]) - c(mu))
      S_n <- solve(S_0 + S_mu)
      Sigma <- solve(rwish(1, v0 + n, S_n))
      
      # Imputation step
      Y[,,s] <- sapply(1:nrow(Z), function(x)
        tnorm.impute(Z[x, ], mu, Sigma, cens[x, ])
      ) |> t()
    }
  } else {
    # Case: censoring limits common across observations
    Y[,,1] <- apply(Z, 1, tnorm.impute, mu = mu, Sigma = Sigma, cens = cens) |> t()
    
    for(s in 2:R){
      
      # Update mu
      ybar <- apply(Y[,,s-1], 2, mean)
      A_n <- solve(Omega) + n * solve(Sigma)
      b_n <- solve(Omega) %*% mu_0 + n * solve(Sigma) %*% ybar
      mu <- rmvnorm(1, solve(A_n) %*% b_n, solve(A_n))
      
      # Update Sigma
      S_mu <- (t(Y[,,s-1]) - c(mu)) %*% t(t(Y[,,s-1]) - c(mu))
      S_n <- solve(S_0 + S_mu)
      Sigma <- solve(rwish(1, v0 + n, S_n))
      
      # Imputation step
      Y[,,s] <- apply(Z, 1, tnorm.impute, mu = mu, Sigma = Sigma, cens = cens) |> t()
    }
  }
  
  output <- list(Y = Y, mu = mu, Sigma = Sigma)
  return(output)
}


# DPML Model --------------------------------------------------------------

gibbs_DPML <- function(Z, mu_0, Sigma_0, v0, S_0, c, cens, R, acc = T, 
                       ag = "gibbs") {
 
  # Auxiliary function to sample from theta via Pólya urn
  # Returns: 
  # a vector of theta for each observation
  
  update_individual <- function(i, theta, Y, Sigma, Sigma_0, mu_0, q, n, s, 
                                logc, logcn_1, solve_Sigma_0, log_det_Sigma_0, 
                                solve_Sigma, log_det_Sigma) {
    thetai <- theta[i, ]
    theta_meno_i <- theta[-i, ]
    Yi <- Y[i, , s - 1]
    COUNT <- plyr::count(theta_meno_i)
    theta_meno_i_star <- as.matrix(COUNT[, 1:q])
    n_meno_i_star <- COUNT$freq
    q_meno_i_star <- nrow(COUNT)
    log_pi_vec_temp <- rep(0, q_meno_i_star + 1)
    
    Sigma_i <- solve(solve_Sigma_0 + solve_Sigma)
    mu_i <- Sigma_i %*% (solve_Sigma %*% Yi + solve_Sigma_0 %*% mu_0)
    
    log_pi_vec_temp[1] <- logc - logcn_1 - (q / 2) * log(2 * pi) +
      (0.5) * (log(det(Sigma_i)) - log_det_Sigma - log_det_Sigma_0) -
      (0.5) * (t(Yi) %*% solve_Sigma %*% Yi) -
      (0.5) * (t(mu_0) %*% solve_Sigma_0 %*% mu_0) +
      (0.5) * (t(mu_i) %*% solve(Sigma_i) %*% mu_i)
    
    log_pi_vec_temp[2:(q_meno_i_star + 1)] <- log(n_meno_i_star) - logcn_1 +
      sapply(1:nrow(theta_meno_i_star), function(j)
        dmvnorm(Yi, theta_meno_i_star[j, ], Sigma, log = TRUE))
    
    pi_vec <- exp(log_pi_vec_temp - max(log_pi_vec_temp))
    pi_vec <- pi_vec / sum(pi_vec)
    
    index <- sample.int(q_meno_i_star + 1, size = 1, prob = pi_vec)
    if (index == 1) {
      theta[i, ] <- rmvnorm(1, mu_i, Sigma_i)
    } else {
      theta[i, ] <- theta_meno_i_star[index - 1, ]
    }
    return(theta[i, ])
  }
  
  fix_start_value <- function(cens_val) {
    if (cens_val >= 0) {
      return(cens_val / sqrt(2))
    } else {
      return(cens_val - abs(cens_val) / sqrt(2))
    }
  }
  
  # Auxiliary function to impute values below LODs in the DPML model 
  # Returns:
  # Vector (row of the Y matrix) with imputed values 
  
  tnorm.mixture.impute <- function(z, theta, Sigma, cens) {
    id.na <- which(is.na(z))
    start.val <- sapply(cens[id.na], fix_start_value)
    
    if (length(id.na) > 0) {
      Sigma_inv <- solve(Sigma[-id.na, -id.na])
      theta.cond <- theta[id.na] + Sigma[id.na, -id.na] %*% 
        Sigma_inv %*% (z[-id.na] - theta[-id.na])
      Sigma.cond <- Sigma[id.na, id.na] - Sigma[id.na, -id.na] %*% 
        Sigma_inv %*% Sigma[-id.na, id.na]
      z[id.na] <- rtmvnorm(n = 1, 
                           lower = rep(-Inf, length(theta.cond)),
                           upper = cens[id.na], 
                           mean = as.numeric(theta.cond), 
                           sigma = Sigma.cond,
                           thinning = 2, algorithm = ag, 
                           start.value = start.val)
    }
    return(z)
  }
  
  n <- nrow(Z)
  q <- ncol(Z)
  
  Y <- array(dim = c(n, q, R))
  theta_store <- array(dim = c(n, q, R))
  Sigma_store <- array(dim = c(q, q, R))
  
  # Initialization
    theta <- matrix(rnorm(n * q), n, q)
    Sigma <- diag(runif(1, 0, 2), q)
  
  if (is.matrix(cens)) {  # i.e., individual-specific LODs
    Y[,,1] <- sapply(1:n, function(x)
      tnorm.mixture.impute(Z[x,], theta = theta[x,], 
                           Sigma = Sigma, cens = cens[x,])
    ) |> t()
    theta_store[,,1] <- theta
    Sigma_store[,,1] <- Sigma
    
    logc <- log(c)
    logcn_1 <- log(c + n - 1)
    solve_Sigma_0 <- solve(Sigma_0)
    log_det_Sigma_0 <- log(det(Sigma_0))
    
    # Gibbs steps 
    for (s in 2:R) {
      solve_Sigma <- solve(Sigma)
      log_det_Sigma <- log(det(Sigma))
      
      # Update theta
      for (i in 1:n) {
        theta[i, ] <- update_individual(i, theta, Y, Sigma, Sigma_0, 
                                        mu_0, q, n, s,
                                        logc, logcn_1, solve_Sigma_0, 
                                        log_det_Sigma_0, solve_Sigma, 
                                        log_det_Sigma)
      }
      
      # Update Sigma
      S_theta <- t(Y[,,s - 1] - theta) %*% (Y[,,s - 1] - theta)
      S_n <- solve(S_0 + S_theta)
      Sigma <- solve(rwish(1, v0 + n, S_n))
      
      # Imputation step
      Y[,,s] <- sapply(1:n, function(x)
        tnorm.mixture.impute(Z[x,], theta = theta[x,], 
                             Sigma = Sigma, cens = cens[x,])
      ) |> t()
      
      # Acceleration step
      if (acc == T) {
        COUNT <- plyr::count(theta[,])
        theta_star <- as.matrix(COUNT[, 1:q])
        n_star <- COUNT$freq
        q_star <- nrow(theta_star)
        
        for (j in 1:q_star) {
          thetaj_star <- theta_star[j, ]
          indexj <- which(apply(theta[,], 1, function(x) 
            identical(x, thetaj_star)))
          indexj <- ifelse(length(indexj) == 0, j, indexj)
          ybarj <- if (length(indexj) > 1) apply(Y[indexj, , s], 2, mean) 
          else Y[indexj, , s]
          Cj <- length(indexj)
          
          A_n <- solve(Sigma_0) + Cj * solve(Sigma)
          b_n <- solve(Sigma_0) %*% mu_0 + Cj * solve(Sigma) %*% ybarj
          
          thetaj_star <- rmvnorm(1, solve(A_n) %*% b_n, solve(A_n))
          theta[indexj, ] <- thetaj_star
        }
      }
      
      # Storage
      theta_store[,,s] <- theta
      Sigma_store[,,s] <- Sigma
    }
    
  } else {  # i.e., LODs varying only across the exposures
    Y[,,1] <- sapply(1:n, function(x)
      tnorm.mixture.impute(Z[x,], theta = theta[x,], Sigma = Sigma, 
                           cens = cens)
    ) |> t()
    
    logc <- log(c)
    logcn_1 <- log(c + n - 1)
    solve_Sigma_0 <- solve(Sigma_0)
    log_det_Sigma_0 <- log(det(Sigma_0))
    
    
    # Gibbs steps
    for (s in 2:R) {
      solve_Sigma <- solve(Sigma)
      log_det_Sigma <- log(det(Sigma))
      
      # Update theta
      for (i in 1:n) {
        theta[i, ] <- update_individual(i, theta, Y, Sigma, Sigma_0, mu_0,
                                        q, n, s, logc, logcn_1, 
                                        solve_Sigma_0, log_det_Sigma_0, 
                                        solve_Sigma, log_det_Sigma)
      }
      
      # Update Sigma
      S_theta <- t(Y[,,s - 1] - theta) %*% (Y[,,s - 1] - theta)
      S_n <- solve(S_0 + S_theta)
      Sigma <- solve(rwish(1, v0 + n, S_n))
      
      # Imputation step
      Y[,,s] <- sapply(1:n, function(x)
        tnorm.mixture.impute(Z[x,], theta = theta[x,], 
                             Sigma = Sigma, cens = cens)
      ) |> t()
      
      # Acceleration step
      if (acc == T) {
        COUNT <- plyr::count(theta[,])
        theta_star <- as.matrix(COUNT[, 1:q])
        n_star <- COUNT$freq
        q_star <- nrow(theta_star)
        
        for (j in 1:q_star) {
          thetaj_star <- theta_star[j, ]
          indexj <- which(apply(theta[,], 1, function(x) identical(x, thetaj_star)))
          indexj <- ifelse(length(indexj) == 0, j, indexj)
          ybarj <- if (length(indexj) > 1) apply(Y[indexj, , s], 2, mean) else Y[indexj, , s]
          Cj <- length(indexj)
          
          A_n <- solve(Sigma_0) + Cj * solve(Sigma)
          b_n <- solve(Sigma_0) %*% mu_0 + Cj * solve(Sigma) %*% ybarj
          
          thetaj_star <- rmvnorm(1, solve(A_n) %*% b_n, solve(A_n))
          theta[indexj, ] <- thetaj_star
        }
      }
      
      # Storage
      theta_store[,,s] <- theta
      Sigma_store[,,s] <- Sigma
      
      if (s %% 100 == 0) {
        cat("end iter", s, "\n")
      }
    }
  }
  
  return(list(Y = Y, theta_store = theta_store, Sigma_store = Sigma_store))
}





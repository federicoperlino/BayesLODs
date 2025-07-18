
# Function for performing prior elicitation as in reference paper ---------

standard_prior_elicitation <- function(model, Z, Ek = 1, scale = T){
  
  # Only for the DPML model 
  c_value <- function(Ek, N){
    
    # Ek is the expected number of cluster. 
    # Give a guess on this based on the characteristic of the problem
    
    objective_function <- function(c, Ek, N) {
      sum_result <- sum(c / (1:(N)-1 + c))
      return(abs(sum_result - Ek))
    }
    
    result <- optim(par = 1, fn = objective_function, method = "L-BFGS-B", lower = 1e-6, upper = 10^9, N = N, Ek= Ek)
    return(result$par)
  }
  
  # Function for "standard" (as in the main paper) prior elicitation of the MI MVN and DPML models
  
  # IW
  q <- ncol(Z)
  v0 <- q+4
  S0 <- (v0 - q - 1)*cov(Z, use  = "complete.obs")
  
  # MVN
  if(scale == T){
    mu <- rep(0,q)
    Omega <- diag(2.5, q)
  }
  if(scale == F){
    mu <- apply(Z, 2, mean, na.rm =T)
    Omega <- diag(apply(Z, 2, var, na.rm = T), q)
  }
  
  
  if(model == "DPML"){
    cvalue<-c_value(Ek, nrow(Z))
    out <- list(mu_G0 = mu, Sigma_G0 = Omega, c = cvalue, v0 = v0, S0 = S0   )
  }
  
  if(model == "MVN"){
    out <- list(mu = mu, Omega = Omega, v0 = v0, S0 = S0 )
  }
  
  return(out)
}

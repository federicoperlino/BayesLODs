##### Toy Example #####

set.seed(1234)
Y <- rmvnorm(n = 200, mu = rep(1, 8), Sigma = diag(1,8))
LODs <- apply(Y, 2, quantile, probs = 0.25)
IS_LODs <- t(sapply(1:nrow(Y), function(x) apply(Y, 2, quantile, probs = 0.25))) +
  rnorm(nrow(Y))

Z_LODs<-sapply(1:ncol(Y), function(j) ifelse(Y[ ,j  ] <= LODs[j], NA, Y[,j]  )) 


Z_IS_LODs_full <- ifelse(Y <= IS_LODs, NA, Y)

# Indexes of rows with all NA
keep <- !apply(Z_IS_LODs_full, 1, function(r) all(is.na(r)))

# Remove rows having all NA

Z_IS_LODs   <- Z_IS_LODs_full[keep, ]
IS_LODs <- IS_LODs[keep, ]

## MVN Model with standard LODs

mvn_priors_LODs <- standard_prior_elicitation("MVN", Z = Z_LODs)
Z_imputed_MVN_LODs <- gibbs_mvn(Z = Z_LODs,
                                mu_0 = mvn_priors_LODs[["mu"]],
                                Omega = mvn_priors_LODs[["Omega"]],
                                v0 = mvn_priors_LODs[["v0"]],
                                S_0 = mvn_priors_LODs[["S0"]],
                                cens = LODs, R = 10)  

## DPML Model with standard LODs 

dpml_priors_LODs <- standard_prior_elicitation("DPML", Z = Z_LODs, Ek = 2) 

Z_imputed_DPML_LODs <- gibbs_DPML(Z = Z_LODs,
                                  mu_0 = dpml_priors_LODs[["mu_G0"]],
                                  Sigma_0 = dpml_priors_LODs[["Sigma_G0"]],
                                  c = dpml_priors_LODs[["c"]],
                                  v0 = dpml_priors_LODs[["v0"]],
                                  S_0 = dpml_priors_LODs[["S0"]],
                                  cens = LODs,
                                  R = 10)


## MVN Model with individual-specific LODs

mvn_priors_IS_LODs <- standard_prior_elicitation("MVN", Z = Z_IS_LODs)
Z_imputed_MVN_IS_LODs <- gibbs_mvn(Z = Z_IS_LODs,
                                   mu_0 = mvn_priors_IS_LODs[["mu"]],
                                   Omega = mvn_priors_IS_LODs[["Omega"]],
                                   v0 = mvn_priors_IS_LODs[["v0"]],
                                   S_0 = mvn_priors_IS_LODs[["S0"]],
                                   cens = IS_LODs, R = 10)

## DPML Model with individual-specific LODs

dpml_priors_IS_LODs <- standard_prior_elicitation("DPML", Z = Z_IS_LODs)

Z_imputed_DPML_IS_LODs <- gibbs_DPML(
  Z = Z_IS_LODs,mu_0 = dpml_priors_IS_LODs[["mu_G0"]],
  Sigma_0 = dpml_priors_IS_LODs[["Sigma_G0"]],
  c = dpml_priors_IS_LODs[["c"]],
  v0 = dpml_priors_IS_LODs[["v0"]],
  S_0 = dpml_priors_IS_LODs[["S0"]],
  cens = IS_LODs, R = 10)


# Basic checks 

check_imputation_LOD <- function(Z_orig = Z_LODs, Z_imp, LOD, tol = .Machine$double.eps^0.5) {
  
  # initial checks
  if (!all(dim(Z_orig) == dim(Z_imp))) {
    stop("Z_orig e Z_imp must have the same dimensions [n x q].")
  }
  n <- nrow(Z_orig); q <- ncol(Z_orig)
  
  # LOD matrix if LOD is a  vector
  if (is.vector(LOD) && length(LOD) == q) {
    lod_mat <- matrix(LOD, nrow = n, ncol = q, byrow = TRUE)
  } else if (is.matrix(LOD) && all(dim(LOD) == c(n, q))) {
    lod_mat <- LOD
  } else {
    stop("LOD must be a vector of length q or a matrix of length n x q")
  }
  
  # Detect original NA positions
  miss_idx <- is.na(Z_orig)
  
  # Check: TRUE if Z_imp > LOD + tol, FALSE otherwise 
  comp <- Z_imp > (lod_mat + tol)
  
  # Consider only imputed positions
  viol_idx <- miss_idx & comp
  
  # Summaries
  total_viol     <- sum(viol_idx,      na.rm = TRUE)
  any_viol       <- total_viol > 0L
  viol_by_col    <- colSums(viol_idx,  na.rm = TRUE)
  
  res <- list(
    any_violation        = any_viol,
    total_violations     = total_viol,
    violations_by_column = viol_by_col,
    violations_matrix    = viol_idx
  )
  
  if (res$any_violation) {
    warning("Found ", res$total_violations, 
            " imputations > LOD (see res$violations_by_column).")
  } else {
    message("All imputations are <= LOD.")
  }
  
  return(res)
}


# General LODs (GL)
check_mvn_GL <- check_imputation_LOD(Z_imp = Z_imputed_MVN_LODs[[1]][, ,sample(1:dim(Z_imputed_MVN_LODs[[1]])[3], 1)],
                 LOD = LODs)

check_DPML_GL <- check_imputation_LOD(Z_imp = Z_imputed_DPML_LODs[[1]][, ,sample(1:dim(Z_imputed_DPML_LODs[[1]])[3], 1)],
                             LOD = LODs)

# Individual-specific LODs

check_mvn_IsL <- check_imputation_LOD(Z_orig = Z_IS_LODs, Z_imp = Z_imputed_MVN_IS_LODs[[1]][, ,sample(1:dim(Z_imputed_MVN_IS_LODs[[1]])[3], 1)],
                                      LOD = IS_LODs)
  
check_DPML_IsL <- check_imputation_LOD(Z_orig = Z_IS_LODs, Z_imp = Z_imputed_DPML_IS_LODs[[1]][, ,sample(1:dim(Z_imputed_DPML_IS_LODs[[1]])[3], 1)],
                                       LOD = IS_LODs)










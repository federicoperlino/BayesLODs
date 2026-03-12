######### Utility functions for Bayesian LOD Imputation #####################
#
# Shared helpers for the MVN and DPML Gibbs samplers.
#
# References
# ----------
# Hoff, P. D. (2009). A First Course in Bayesian Statistical Methods.
# Tierney, L. (1994). Markov chains for exploring posterior distributions.
# Geweke, J. (1991). Efficient simulation from the multivariate normal
#   and Student-t distributions subject to linear constraints.

library(tmvtnorm)

# -------------------------------------------------------------------------
# Basic samplers (Hoff, 2009)
# -------------------------------------------------------------------------

#' Sample from a multivariate normal distribution
#' @param n  Number of draws
#' @param mu Mean vector (length q)
#' @param Sigma Covariance matrix (q x q)
#' @return Matrix of dimension n x q
rmvnorm <- function(n, mu, Sigma) {
  p <- length(mu)
  res <- matrix(0, nrow = n, ncol = p)
  if (n > 0 && p > 0) {
    E <- matrix(rnorm(n * p), n, p)
    res <- t(t(E %*% chol(Sigma)) + c(mu))
  }
  res
}

#' Sample from a Wishart distribution
#' @param n   Number of draws
#' @param nu0 Degrees of freedom
#' @param S0  Scale matrix (q x q)
#' @return Array of dimension q x q x n
rwish <- function(n, nu0, S0) {
  sS0 <- chol(S0)
  S <- array(dim = c(dim(S0), n))
  for (i in 1:n) {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[, , i] <- t(Z) %*% Z
  }
  S[, , 1:n, drop = FALSE]
}

# -------------------------------------------------------------------------
# Helpers for truncated sampling
# -------------------------------------------------------------------------

#' Deterministic feasible start value (used only at initialisation)
#'
#' Returns a point strictly below the censoring limit, used as a fallback
#' when no previous chain state is available (iteration s = 1).
#' @param cens_val Scalar censoring limit
#' @return Scalar below cens_val
fix_start_value <- function(cens_val) {
  if (cens_val >= 0) {
    cens_val / sqrt(2)
  } else {
    cens_val - abs(cens_val) / sqrt(2)
  }
}

#' Single draw from a truncated multivariate normal
#'
#' Thin wrapper around tmvtnorm::rtmvnorm that passes through the
#' algorithm choice and, for the Gibbs algorithm, the warm-start state.
#'
#' @param mean        Conditional mean vector
#' @param sigma       Conditional covariance matrix
#' @param upper       Upper truncation limits (vector, same length as mean)
#' @param ag          Algorithm: "gibbs" (default), "rejection", "gibbsR"
#' @param start.value Warm-start point for algorithm = "gibbs"
#' @param burn.in.samples  Burn-in sweeps before collecting (default 0)
#' @param thinning    Number of sweeps between returned samples (default 1)
#' @return Numeric vector of length = length(mean)
tmvtnorm_draw <- function(mean, sigma, upper, ag = "gibbs",
                          start.value = NULL,
                          burn.in.samples = 0,
                          thinning = 1) {
  args <- list(
    n     = 1,
    mean  = as.numeric(mean),
    sigma = sigma,
    lower = rep(-Inf, length(mean)),
    upper = upper,
    algorithm = ag
  )

  if (ag == "gibbs") {
    args$start.value     <- as.numeric(start.value)
    args$burn.in.samples <- burn.in.samples
    args$thinning        <- thinning
  }

  as.numeric(do.call(tmvtnorm::rtmvnorm, args))
}

# -------------------------------------------------------------------------
# Imputation kernels (warm-started)
# -------------------------------------------------------------------------

#' Impute censored components for the MVN model
#'
#' Given a row of observed data (with NA for censored entries), the current
#' chain state for that row, and the current parameter values, draw the
#' missing components from their truncated MVN full conditional.
#'
#' @param obs_row   Numeric vector: original data row (NA = censored)
#' @param curr_row  Numeric vector: current chain state for the full row
#' @param mu        Current mean vector (length q)
#' @param Sigma     Current covariance matrix (q x q)
#' @param cens_row  Censoring limits for this row (length q)
#' @param ag        Algorithm for tmvtnorm ("gibbs", "rejection", ...)
#' @param burn.in.samples  Burn-in for internal Gibbs (default 0)
#' @param thinning  Thinning for internal Gibbs (default 1)
#' @return Complete row (numeric vector of length q)
tnorm_impute_state <- function(obs_row, curr_row, mu, Sigma, cens_row,
                               ag = "gibbs",
                               burn.in.samples = 0,
                               thinning = 1) {
  miss  <- is.na(obs_row)
  id.na <- which(miss)

  if (length(id.na) == 0) return(obs_row)

  # Build current state: observed components fixed, missing from chain

  z        <- curr_row
  z[!miss] <- obs_row[!miss]

  # Fallback for first iteration (curr_row still contains NA)
  if (any(is.na(z[id.na]))) {
    z[id.na] <- sapply(cens_row[id.na], fix_start_value)
  }

  # Conditional distribution parameters
  Sigma_oo <- Sigma[-id.na, -id.na, drop = FALSE]
  Sigma_mo <- Sigma[ id.na, -id.na, drop = FALSE]
  Sigma_mm <- Sigma[ id.na,  id.na, drop = FALSE]

  Sigma_oo_inv <- solve(Sigma_oo)
  mu.cond    <- mu[id.na] + Sigma_mo %*% Sigma_oo_inv %*% (z[-id.na] - mu[-id.na])
  Sigma.cond <- Sigma_mm  - Sigma_mo %*% Sigma_oo_inv %*% t(Sigma_mo)

  # Warm-start from current chain state of the missing components
  z[id.na] <- tmvtnorm_draw(
    mean          = mu.cond,
    sigma         = Sigma.cond,
    upper         = cens_row[id.na],
    ag            = ag,
    start.value   = z[id.na],
    burn.in.samples = burn.in.samples,
    thinning      = thinning
  )

  z
}

#' Impute censored components for the DPML model
#'
#' Same logic as \code{tnorm_impute_state} but uses the observation-specific
#' location parameter theta_i instead of the global mean mu.
#'
#' @param obs_row    Numeric vector: original data row (NA = censored)
#' @param curr_row   Numeric vector: current chain state for the full row
#' @param theta_row  Current location parameter for this observation
#' @param Sigma      Current shared covariance matrix (q x q)
#' @param cens_row   Censoring limits for this row
#' @param ag         Algorithm for tmvtnorm
#' @param burn.in.samples  Burn-in for internal Gibbs (default 0)
#' @param thinning   Thinning for internal Gibbs (default 1)
#' @return Complete row (numeric vector of length q)
tnorm_mixture_impute_state <- function(obs_row, curr_row, theta_row, Sigma,
                                       cens_row,
                                       ag = "gibbs",
                                       burn.in.samples = 0,
                                       thinning = 1) {
  miss  <- is.na(obs_row)
  id.na <- which(miss)

  if (length(id.na) == 0) return(obs_row)

  z        <- curr_row
  z[!miss] <- obs_row[!miss]

  if (any(is.na(z[id.na]))) {
    z[id.na] <- sapply(cens_row[id.na], fix_start_value)
  }

  Sigma_oo <- Sigma[-id.na, -id.na, drop = FALSE]
  Sigma_mo <- Sigma[ id.na, -id.na, drop = FALSE]
  Sigma_mm <- Sigma[ id.na,  id.na, drop = FALSE]

  Sigma_oo_inv <- solve(Sigma_oo)
  theta.cond  <- theta_row[id.na] +
    Sigma_mo %*% Sigma_oo_inv %*% (z[-id.na] - theta_row[-id.na])
  Sigma.cond  <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% t(Sigma_mo)

  z[id.na] <- tmvtnorm_draw(
    mean          = theta.cond,
    sigma         = Sigma.cond,
    upper         = cens_row[id.na],
    ag            = ag,
    start.value   = z[id.na],
    burn.in.samples = burn.in.samples,
    thinning      = thinning
  )

  z
}

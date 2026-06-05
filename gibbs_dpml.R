# =========================================================================
# DPML Imputation Model 
# =========================================================================
#
# Gibbs sampler for the nonparametric (location Dirichlet process mixture)
# imputation model described in Section 3.2 of the paper.
#
# At each outer iteration s the sampler:
#   1. Updates theta_i via the generalised Polya urn (eq. 9)  [exact]
#   2. Updates Sigma ~ InvWishart(...)                         [exact]
#   3. Imputes y_i1  ~ TN(theta_{i,1|2}, Sigma_{1|2}; zeta)
#      using a componentwise Gibbs kernel warm-started from Y[i, , s-1]
#   4. (Optional) Acceleration step: re-samples unique theta* values
#
# See R/utils.R for the imputation kernel and R/gibbs_mvn.R and 
# R/gibbs_mfm for the analogous samplers.

library(plyr)

# -------------------------------------------------------------------------
# Internal: univariate Normal truncated 
# -------------------------------------------------------------------------
.rtruncnorm_upper <- function(mean, sd, upper) {
  sd <- max(sd, 1e-12)
  p_upper <- stats::pnorm((upper - mean) / sd)
  p_upper <- min(max(p_upper, .Machine$double.xmin), 1 - .Machine$double.eps)
  u <- stats::runif(1, min = 0, max = p_upper)
  mean + sd * stats::qnorm(u)
}

# -------------------------------------------------------------------------
# Internal: Gibbs update for a row with all entries censored / missing
#   Y ~ MVN(theta_row, Sigma) truncated componentwise above by cens_row
#   Warm-started from curr_row
# -------------------------------------------------------------------------
.impute_all_na_row <- function(curr_row, theta_row, Sigma, cens_row,
                               burn.in.samples = 0,
                               thinning = 1, 
                               jitter = 1e-10,
                               max_tries = 6,
                               rel = TRUE) {
  q <- length(theta_row)
  y <- as.numeric(curr_row)
  cens_row <- as.numeric(cens_row)
  
  # Safe warm start
  if (anyNA(y)) {
    miss0 <- is.na(y)
    y[miss0] <- sapply(cens_row[miss0], fix_start_value)
  }
  
  n_sweeps <- burn.in.samples + max(1L, thinning)
  
  for (sweep in seq_len(n_sweeps)) {
    for (k in seq_len(q)) {
      idx <- setdiff(seq_len(q), k)
      
      if (length(idx) == 0) {
        cond_mean <- theta_row[k]
        cond_var  <- Sigma[k, k]
      } else {
        Sigma_22 <- Sigma[idx, idx, drop = FALSE]
        Sigma_12 <- matrix(Sigma[k, idx], nrow = 1)
        Sigma_21 <- matrix(Sigma[idx, k], ncol = 1)
        
        diff_y <- y[idx] - theta_row[idx]
        
        cond_mean <- theta_row[k] +
          as.numeric(Sigma_12 %*% solve_spd(Sigma_22, diff_y, jitter = jitter, max_tries = max_tries, rel = rel))
        
        cond_var <- Sigma[k, k] -
          as.numeric(Sigma_12 %*% solve_spd(Sigma_22, Sigma_21, jitter = jitter, max_tries = max_tries, rel = rel))
      }
      
      cond_var <- max(cond_var, 1e-12)
      
      y[k] <- .rtruncnorm_upper(
        mean  = cond_mean,
        sd    = sqrt(cond_var),
        upper = cens_row[k]
      )
    }
  }
  
  y
}

# -------------------------------------------------------------------------
# Internal: row imputation wrapper
#   - if nothing missing: return observed row
#   - if everything missing: use dedicated Gibbs kernel above
#   - otherwise: use your original tnorm_mixture_impute_state()
# -------------------------------------------------------------------------
.impute_row_dpml <- function(obs_row, curr_row, theta_row, Sigma, cens_row,
                             ag = "gibbs",
                             burn.in.samples = 0,
                             thinning = 1, 
                             jitter = 1e-10,
                             max_tries = 6,
                             rel = TRUE) {
  miss <- is.na(obs_row)
  
  if (!any(miss)) {
    return(as.numeric(obs_row))
  }
  
  if (all(miss)) {
    return(.impute_all_na_row(
      curr_row         = curr_row,
      theta_row        = theta_row,
      Sigma            = Sigma,
      cens_row         = cens_row,
      burn.in.samples  = burn.in.samples,
      thinning         = thinning, 
      jitter           = jitter,
      max_tries        = max_tries,
      rel              = rel
    ))
  }
  
  as.numeric(
    tnorm_mixture_impute_state(
      obs_row          = obs_row,
      curr_row         = curr_row,
      theta_row        = theta_row,
      Sigma            = Sigma,
      cens_row         = cens_row,
      ag               = ag,
      burn.in.samples  = burn.in.samples,
      thinning         = thinning, 
      jitter           = jitter,
      max_tries        = max_tries,
      rel              = rel
    )
  )
}

# -------------------------------------------------------------------------
# Internal: Polya urn update for theta_i
# -------------------------------------------------------------------------

.update_theta_i <- function(i, theta, Y, Sigma, Sigma_0, mu_0, q, n, s,
                            logc, logcn_1, solve_Sigma_0, log_det_Sigma_0,
                            solve_Sigma, log_det_Sigma, 
                            jitter = 1e-10,
                            max_tries = 6,
                            rel = TRUE) {
  
  theta_minus_i <- theta[-i, , drop = FALSE]
  Yi <- Y[i, , s - 1]
  
  COUNT <- plyr::count(as.data.frame(theta_minus_i))
  theta_star <- as.matrix(COUNT[, 1:q, drop = FALSE])
  n_star     <- COUNT$freq
  k_star     <- nrow(COUNT)
  
  log_pi <- numeric(k_star + 1)
  
  # Weight for a new cluster
  Sigma_i <- inv_spd(solve_Sigma_0 + solve_Sigma, 
                     jitter = jitter, max_tries = max_tries, rel = rel)
  mu_i    <- Sigma_i %*% (solve_Sigma %*% Yi + solve_Sigma_0 %*% mu_0)
  
  log_pi[1] <- logc - logcn_1 - (q / 2) * log(2 * pi) +
    0.5 * (log(det(Sigma_i)) - log_det_Sigma - log_det_Sigma_0) -
    0.5 * (t(Yi) %*% solve_Sigma %*% Yi) -
    0.5 * (t(mu_0) %*% solve_Sigma_0 %*% mu_0) +
    0.5 * (t(mu_i) %*% inv_spd(Sigma_i, jitter = jitter, max_tries = max_tries, rel = rel) %*% mu_i)
  
  # Weights for existing clusters
  log_pi[2:(k_star + 1)] <- log(n_star) - logcn_1 +
    sapply(1:k_star, function(j)
      mvtnorm::dmvnorm(Yi, theta_star[j, ], Sigma, log = TRUE))
  
  # Normalise in log-space
  pi_vec <- exp(log_pi - max(log_pi))
  pi_vec <- pi_vec / sum(pi_vec)
  
  idx <- sample.int(k_star + 1, size = 1, prob = pi_vec)
  if (idx == 1) {
    theta[i, ] <- rmvnorm(1, mu_i, Sigma_i)
  } else {
    theta[i, ] <- theta_star[idx - 1, ]
  }
  
  theta[i, ]
}

# -------------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------------

#' Gibbs sampler for the DPML imputation model
#'
#' @param Z         Data matrix (n x q) with NA for censored entries
#' @param mu_0      Base measure mean (length q)
#' @param Sigma_0   Base measure covariance (q x q)
#' @param v0        Prior df for Sigma
#' @param S_0       Prior scale for Sigma (q x q)
#' @param c         DP concentration parameter
#' @param cens      Censoring limits: vector (length q) or matrix (n x q)
#' @param R         Number of Gibbs iterations
#' @param acc       Logical: use the acceleration step? (default TRUE)
#' @param ag        Algorithm for truncated sampling (default "gibbs")
#' @param burn.in.samples  Burn-in for the inner TMV sampler (default 0)
#' @param thinning  Thinning for the inner TMV sampler (default 1)
#' @param jitter Initial diagonal jitter for SPD operations
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return List with components:
#'   \item{Y}{Array n x q x R of imputed latent data}
#'   \item{theta_store}{Array n x q x R of location parameters}
#'   \item{Sigma_store}{Array q x q x R of covariance matrices}
gibbs_DPML <- function(Z, mu_0, Sigma_0, v0, S_0, c, cens, R,
                       acc = TRUE,
                       ag = "gibbs",
                       burn.in.samples = 0,
                       thinning = 1, 
                       jitter = 1e-10,
                       max_tries = 6,
                       rel = TRUE) {
  
  n <- nrow(Z)
  q <- ncol(Z)
  cens_is_matrix <- is.matrix(cens)
  
  Y           <- array(NA_real_, dim = c(n, q, R))
  theta_store <- array(NA_real_, dim = c(n, q, R))
  Sigma_store <- array(NA_real_, dim = c(q, q, R))
  
  # --- Initialisation -----------------------------------------------------
  theta <- matrix(rnorm(n * q), n, q)
  Sigma <- diag(runif(1, 0, 2), q)
  
  Y_curr <- Z
  for (i in seq_len(n)) {
    miss <- is.na(Y_curr[i, ])
    if (any(miss)) {
      cens_i <- if (cens_is_matrix) cens[i, miss] else cens[miss]
      Y_curr[i, miss] <- sapply(cens_i, fix_start_value)
    }
  }
  
  # First imputation (s = 1)
  Y[, , 1] <- t(vapply(seq_len(n), function(i) {
    .impute_row_dpml(
      obs_row   = Z[i, ],
      curr_row  = Y_curr[i, ],
      theta_row = theta[i, ],
      Sigma     = Sigma,
      cens_row  = if (cens_is_matrix) cens[i, ] else cens,
      ag        = ag,
      burn.in.samples = burn.in.samples,
      thinning  = thinning, 
      jitter = jitter, 
      max_tries = max_tries, 
      rel = rel
    )
  }, FUN.VALUE = numeric(q)))
  
  theta_store[, , 1] <- theta
  Sigma_store[, , 1] <- Sigma
  
  # Pre-compute constants
  logc             <- log(c)
  logcn_1          <- log(c + n - 1)
  solve_Sigma_0    <- inv_spd(Sigma_0, jitter = jitter, max_tries = max_tries, rel = rel)
  log_det_Sigma_0  <- log(det(Sigma_0))
  
  # --- Gibbs iterations ---------------------------------------------------
  for (s in 2:R) {
    solve_Sigma   <- inv_spd(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
    log_det_Sigma <- log(det(Sigma))
    
    # 1. Update theta (Polya urn)
    for (i in seq_len(n)) {
      theta[i, ] <- .update_theta_i(
        i, theta, Y, Sigma, Sigma_0, mu_0, q, n, s,
        logc, logcn_1, solve_Sigma_0, log_det_Sigma_0,
        solve_Sigma, log_det_Sigma, 
        jitter = jitter, 
        max_tries = max_tries, 
        rel = rel
      )
    }
    
    # 2. Update Sigma
    S_theta <- t(Y[, , s - 1] - theta) %*% (Y[, , s - 1] - theta)
    S_post <- symmetrize(S_0 + S_theta)
    S_n  <- inv_spd(S_post, jitter = jitter, max_tries = max_tries, rel = rel)
    Sigma <- rinvwish1(
      v0 = v0 + n,
      S0 = S_n,
      jitter = jitter,
      max_tries = max_tries,
      rel = rel
    )
    
    # 3. Imputation step: warm-start from Y[, , s-1]
    Y[, , s] <- t(vapply(seq_len(n), function(i) {
      .impute_row_dpml(
        obs_row          = Z[i, ],
        curr_row         = Y[i, , s - 1],
        theta_row        = theta[i, ],
        Sigma            = Sigma,
        cens_row         = if (cens_is_matrix) cens[i, ] else cens,
        ag               = ag,
        burn.in.samples  = burn.in.samples,
        thinning         = thinning, 
        jitter           = jitter,
        max_tries        = max_tries,
        rel              = rel
      )
    }, FUN.VALUE = numeric(q)))
    
    # 4. Acceleration step (optional)
    if (acc) {
      COUNT      <- plyr::count(as.data.frame(theta))
      theta_star <- as.matrix(COUNT[, 1:q, drop = FALSE])
      q_star     <- nrow(theta_star)
      
      for (j in seq_len(q_star)) {
        thetaj <- theta_star[j, ]
        idx_j  <- which(apply(theta, 1, function(x)
          identical(as.numeric(x), as.numeric(thetaj))))
        if (length(idx_j) == 0) idx_j <- j
        
        ybarj <- if (length(idx_j) > 1) colMeans(Y[idx_j, , s, drop = FALSE])
        else Y[idx_j, , s]
        Cj    <- length(idx_j)
        
        A_n <- inv_spd(Sigma_0, jitter = jitter, max_tries = max_tries, rel = rel) +
          Cj * inv_spd(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
        b_n <- inv_spd(Sigma_0, jitter = jitter, max_tries = max_tries, rel = rel) %*% 
          mu_0 + Cj * inv_spd(Sigma, jitter = jitter, max_tries = max_tries, rel = rel) %*% ybarj
        
        new_theta <- rmvnorm(1, inv_spd(A_n, jitter = jitter, max_tries = max_tries, rel = rel) %*% b_n, 
                             inv_spd(A_n, jitter = jitter, max_tries = max_tries, rel = rel))
        theta[idx_j, ] <- matrix(new_theta,
                                 nrow = length(idx_j), ncol = q,
                                 byrow = TRUE)
      }
    }
    
    theta_store[, , s] <- theta
    Sigma_store[, , s] <- Sigma
    
    if (s %% 100 == 0) cat("DPML iteration", s, "/", R, "\n")
  }
  
  list(Y = Y, theta_store = theta_store, Sigma_store = Sigma_store)
}
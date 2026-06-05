# ======================================================================
# Data-generating process for the toy example:
#   - Y generated either from a MVN distribution or a finite Gaussian mixture
#   - Common covariance matrix across mixture components
#   - LOD applied to the standardized latent variables
# ======================================================================
 
# ----------------------------
# Helpers
# ----------------------------

is_spd <- function(S) {
  if (!is.matrix(S)) return(FALSE)
  if (nrow(S) != ncol(S)) return(FALSE)
  if (!isTRUE(all.equal(S, t(S), tolerance = 1e-10))) return(FALSE)
  out <- try(chol(S), silent = TRUE)
  !inherits(out, "try-error")
}


rmvnorm_mat <- function(n, mean, Sigma) {
  if (!is_spd(Sigma)) {
    stop("Sigma must be symmetric positive definite.")
  }
  
  p <- length(mean)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- Z %*% chol(Sigma)
  sweep(Y, 2, mean, FUN = "+")
}


safe_scale <- function(Y) {
  mu <- colMeans(Y)
  sdv <- apply(Y, 2, sd)
  sdv[sdv <= 0] <- 1
  
  Ys <- scale(Y, center = mu, scale = sdv)
  
  list(
    Ys = Ys,
    center = mu,
    scale = sdv
  )
}


make_lod_from_censoring <- function(Y_std, censoring) {
  q <- ncol(Y_std)
  
  if (length(censoring) == 1L) {
    censoring <- rep(censoring, q)
  }
  if (length(censoring) != q) {
    stop("censoring must have length 1 or q.")
  }
  if (any(censoring < 0 | censoring >= 1)) {
    stop("Each censoring value must lie in [0, 1).")
  }
  
  lod <- rep(-Inf, q)
  
  for (j in seq_len(q)) {
    if (censoring[j] > 0) {
      lod[j] <- as.numeric(
        quantile(Y_std[, j], probs = censoring[j], type = 8)
      )
    }
  }
  
  lod
}

# ----------------------------
# Main DGP
# ----------------------------

simulate_lod_dgp <- function(
    n = 300,
    
    # Latent y_distribution
    y_distribution = c("mixture", "mvn"),
    
    # MVN specification
    mu_raw = rep(0, 5),
    
    # Mixture specification
    pi_mix = c(0.3, 0.7),
    mu_components_raw = rbind(
      comp1 = c(-2.1, -2.0, -1.4, -1.4, -1.2),
      comp2 = c( 0.35, 0.55, 0.45, 0.45, 0.40)
    ),
    
    # Covariance matrix
    Sigma_Y = NULL,
    
    # Censoring applied after standardization
    censoring = c(0.25, 0.25, 0.20, 0.15, 0.10),
    
    seed = NULL
) {
  
  y_distribution <- match.arg(y_distribution)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ------------------------------------------------------------
  # 1. Determine dimension q
  # ------------------------------------------------------------
  
  if (y_distribution == "mvn") {
    
    mu_raw <- as.numeric(mu_raw)
    q <- length(mu_raw)
    K <- 1L
    
  } else {
    
    mu_components_raw <- as.matrix(mu_components_raw)
    
    K <- nrow(mu_components_raw)
    q <- ncol(mu_components_raw)
    
    if (K < 1) {
      stop("mu_components_raw must contain at least one component.")
    }
    
    if (length(pi_mix) != K) {
      stop("pi_mix must have length equal to the number of mixture components.")
    }
    if (any(pi_mix < 0)) {
      stop("pi_mix must contain non-negative probabilities.")
    }
    if (sum(pi_mix) <= 0) {
      stop("pi_mix must have positive sum.")
    }
    
    pi_mix <- pi_mix / sum(pi_mix)
  }
  
  # ------------------------------------------------------------
  # 2. Censoring checks
  # ------------------------------------------------------------
  
  if (length(censoring) == 1L) {
    censoring <- rep(censoring, q)
  }
  if (length(censoring) != q) {
    stop("censoring must have length 1 or q.")
  }
  if (any(censoring < 0 | censoring >= 1)) {
    stop("Each censoring value must lie in [0, 1).")
  }
  
  # ------------------------------------------------------------
  # 3. Covariance matrix
  # ------------------------------------------------------------
  
  if (is.null(Sigma_Y)) {
    stop("Sigma_Y must be provided.")
  }
  
  Sigma_Y <- as.matrix(Sigma_Y)
  
  if (!all(dim(Sigma_Y) == c(q, q))) {
    stop("Sigma_Y must be q x q.")
  }
  if (!is_spd(Sigma_Y)) {
    stop("Sigma_Y must be symmetric positive definite.")
  }
  
  # ------------------------------------------------------------
  # 4. Generate latent complete data
  # ------------------------------------------------------------
  
  if (y_distribution == "mvn") {
    
    Y_raw <- rmvnorm_mat(
      n = n,
      mean = mu_raw,
      Sigma = Sigma_Y
    )
    
    w_true <- rep(1L, n)
    pi_mix_true <- 1
    
    mu_components_raw_true <- matrix(mu_raw, nrow = 1)
    rownames(mu_components_raw_true) <- "mvn"
    colnames(mu_components_raw_true) <- paste0("X", seq_len(q))
    
  } else {
    
    w_true <- sample(
      x = seq_len(K),
      size = n,
      replace = TRUE,
      prob = pi_mix
    )
    
    Y_raw <- matrix(NA_real_, nrow = n, ncol = q)
    
    for (k in seq_len(K)) {
      idx_k <- which(w_true == k)
      n_k <- length(idx_k)
      
      if (n_k > 0) {
        Y_raw[idx_k, ] <- rmvnorm_mat(
          n = n_k,
          mean = mu_components_raw[k, ],
          Sigma = Sigma_Y
        )
      }
    }
    
    pi_mix_true <- pi_mix
    mu_components_raw_true <- mu_components_raw
  }
  
  colnames(Y_raw) <- paste0("X", seq_len(q))
  
  # ------------------------------------------------------------
  # 5. Standardization
  # ------------------------------------------------------------
  
  Y_scal <- safe_scale(Y_raw)
  Y_std_true <- as.matrix(Y_scal$Ys)
  colnames(Y_std_true) <- paste0("X", seq_len(q))
  
  # ------------------------------------------------------------
  # 6. Apply LOD on standardized scale
  # ------------------------------------------------------------
  
  LOD_std <- make_lod_from_censoring(
    Y_std = Y_std_true,
    censoring = censoring
  )
  
  Z_out <- Y_std_true
  
  for (j in seq_len(q)) {
    if (is.finite(LOD_std[j])) {
      Z_out[Z_out[, j] < LOD_std[j], j] <- NA_real_
    }
  }
  
  # ------------------------------------------------------------
  # 7. Output
  # ------------------------------------------------------------
  
  list(
    y_distribution = y_distribution,
    
    Z_out = Z_out,
    Y_true = Y_std_true,
    
    w_true = w_true,
    pi_mix_true = pi_mix_true,
    
    LOD_std = LOD_std,
    censoring = censoring,
    
    center_Y = Y_scal$center,
    scale_Y  = Y_scal$scale,
    
    mu_components_raw_true = mu_components_raw_true,
    Sigma_Y_true = Sigma_Y
  )
}


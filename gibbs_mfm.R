# =========================================================================
# Dynamic MFM Imputation Model 
# =========================================================================
#
# Gibbs sampler for the imputation model based on a finite mixture
# with an unknown number of components. 
# This specification is closely related to the nonparametric mixture formulation
# described in Section 3.2 of Perlino et al. (2025), but replaces the DP mixture
# with an MFM prior on the number of occupied components.
#
# Hierarchical Model:
#   Y | w, mu_k, Sigma ~ MVN
#   Number of mixture components: K ~ p_K 
#
# Dynamic Mixture of Finite Mixtures (MFM) Model:
#   K ~ p_K => K - 1 ~ Poisson(lambda = 2) or K ~ Uniform(K_+, K_max)
#   p_i | K ~ Dirichlet(c/K, ..., c/K), with c = concentration parameter
#   w_i | p_i, K ~ Categorical(p_i)
#   Y_i | w_i=k  ~ MVN(mu_k, Sigma)
#   mu_k  ~ MVN(m0, Omega), with Omega in common for all mu_k
#   Sigma ~ Inv-Wishart(v0, S0^(-1))
#
# MCMC (s -> s+1):
#   (1) w_i | Y, p_i, mu, Sigma       [Gibbs: Multinomial]
#   (2) K_+ from w labels             [deterministic]
#   (3) K | w (K_+, n)                [Gibbs on a grid]
#   (4) p_i | w                       [Gibbs: Dirichlet]
#   (5) mu_k | Y, w, Sigma            [Gibbs: Normal]
#   (6) Sigma | Y, w, mu              [Gibbs: Inverse-Wishart]
#   (7) Y_cens | ...                  [Gibbs: TMVN]
#
# See R/utils.R for the imputation kernel and R/gibbs_mvn.R and 
# R/gibbs_dpml.R for the analogous parametric and non-parametric samplers.
#
# References
# ----------
# Miller, J. W. and Harrison, M. T. (2018). Mixture Models With a Prior on the 
#   Number of Components.
# Frühwirth-Schnatter, S., Malsiner-Walli, G. and Grün, B. (2021). Generalized 
#   mixtures of finite mixtures and telescoping sampling.
# Perlino, F. L., Nipoti, B., Williams, F. L., and Bellavia A. (2025). A Bayesian
#   Parametric and Nonparametric Approach for the Imputation of Multivariate 
#   Left-Censored Data Due to Limit of Detection.

# ------------------------------------------------------------------------------
# MFM block: K update, allocations, p, mu, Sigma
# ------------------------------------------------------------------------------

sample_K_mfm_dynamic <- function(w, K_max, c_total,
                                 prior_K = c("shifted_poisson", "uniform"),
                                 lambda_K = 2) {
  prior_K <- match.arg(prior_K)
  
  if (anyNA(w)) stop("sample_K_mfm_dynamic(): 'w' contains NA values.")
  if (!is.numeric(K_max) || length(K_max) != 1 || K_max < 1)
    stop("sample_K_mfm_dynamic(): 'K_max' must be a positive scalar.")
  if (!is.numeric(c_total) || length(c_total) != 1 || c_total <= 0)
    stop("sample_K_mfm_dynamic(): 'c_total' must be positive.")
  if (!is.numeric(lambda_K) || length(lambda_K) != 1 || lambda_K <= 0)
    stop("sample_K_mfm_dynamic(): 'lambda_K' must be positive.")
  
  w_comp <- match(w, sort(unique(w)))
  K_plus <- length(unique(w_comp))
  if (K_max < K_plus) stop("sample_K_mfm_dynamic(): 'K_max' must be >= K_plus.")
  
  nk <- tabulate(w_comp, nbins = K_plus)
  grid <- K_plus:K_max
  
  log_pk <- switch(
    prior_K,
    shifted_poisson = dpois(grid - 1, lambda = lambda_K, log = TRUE),
    uniform = rep(-log(length(grid)), length(grid))
  )
  
  # log{ K! / (K-K_plus)! }
  log_fall <- lgamma(grid + 1) - lgamma(grid - K_plus + 1)
  
  # sum_k log Gamma(nk + c/K) - log Gamma(c/K)
  log_prod <- vapply(grid, function(K) {
    a <- c_total / K
    sum(lgamma(nk + a) - lgamma(a))
  }, numeric(1))
  
  log_w <- log_pk + log_fall + log_prod
  prob <- softmax_logw(log_w)
  
  if (anyNA(prob) || any(!is.finite(prob)) || abs(sum(prob) - 1) > 1e-8) {
    stop("sample_K_mfm_dynamic(): invalid probabilities for K.")
  }
  
  sample(grid, size = 1, prob = prob)
}

sample_w <- function(Y, p_i, mu, Sigma,
                     jitter = 1e-10, max_tries = 6, rel = TRUE) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  q <- ncol(Y)
  
  p_i <- as.numeric(p_i)
  K <- length(p_i)
  
  if (!is.matrix(mu) || nrow(mu) != K || ncol(mu) != q)
    stop("sample_w(): mu must be [K x q] with K=length(p_i) and q=ncol(Y).")
  
  preS <- mvn_precompute(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
  log_p_i <- log(pmax(p_i, 1e-300))
  
  w <- integer(n)
  for (i in seq_len(n)) {
    yi <- Y[i, ]
    logw <- numeric(K)
    for (k in seq_len(K)) {
      logw[k] <- log_p_i[k] + log_dmvnorm_pre(y = yi, mean = mu[k, ], pre = preS)
    }
    prob <- softmax_logw(logw)
    
    if (anyNA(prob) || any(!is.finite(prob)) || abs(sum(prob) - 1) > 1e-8) {
      stop("sample_w(): invalid allocation probabilities.")
    }
    
    w[i] <- sample.int(K, size = 1, prob = prob)
  }
  
  w
}

update_p_dirichlet <- function(w, K, c_total) {
  nk <- tabulate(w, nbins = K)
  alpha <- nk + c_total / K
  rdirichlet1(alpha)
}

update_mu_all <- function(Y, w, K,
                          m0, Omega, Sigma,
                          jitter = 1e-10, max_tries = 6, rel = TRUE) {
  
  Y <- as.matrix(Y)
  n <- nrow(Y); q <- ncol(Y)
  
  m0 <- as.numeric(m0)
  if (length(m0) != q) stop("update_mu_all(): length(m0) must equal ncol(Y).")
  
  Omega <- symmetrize(as.matrix(Omega))
  if (!all(dim(Omega) == c(q, q))) stop("update_mu_all(): Omega must be [q x q].")
  
  Sigma <- symmetrize(as.matrix(Sigma))
  if (!all(dim(Sigma) == c(q, q))) stop("update_mu_all(): Sigma must be [q x q].")
  
  nk <- tabulate(w, nbins = K)
  
  # Rowsum returns only present groups; reinsert into full K x q
  sumY_occ <- rowsum(Y, w, reorder = FALSE)
  sumY <- matrix(0, nrow = K, ncol = q)
  if (nrow(sumY_occ) > 0) {
    idx <- as.integer(rownames(sumY_occ))
    sumY[idx, ] <- as.matrix(sumY_occ)
  }
  
  Sigma_inv <- inv_spd(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
  Omega_inv <- inv_spd(Omega, jitter = jitter, max_tries = max_tries, rel = rel)
  
  mu <- matrix(0, nrow = K, ncol = q)
  
  for (k in seq_len(K)) {
    if (nk[k] == 0) {
      mu[k, ] <- as.numeric(rmvnorm(1, m0, Omega,
                                    jitter = jitter, max_tries = max_tries, rel = rel))
      next
    }
    
    Qk <- symmetrize(Omega_inv + nk[k] * Sigma_inv)
    tmpQ <- safe_SPD(Qk, jitter = jitter, max_tries = max_tries, rel = rel)
    Vk <- chol2inv(tmpQ$R)
    
    bvec <- as.numeric(Omega_inv %*% m0 + Sigma_inv %*% sumY[k, ])
    mk <- as.numeric(Vk %*% bvec)
    
    preVk <- mvn_precompute(Vk, jitter = jitter, max_tries = max_tries, rel = rel)
    mu[k, ] <- mk + as.numeric(t(preVk$R) %*% rnorm(q))
  }
  
  mu
}

update_Sigma_invwish <- function(Y, w, mu,
                                 v0, S0,
                                 jitter = 1e-10, max_tries = 6, rel = TRUE) {
  
  Y <- as.matrix(Y)
  n <- nrow(Y); q <- ncol(Y)
  
  mu <- as.matrix(mu)
  if (ncol(mu) != q) stop("update_Sigma_invwish(): mu must have ncol = ncol(Y).")
  
  S0 <- symmetrize(as.matrix(S0))
  if (!all(dim(S0) == c(q, q))) stop("update_Sigma_invwish(): S0 must be [q x q].")
  
  M <- mu[w, , drop = FALSE]
  R <- Y - M
  S <- crossprod(R)
  
  nu_post <- v0 + n
  
  S_post <- symmetrize(S0 + S)
  # rinvwish1 samples IW(nu, Sarg^{-1}); to target IW(nu_post, S_post),
  # pass Sarg = S_post^{-1}.
  S_arg <- inv_spd(S_post, jitter = jitter, max_tries = max_tries, rel = rel)
  
  rinvwish1(v0 = nu_post, S0 = S_arg,
            jitter = jitter, max_tries = max_tries, rel = rel)
}

# Reorder components so occupied labels come first (1..K_plus), empties at tail
reorder_by_occupancy <- function(w, p, mu) {
  K <- length(p)
  used <- sort(unique(w))
  K_plus <- length(used)
  
  w_new <- match(w, used)
  empty <- setdiff(seq_len(K), used)
  
  p_new  <- c(p[used],  if (length(empty) > 0) p[empty] else numeric(0))
  mu_new <- rbind(mu[used, , drop = FALSE],
                  if (length(empty) > 0) mu[empty, , drop = FALSE] else mu[0, , drop = FALSE])
  
  list(w = w_new, p = p_new, mu = mu_new, K_plus = K_plus)
}

# ------------------------------------------------------------------------------
# Fixed-K finite mixture option
# ------------------------------------------------------------------------------
# Optional utilities for running the sampler as a finite mixture model with a
# practitioner-specified number of occupied components.
#
# NOTE:
# When this option is used, the model is no longer a dynamic MFM. The number of
# occupied components is fixed by the practitioner, and the sampler should be 
# interpreted as a finite mixture sampler rather than as an MFM sampler.
# ------------------------------------------------------------------------------

sample_w_fixed_Kplus <- function(Y, p_i, mu, Sigma, K_plus_fixed,
                                 max_resample = 1000,
                                 jitter = 1e-10, max_tries = 6, rel = TRUE) {
  
  K <- length(p_i)
  
  if (K != K_plus_fixed) {
    stop("sample_w_fixed_Kplus(): K must equal K_plus_fixed.")
  }
  
  for (attempt in seq_len(max_resample)) {
    w <- sample_w(
      Y = Y,
      p_i = p_i,
      mu = mu,
      Sigma = Sigma,
      jitter = jitter,
      max_tries = max_tries,
      rel = rel
    )
    
    if (length(unique(w)) == K_plus_fixed) {
      return(w)
    }
  }
  
  stop("sample_w_fixed_Kplus(): failed to sample allocations with all clusters occupied.")
}

init_w_occupied <- function(n, K, p = NULL) {
  
  if (n < K) {
    stop("init_w_occupied(): n must be >= K.")
  }
  
  if (is.null(p)) {
    p <- rep(1 / K, K)
  }
  
  w <- c(seq_len(K), sample.int(K, size = n - K, replace = TRUE, prob = p))
  sample(w)
}

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

#' Gibbs sampler for the MFM imputation model
#'
#' If `K_plus_fixed` is provided, the dynamic MFM update of K is disabled and
#' the sampler is run as a fixed-K finite mixture model.
#'
#' @param Z             Data matrix (n x q) with NA for censored entries
#' @param cens          Censoring limits: vector (length q) or matrix (n x q)
#' @param R             Number of Gibbs iterations (default 1000)
#' @param burn          Burn-in period for Gibbs MCMC (default 0)
#' @param thin          Thinning parameter (default 1)
#' @param K_max         Maximum number of clusters allowed (default 10)
#' @param K_init        Warm start for the number of clusters (default 3)
#' @param K_plus_fixed  Optional integer. If provided, the sampler is run as a
#'   finite mixture model with a fixed number of occupied components. In this
#'   case, the dynamic MFM update of K is disabled and the model should not be
#'   interpreted as an MFM (default NULL).
#' @param c_total       Concentration parameter (default 1)
#' @param prior_K       Prior distribution for the number of components (default "shifted_poisson")
#' @param lambda_K      Expected values for the shifted-poisson distribution (default 2)
#' @param m0            Prior mean for mu (length q)
#' @param Omega         Prior Covariance matrix for mu (q x q)
#' @param v0            Prior df for Sigma
#' @param S0            Prior scale for Sigma (q x q)
#' @param ag            Algorithm for truncated sampling (default "gibbs")
#' @param thinning_tmv  Thinning for the inner TMV sampler (default 1)
#' @param jitter        Jitter value for stable algebra (default 1e-10)
#' @param max_tries     Max number of attempts with increasing jitter (default 6)
#' @param rel           Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of a square matrix
#' @param verbose_step  Number of Gibbs iterations between progress messages (default 100).
#' @param verbose       Logical; if TRUE, print the current iteration every `verbose_step` steps (default TRUE)
#' @return A list with posterior draws saved after burn-in and thinning:
#'   \item{K}: integer vector of length n_save containing the sampled number of mixture components
#'   \item{K_plus}: integer vector of length n_save of effective number of components occupied
#'   \item{p}: List of length n_save; each element is a numeric vector of mixture weights of length K^(s) 
#'   \item{w}: List of length n_save; each element is an integer vector of length n with cluster allocations 
#'   \item{mu}: List of length n_save; each element is a matrix of dimension `K^(s) x q` containing the component-specific means
#'   \item{Sigma}: Array q x q x n_save containing posterior draws of the covariance matrix
#'   \item{Y}: Array n x q x n_save of imputed latent data
gibbs_MFM <- function(Z, cens,
                      R = 1000, burn = 0, thin = 1,
                      # MFM controls 
                      K_max = 10, K_init = 3,
                      K_plus_fixed = NULL,
                      c_total = 1,
                      prior_K = c("shifted_poisson", "uniform"),
                      lambda_K = 2,
                      # Mixture priors
                      m0, Omega,
                      v0, S0,
                      # TMVN controls 
                      ag = "gibbs", thinning_tmv = 1, 
                      # Numerics
                      jitter = 1e-10, max_tries = 6, rel = TRUE, 
                      verbose_step = 100,
                      verbose = TRUE) {
  
  # ---- args ----
  prior_K <- match.arg(prior_K)
  
  if (!is.numeric(K_max) || length(K_max) != 1 || abs(K_max - round(K_max)) > 1e-8 || K_max < 1) {
    stop("gibbs_MFM(): K_max must be a positive integer.")
  }
  K_max <- as.integer(round(K_max))
  
  if (!is.numeric(K_init) || length(K_init) != 1 || abs(K_init - round(K_init)) > 1e-8 || K_init < 1) {
    stop("gibbs_MFM(): K_init must be a positive integer.")
  }
  K_init <- as.integer(round(K_init))
  
  if (!is.null(K_plus_fixed)) {
    if (!is.numeric(K_plus_fixed) || length(K_plus_fixed) != 1 ||
        abs(K_plus_fixed - round(K_plus_fixed)) > 1e-8 || K_plus_fixed < 1) {
      stop("gibbs_MFM(): K_plus_fixed must be a positive integer or NULL.")
    }
    
    K_plus_fixed <- as.integer(round(K_plus_fixed))
    
    if (K_plus_fixed > nrow(as.matrix(Z))) {
      stop("gibbs_MFM(): K_plus_fixed cannot be larger than n.")
    }
    
    K_init <- K_plus_fixed
    K_max  <- K_plus_fixed
  }
  
  if (!is.numeric(c_total) || length(c_total) != 1 || !is.finite(c_total) || c_total <= 0) {
    stop("gibbs_MFM(): c_total must be a positive finite scalar.")
  }
  
  if (!is.numeric(lambda_K) || length(lambda_K) != 1 || !is.finite(lambda_K) || lambda_K <= 0) {
    stop("gibbs_MFM(): lambda_K must be a positive finite scalar.")
  }
  
  # ---- data ----
  Y_obs <- as.matrix(Z)
  n <- nrow(Y_obs); q <- ncol(Y_obs)
  
  if (!is.numeric(v0) || length(v0) != 1 || abs(v0 - round(v0)) > 1e-8 || v0 < q + 1L) {
    stop("gibbs_MFM(): v0 must be an integer >= q + 1.")
  }
  v0 <- as.integer(round(v0))
  
  if (burn >= R) stop("burn must be < R.")
  n_save <- floor((R - burn) / thin)
  if (n_save <= 0) stop("n_save computed as 0. Check burn/thin/R.")
  
  # -------------------------
  # 0) Initializations
  # -------------------------
  
  K <- K_init
  if (K > K_max) stop("K_init must satisfy 1 <= K_init <= K_max.")
  
  mu <- rmvnorm(K, mu = m0, Sigma = Omega,
                jitter = jitter, max_tries = max_tries, rel = rel)
  
  S0_arg <- inv_spd(S0, jitter = jitter, max_tries = max_tries, rel = rel)
  Sigma  <- rinvwish1(v0 = v0, S0 = S0_arg,
                      jitter = jitter, max_tries = max_tries, rel = rel)
  
  p <- rdirichlet1(rep(c_total / K, K))
  if (is.null(K_plus_fixed)) {
    w <- sample.int(K, size = n, replace = TRUE, prob = p)
  } else {
    w <- init_w_occupied(n = n, K = K, p = p)
  }
  K_plus <- length(unique(w))
  
  Y_state <- init_Y_state(Y_obs, w, mu, Sigma, cens,
                          ag = ag, thinning = thinning_tmv,
                          jitter = jitter, max_tries = max_tries, rel = rel)
  
  # -------------------------
  # Storage
  # -------------------------
  
  K_draw <- integer(n_save)
  Kplus_draw <- integer(n_save)
  
  p_draw  <- vector("list", n_save)
  w_draw  <- vector("list", n_save)
  mu_draw <- vector("list", n_save)
  Sigma_draw  <- array(NA_real_, dim = c(q, q, n_save))
  Y_full_draw <- array(NA_real_, dim = c(n, q, n_save))
  
  # -------------------------
  # MCMC loop
  # -------------------------
  save_idx <- 0L
  
  for (iter in seq_len(R)) {
    
    if (verbose && (iter %% verbose_step == 0)) {
      model_name <- if (is.null(K_plus_fixed)) "MFM" else "Fixed-K finite mixture"
      cat(model_name, "iteration", iter, "/", R, "\n")
    }
    if (anyNA(Y_state) || any(!is.finite(Y_state))) 
      stop("Y_state has NA/Inf before Y-model updates.")
    

    # (1) w | Y, p, mu, Sigma
    if (is.null(K_plus_fixed)) {
      w <- sample_w(
        Y = Y_state, 
        p_i = p, 
        mu = mu, 
        Sigma = Sigma, 
        jitter = jitter, 
        max_tries = max_tries, 
        rel = rel
      )
    } else {
      w <- sample_w_fixed_Kplus(
        Y = Y_state,
        p_i = p,
        mu = mu,
        Sigma = Sigma,
        K_plus_fixed = K_plus_fixed,
        jitter = jitter,
        max_tries = max_tries,
        rel = rel
      )
    }
    # Reorder so occupied labels come first 
    tmp <- reorder_by_occupancy(w, p, mu)
    w <- tmp$w; p <- tmp$p; mu <- tmp$mu
    
    # (2) Effective number of clusters occupied
    K_plus <- tmp$K_plus 
    
    # (3) K | w
    if (is.null(K_plus_fixed)) {
      
      K_new <- sample_K_mfm_dynamic(
        w, 
        K_max = K_max, 
        c_total = c_total,
        prior_K = prior_K, 
        lambda_K = lambda_K
      )
      if (K_new != K) {
        if (K_new > K) {
          add <- K_new - K
          mu_add <- rmvnorm(add, mu = m0, Sigma = Omega,
                            jitter = jitter, max_tries = max_tries, rel = rel)
          mu <- rbind(mu, mu_add)
          p <- c(p, rep(1 / K_new, add))
          p <- p / sum(p)
        } else {
          mu <- mu[seq_len(K_new), , drop = FALSE]
          p  <- p[seq_len(K_new)]
          p <- p / sum(p)
        }
        K <- K_new
      }
      
    } else {
      K <- K_plus_fixed
      if (K_plus != K_plus_fixed) {
        stop("gibbs_MFM(): K_plus changed despite K_plus_fixed constraint.")
      }
    }
    
    # (4) p | w, K
    p <- update_p_dirichlet(w, K = K, c_total = c_total)

    
    # (5) mu | Y, w, Sigma
    mu <- update_mu_all(Y_state, w, K,
                        m0 = m0, Omega = Omega, Sigma = Sigma,
                        jitter = jitter, max_tries = max_tries, rel = rel)
    
    # (6) Sigma | Y, w, mu
    Sigma <- update_Sigma_invwish(Y_state, w, mu,
                                  v0 = v0, S0 = S0,
                                  jitter = jitter, max_tries = max_tries, rel = rel)
    
    if (any(!is.finite(Sigma))) {
      stop(paste("Sigma has non-finite values at iter", iter))
    }
    
    tmp_check_Sigma <- tryCatch(
      safe_SPD(Sigma, jitter = jitter, max_tries = max_tries, rel = rel),
      error = function(e) NULL
    )
    
    if (is.null(tmp_check_Sigma)) {
      stop(paste("Sigma is not SPD at iter", iter))
    }
    
    # (7) Y_cens | w, mu, Sigma
    Y_state <- update_Y_cens(
      Y_state = Y_state,
      Y_obs = Y_obs,
      w = w,
      mu = mu,
      Sigma = Sigma,
      cens = cens,
      ag = ag,
      thinning = thinning_tmv,
      burn_in = 0,
      jitter = jitter,
      max_tries = max_tries,
      rel = rel
    )
    
    if (any(!is.finite(Y_state))) {
      stop(paste("Y_state has non-finite values at iter", iter))
    }
    
    idx_obs <- which(!is.na(Y_obs), arr.ind = TRUE)
    if (nrow(idx_obs) > 0) {
      if (any(abs(Y_state[idx_obs] - Y_obs[idx_obs]) > 1e-12)) {
        stop(paste("Observed Y changed during update at iter", iter))
      }
    }
    
    idx_cens <- which(is.na(Y_obs), arr.ind = TRUE)
    if (nrow(idx_cens) > 0) {
      lod_check <- as_cens_matrix(cens, n = n, q = q)
      if (any(Y_state[idx_cens] >= lod_check[idx_cens])) {
        stop(paste("Censored Y not below cens at iter", iter))
      }
    }
    
    # Save draws
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      save_idx <- save_idx + 1L
      
      K_draw[save_idx] <- K
      Kplus_draw[save_idx] <- K_plus
      
      p_draw[[save_idx]]  <- p
      w_draw[[save_idx]]  <- w
      mu_draw[[save_idx]] <- mu
      Sigma_draw[, , save_idx]  <- Sigma
      Y_full_draw[, , save_idx] <- Y_state
    }
  }
    
  list(
      K = K_draw,
      K_plus = Kplus_draw,
      p = p_draw,
      w = w_draw,
      mu = mu_draw,
      Sigma = Sigma_draw,
      Y = Y_full_draw
      )
}


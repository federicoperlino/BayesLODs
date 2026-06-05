# =========================================================================
# Utility functions for Bayesian LOD Imputation 
# =========================================================================
#
# Shared helpers for the MVN, MFM and DPML Gibbs samplers.
#
# References
# ----------
# Hoff, P. D. (2009). A First Course in Bayesian Statistical Methods.
# Tierney, L. (1994). Markov chains for exploring posterior distributions.
# Geweke, J. (1991). Efficient simulation from the multivariate normal
#   and Student-t distributions subject to linear constraints.

library(tmvtnorm)

# ------------------------------------------------------------------------------
# Linear algebra utilities
# ------------------------------------------------------------------------------

#' Symmetrize a square matrix
#' 
#' @param A Matrix of dimension q x q
#' @return Symmetric matrix of dimension q x q
symmetrize <- function(A) 0.5 * (A + t(A))

#' Cholesky decomposition of a square matrix, with diagonal jitter if needed
#' 
#' @param A Matrix of dimension q x q
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `A`
#' @return A list with components:
#'   \item{R}: Upper triangular Cholesky factor of the adjusted matrix
#'   \item{eps}: Diagonal jitter added to `A` to obtain a successful decomposition
safe_SPD <- function(A, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  A <- symmetrize(as.matrix(A))
  p <- nrow(A)
  if (p != ncol(A)) stop("safe_SPD(): A must be square.")
  if (p == 0) stop("safe_SPD(): empty matrix.")
  if (any(!is.finite(A))) stop("safe_SPD(): A contains non-finite values.")
  
  scale_A <- 1
  if (rel) {
    d <- diag(A)
    if (any(!is.finite(d))) stop("safe_SPD(): non-finite diagonal.")
    scale_A <- max(1, mean(abs(d)))
  }
  
  out <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(out)) return(list(R = out, eps = 0))
  
  I <- diag(p)
  eps <- jitter * scale_A
  for (i in seq_len(max_tries)) {
    out <- tryCatch(chol(A + eps * I), error = function(e) NULL)
    if (!is.null(out)) return(list(R = out, eps = eps))
    eps <- eps * 10
  }
  stop("safe_SPD(): matrix not SPD even after jitter.")
}

#' Safe inverse of a square matrix via Cholesky decomposition
#'
#' @param A Matrix of dimension q x q
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `A`
#' @return Inverse matrix obtained via Cholesky decomposition, possibly after
#'   adding a small diagonal perturbation to `A`
inv_spd <- function(A, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  tmp <- safe_SPD(A, jitter = jitter, max_tries = max_tries, rel = rel)
  chol2inv(tmp$R)
}

#' Return an SPD version of a square matrix and its Cholesky factor
#' 
#' @param A Matrix of dimension q x q
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `A`
#' @return A list with components:
#'   \item{A}: Symmetric positive definite matrix obtained from `A`, possibly after
#'   adding a diagonal perturbation.
#'   \item{R}: Upper triangular Cholesky factor of the adjusted matrix
#'   \item{eps}: Diagonal jitter added to `A`
make_spd <- function(A, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  A <- symmetrize(as.matrix(A))
  tmp <- safe_SPD(A, jitter = jitter, max_tries = max_tries, rel = rel)
  M <- crossprod(tmp$R)  # = A + eps*I
  list(A = M, R = tmp$R, eps = tmp$eps)
}

#' Solve a linear system with an SPD matrix via Cholesky decomposition
#' 
#' @param A Matrix of dimension q x q
#' @param B Right-hand side matrix or vector
#' @param jitter Positive real number used as initial diagonal perturbation.
#' @param max_tries Maximum number of attempts with increasing jitter.
#' @param rel Logical; if TRUE, scale the jitter according to the magnitude
#'   of the diagonal entries of `A`
#' @return Solution of the linear system `A %*% X = B`, possibly after
#'   adding a small diagonal perturbation to `A`
solve_spd <- function(A, B, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  A <- symmetrize(as.matrix(A))
  B <- as.matrix(B)
  tmp <- safe_SPD(A, jitter = jitter, max_tries = max_tries, rel = rel)
  R <- tmp$R
  backsolve(R, forwardsolve(t(R), B))
}

#' Compute normalized weights from log-weights using a stable softmax
#' 
#' @param logw Numeric vector of log-weights
#' @return Numeric vector of normalized probabilities summing to 1
softmax_logw <- function(logw) {
  logw <- as.numeric(logw)
  if (length(logw) == 0) return(numeric(0))
  logw[is.na(logw)] <- -Inf
  
  m <- max(logw)
  if (!is.finite(m)) {
    stop("softmax_logw(): all log-weights are -Inf.")
  }
  
  w <- exp(logw - m)
  s <- sum(w)
  
  if (!is.finite(s) || s <= 0) {
    stop("softmax_logw(): invalid normalized weights.")
  }
  
  out <- w / s
  
  if (anyNA(out) || any(!is.finite(out)) || abs(sum(out) - 1) > 1e-8) {
    stop("softmax_logw(): probabilities are invalid after normalization.")
  }
  
  out
}


# ------------------------------------------------------------------------------
# Basic samplers (Hoff, 2009)
# ------------------------------------------------------------------------------

#' Precompute quantities for a multivariate normal covariance matrix
#' 
#' @param Sigma Covariance matrix of dimension q x q.
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `Sigma`
#' @return A list with components:
#'   \item{R}: Upper triangular Cholesky factor of the adjusted covariance matrix.
#'   \item{logdet}: Log-determinant of the adjusted covariance matrix.
#'   \item{Sigma}: Symmetric positive definite covariance matrix, possibly after
#'   adding a diagonal perturbation.
#'   \item{eps}: Diagonal jitter added to `Sigma`.
#' }
mvn_precompute <- function(Sigma, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  spd <- make_spd(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
  list(
    R = spd$R,
    logdet = 2 * sum(log(diag(spd$R))),
    Sigma = spd$A,
    eps = spd$eps
  )
}

#' Sample from a multivariate normal distribution
#' 
#' @param n  Number of draws
#' @param mu Mean vector (length q)
#' @param Sigma Covariance matrix (q x q)
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `Sigma`.
#' @return Matrix of dimension n x q
rmvnorm <- function(n, mu, Sigma, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  mu <- as.numeric(mu)
  p <- length(mu)
  
  if (n <= 0) return(matrix(0, nrow = 0, ncol = p))
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(p, p))) stop("rmvnorm(): dim(Sigma) must be [p x p].")
  
  pre <- mvn_precompute(Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% pre$R
  sweep(X, 2, mu, `+`)
}

#' Log-density of a multivariate normal distribution using Cholesky precomputation
#' 
#' @param y Numeric vector of length q
#' @param mean Mean vector of length q
#' @param pre List of precomputed quantities returned by `mvn_precompute()`
#' @return Log-density of the multivariate normal distribution evaluated at `y`,
#'   including the normalizing constant
log_dmvnorm_pre <- function(y, mean, pre) {
  y <- as.numeric(y)
  mean <- as.numeric(mean)
  p <- length(mean)
  
  z <- forwardsolve(t(pre$R), y - mean)
  quad <- sum(z * z)
  -0.5 * (p * log(2 * pi) + pre$logdet + quad)
}

#' Draw Wishart random matrices using a Gaussian construction
#' 
#' @param n Number of Wishart random matrices to generate
#' @param v0 Degrees of freedom of the Wishart distribution
#' @param S0 Scale matrix of dimension q x q
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `S0`
#' @return A 3-dimensional array of dimension q x q x n containing the sampled
#'   Wishart random matrices. If `n = 0`, an empty array is returned
rwish <- function(n, v0, S0, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  S0 <- as.matrix(S0)
  q  <- nrow(S0)
  
  if (!all(dim(S0) == c(q, q))) stop("rwish(): S0 must be square.")
  if (!is.numeric(v0) || length(v0) != 1 || abs(v0 - round(v0)) > 1e-8 || v0 < q) {
    stop("rwish(): v0 must be an integer >= q.")
  }
  if (!is.numeric(n) || length(n) != 1 || abs(n - round(n)) > 1e-8 || n < 0) {
    stop("rwish(): n must be a nonnegative integer.")
  }
  
  v0 <- as.integer(round(v0))
  n   <- as.integer(round(n))
  
  if (n == 0L) return(array(0, dim = c(q, q, 0)))
  
  pre <- mvn_precompute(S0, jitter = jitter, max_tries = max_tries, rel = rel)
  out <- array(0, dim = c(q, q, n))
  
  for (i in seq_len(n)) {
    Z  <- matrix(rnorm(v0 * q), nrow = v0, ncol = q) %*% pre$R
    Wi <- crossprod(Z)
    out[, , i] <- symmetrize(Wi)
  }
  
  out
}

#' Draw a random matrix from an Inverse-Wishart distribution
#' 
#' @param v0 Degrees of freedom
#' @param S0 Scale matrix of dimension q x q
#' @param jitter Positive real number used as initial diagonal perturbation (default = 1e-10)
#' @param max_tries Maximum number of attempts with increasing jitter (default = 6)
#' @param rel Logical; if TRUE (default), scale the jitter according to the magnitude
#'   of the diagonal entries of `S0`
#' @return A random matrix drawn from an Inverse-Wishart distribution,
#'   obtained by sampling `W ~ Wishart(v0, S0)` and returning `W^{-1}`
rinvwish1 <- function(v0, S0, jitter = 1e-10, max_tries = 6, rel = TRUE) {
  S0 <- symmetrize(as.matrix(S0))
  q <- nrow(S0)
  
  if (!all(dim(S0) == c(q, q))) {
    stop("rinvwish1(): S0 must be square.")
  }
  if (!is.numeric(v0) || length(v0) != 1 || abs(v0 - round(v0)) > 1e-8 || v0 < q + 1L) {
    stop("rinvwish1(): v0 must be an integer >= q + 1.")
  }
  
  v0 <- as.integer(round(v0))
  
  W_arr <- rwish(
    n = 1L,
    v0 = v0,
    S0 = S0,
    jitter = jitter,
    max_tries = max_tries,
    rel = rel
  )
  
  W <- W_arr[, , 1, drop = FALSE]
  W <- matrix(W, nrow = q, ncol = q)
  
  inv_spd(W, jitter = jitter, max_tries = max_tries, rel = rel)
}


#' Draw one sample from a Dirichlet distribution via Gamma variables
#'
#' @param alpha Numeric vector of positive concentration parameters.
#' @return Numeric vector of the same length as `alpha`, containing a single
#'   draw from a Dirichlet distribution and summing to 1.
rdirichlet1 <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (any(!is.finite(alpha)) || any(alpha <= 0))
    stop("rdirichlet1(): alpha must be positive and finite.")
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}
# ------------------------------------------------------------------------------
# Helpers for truncated sampling
# ------------------------------------------------------------------------------

#' Deterministic feasible start value (used only at initialisation)
#'
#' Returns a point strictly below the censoring limit, used as a fallback
#' when no previous chain state is available.
#'
#' @param cens_val Scalar censoring limit
#' @return Scalar below cens_val
fix_start_value <- function(cens_val) {
  if (cens_val >= 0) {
    cens_val / sqrt(2)
  } else {
    cens_val - abs(cens_val) / sqrt(2)
  }
}

#' Coerce censoring limits to an n x q matrix
#'
#' @param cens Either a vector of length q or a matrix n x q
#' @param n Number of observations
#' @param q Number of variables
#' @return Matrix n x q of censoring limits
as_cens_matrix <- function(cens, n, q) {
  if (is.vector(cens) && length(cens) == q) {
    return(matrix(cens, nrow = n, ncol = q, byrow = TRUE))
  }
  if (is.matrix(cens) && all(dim(cens) == c(n, q))) {
    return(cens)
  }
  stop("cens must be a vector of length q or a matrix [n x q].")
}

#' Single draw from a truncated multivariate normal
#'
#' Thin wrapper around tmvtnorm::rtmvnorm that passes through the
#' algorithm choice and, for the Gibbs algorithm, the warm-start state.
#'
#' @param mean Conditional mean vector
#' @param sigma Conditional covariance matrix
#' @param upper Upper truncation limits
#' @param ag Algorithm: "gibbs", "rejection", "gibbsR"
#' @param start.value Warm-start point for algorithm = "gibbs"
#' @param burn.in.samples Burn-in sweeps before collecting
#' @param thinning Number of sweeps between returned samples
#' @return Numeric vector of length = length(mean)
tmvtnorm_draw <- function(mean, sigma, upper, ag = "gibbs",
                          start.value = NULL,
                          burn.in.samples = 0,
                          thinning = 1) {
  args <- list(
    n = 1,
    mean = as.numeric(mean),
    sigma = sigma,
    lower = rep(-Inf, length(mean)),
    upper = as.numeric(upper),
    algorithm = ag
  )
  
  if (ag == "gibbs") {
    args$start.value <- as.numeric(start.value)
    args$burn.in.samples <- burn.in.samples
    args$thinning <- thinning
  }
  
  as.numeric(do.call(tmvtnorm::rtmvnorm, args))
}

# ------------------------------------------------------------------------------
# Conditional Gaussian block
# ------------------------------------------------------------------------------

#' Conditional distribution of missing components under MVN(mu, Sigma)
#'
#' @param z Numeric vector of length q containing observed values and NA in the
#'   censored positions
#' @param mu Mean vector of length q
#' @param Sigma Covariance matrix q x q
#' @param id_miss Integer vector of indices of censored components
#' @param jitter Initial diagonal jitter for SPD solves
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return A list with conditional mean and conditional covariance
mvn_conditional <- function(z, mu, Sigma, id_miss,
                            jitter = 1e-10, max_tries = 6, rel = TRUE) {
  z <- as.numeric(z)
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  
  p <- length(z)
  if (length(mu) != p) stop("mvn_conditional(): length(mu) != length(z).")
  if (!all(dim(Sigma) == c(p, p))) stop("mvn_conditional(): dim(Sigma) must be [p x p].")
  
  id_miss <- sort(unique(id_miss))
  if (length(id_miss) == 0) {
    return(list(mu = numeric(0), Sigma = matrix(0, 0, 0), id_miss = id_miss))
  }
  
  id_obs <- setdiff(seq_len(p), id_miss)
  Sigma <- symmetrize(Sigma)
  
  if (length(id_obs) == 0) {
    Sigma_mm <- symmetrize(Sigma[id_miss, id_miss, drop = FALSE])
    return(list(
      mu = as.numeric(mu[id_miss]),
      Sigma = Sigma_mm,
      id_miss = id_miss
    ))
  }
  
  if (any(is.na(z[id_obs]))) stop("mvn_conditional(): NA present outside id_miss.")
  
  Sigma_oo <- Sigma[id_obs,  id_obs,  drop = FALSE]
  Sigma_mo <- Sigma[id_miss, id_obs,  drop = FALSE]
  Sigma_om <- Sigma[id_obs,  id_miss, drop = FALSE]
  Sigma_mm <- Sigma[id_miss, id_miss, drop = FALSE]
  
  diff_o <- matrix(z[id_obs] - mu[id_obs], ncol = 1)
  sol1 <- solve_spd(Sigma_oo, diff_o, jitter = jitter, max_tries = max_tries, rel = rel)
  mu_cond <- mu[id_miss] + as.numeric(Sigma_mo %*% sol1)
  
  sol2 <- solve_spd(Sigma_oo, Sigma_om, jitter = jitter, max_tries = max_tries, rel = rel)
  Sigma_cond <- Sigma_mm - Sigma_mo %*% sol2
  Sigma_cond <- symmetrize(Sigma_cond)
  
  list(mu = as.numeric(mu_cond), Sigma = Sigma_cond, id_miss = id_miss)
}

# ------------------------------------------------------------------------------
# Row-wise imputation kernel
# ------------------------------------------------------------------------------

#' Impute censored components of one row under a truncated MVN conditional
#'
#' @param z Numeric row vector with NA in censored positions
#' @param mu Mean/location vector for this row
#' @param Sigma Shared covariance matrix
#' @param cens_row Censoring limits for this row
#' @param start_current Optional current state for warm-starting
#' @param ag Algorithm passed to tmvtnorm
#' @param burn_in Burn-in for internal TMVN Gibbs sampler
#' @param thinning Thinning for internal TMVN Gibbs sampler
#' @param jitter Initial diagonal jitter for SPD operations
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return Completed row vector
tnorm_impute <- function(z, mu, Sigma, cens_row,
                         start_current = NULL,
                         ag = "gibbs",
                         burn_in = 1000,
                         thinning = 20,
                         jitter = 1e-10, max_tries = 6, rel = TRUE) {
  
  id.na <- which(is.na(z))
  if (length(id.na) == 0) return(as.numeric(z))
  
  cond <- mvn_conditional(
    z = z,
    mu = mu,
    Sigma = Sigma,
    id_miss = id.na,
    jitter = jitter,
    max_tries = max_tries,
    rel = rel
  )
  
  upper <- as.numeric(cens_row[id.na])
  lower <- rep(-Inf, length(id.na))
  
  if (any(!is.finite(upper))) {
    stop("tnorm_impute(): non-finite upper bounds in cens.")
  }
  
  epsb <- 1e-12
  
  if (!is.null(start_current)) {
    sc <- as.numeric(start_current)
    if (length(sc) == length(z)) sc <- sc[id.na]
    if (length(sc) != length(id.na)) {
      stop("tnorm_impute(): start_current must be length(z) or length(id.na).")
    }
    start.val <- sc
  } else {
    start.val <- vapply(upper, fix_start_value, numeric(1))
  }
  start.val <- pmin(start.val, upper - epsb)
  
  dmiss <- length(id.na)
  
  if (dmiss == 1L) {
    m  <- cond$mu[1]
    v  <- as.numeric(cond$Sigma[1, 1])
    v  <- max(v, 1e-16)
    sd <- sqrt(v)
    
    b <- (upper[1] - m) / sd
    Fb <- pnorm(b)
    Fb <- max(min(Fb, 1 - 1e-16), 1e-16)
    
    u <- runif(1, min = 0, max = Fb)
    x <- m + sd * qnorm(u)
    
    x <- min(x, upper[1] - epsb)
    if (!is.finite(x)) x <- start.val[1]
    
    z[id.na] <- x
    return(as.numeric(z))
  }
  
  pre_cond <- make_spd(cond$Sigma, jitter = jitter, max_tries = max_tries, rel = rel)
  
  draw <- tmvtnorm::rtmvnorm(
    n = 1,
    mean = cond$mu,
    sigma = pre_cond$A,
    lower = lower,
    upper = upper,
    algorithm = ag,
    start.value = start.val,
    burn.in.samples = burn_in,
    thinning = thinning
  )
  draw <- as.numeric(draw)
  
  if (any(!is.finite(draw))) {
    z[id.na] <- pmin(start.val, upper - epsb)
    return(as.numeric(z))
  }
  
  draw <- pmin(draw, upper - epsb)
  z[id.na] <- draw
  as.numeric(z)
}

# ------------------------------------------------------------------------------
# Row-wise wrappers
# ------------------------------------------------------------------------------

#' MVN row-wise imputation wrapper
#'
#' @param obs_row Original row with NA for censored entries
#' @param curr_row Current chain state for the row
#' @param mu Current global mean vector
#' @param Sigma Current covariance matrix
#' @param cens_row Censoring limits for the row
#' @param ag Algorithm for truncated sampling
#' @param burn.in.samples Burn-in for internal Gibbs TMVN
#' @param thinning Thinning for internal Gibbs TMVN
#' @return Completed row vector
tnorm_impute_state <- function(obs_row, curr_row, mu, Sigma, cens_row,
                               ag = "gibbs",
                               burn.in.samples = 0,
                               thinning = 1,
                               jitter = 1e-10, max_tries = 6, rel = TRUE) {
  z <- as.numeric(obs_row)
  miss <- is.na(z)
  
  if (!any(miss)) return(z)
  
  start_row <- as.numeric(curr_row)
  start_row[!miss] <- z[!miss]
  
  if (any(is.na(start_row[miss]))) {
    start_row[miss] <- vapply(cens_row[miss], fix_start_value, numeric(1))
  }
  
  tnorm_impute(
    z = z,
    mu = mu,
    Sigma = Sigma,
    cens_row = cens_row,
    start_current = start_row,
    ag = ag,
    burn_in = burn.in.samples,
    thinning = thinning,
    jitter = jitter,
    max_tries = max_tries,
    rel = rel
  )
}

#' Mixture row-wise imputation wrapper
#'
#' @param obs_row Original row with NA for censored entries
#' @param curr_row Current chain state for the row
#' @param theta_row Current row-specific location vector
#' @param Sigma Shared covariance matrix
#' @param cens_row Censoring limits for the row
#' @param ag Algorithm for truncated sampling
#' @param burn.in.samples Burn-in for internal Gibbs TMVN
#' @param thinning Thinning for internal Gibbs TMVN
#' @param jitter Initial diagonal jitter for SPD operations
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return Completed row vector
tnorm_mixture_impute_state <- function(obs_row, curr_row, theta_row, Sigma,
                                       cens_row,
                                       ag = "gibbs",
                                       burn.in.samples = 0,
                                       thinning = 1,
                                       jitter = 1e-10, max_tries = 6, rel = TRUE) {
  z <- as.numeric(obs_row)
  miss <- is.na(z)
  
  if (!any(miss)) return(z)
  
  start_row <- as.numeric(curr_row)
  start_row[!miss] <- z[!miss]
  
  if (any(is.na(start_row[miss]))) {
    start_row[miss] <- vapply(cens_row[miss], fix_start_value, numeric(1))
  }
  
  tnorm_impute(
    z = z,
    mu = theta_row,
    Sigma = Sigma,
    cens_row = cens_row,
    start_current = start_row,
    ag = ag,
    burn_in = burn.in.samples,
    thinning = thinning,
    jitter = jitter,
    max_tries = max_tries,
    rel = rel
  )
}

# ------------------------------------------------------------------------------
# State update helpers
# ------------------------------------------------------------------------------

#' Update all censored entries row by row under a mixture model
#'
#' @param Y_state Current latent state matrix n x q
#' @param Y_obs Observed matrix with NA in censored positions
#' @param w Allocation vector of length n
#' @param mu Matrix of component means K x q
#' @param Sigma Shared covariance matrix q x q
#' @param cens Censoring limits: vector length q or matrix n x q
#' @param ag Algorithm for truncated sampling
#' @param thinning Thinning for internal TMVN Gibbs sampler
#' @param burn_in Burn-in for internal TMVN Gibbs sampler
#' @param jitter Initial diagonal jitter for SPD operations
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return Updated latent state matrix
update_Y_cens <- function(Y_state, Y_obs, w, mu, Sigma, cens,
                          ag = "gibbs", thinning = 1, burn_in = 0,
                          jitter = 1e-10, max_tries = 6, rel = TRUE) {
  Y_state <- as.matrix(Y_state)
  Y_obs   <- as.matrix(Y_obs)
  
  n <- nrow(Y_state)
  q <- ncol(Y_state)
  
  if (!all(dim(Y_obs) == c(n, q))) {
    stop("update_Y_cens(): Y_state and Y_obs must have the same dimensions.")
  }
  if (length(w) != n) stop("update_Y_cens(): length(w) must equal nrow(Y_state).")
  if (!is.matrix(mu) || ncol(mu) != q) stop("update_Y_cens(): mu must be [K x q].")
  if (any(w < 1) || any(w > nrow(mu))) stop("update_Y_cens(): w must be in 1:nrow(mu).")
  if (any(dim(Sigma) != c(q, q))) stop("update_Y_cens(): Sigma must be [q x q].")
  
  cens_mat  <- as_cens_matrix(cens, n = n, q = q)
  miss_mat <- is.na(Y_obs)
  
  for (i in seq_len(n)) {
    if (any(miss_mat[i, ])) {
      k <- w[i]
      Y_state[i, ] <- tnorm_mixture_impute_state(
        obs_row = Y_obs[i, ],
        curr_row = Y_state[i, ],
        theta_row = mu[k, ],
        Sigma = Sigma,
        cens_row = cens_mat[i, ],
        ag = ag,
        burn.in.samples = burn_in,
        thinning = thinning,
        jitter = jitter,
        max_tries = max_tries,
        rel = rel
      )
    } else {
      Y_state[i, ] <- Y_obs[i, ]
    }
  }
  
  Y_state
}

#' Initialise the latent state row by row
#'
#' @param Y_obs Observed matrix with NA in censored positions
#' @param w Allocation vector of length n
#' @param mu Matrix of component means K x q
#' @param Sigma Shared covariance matrix q x q
#' @param cens Censoring limits: vector length q or matrix n x q
#' @param ag Algorithm for truncated sampling
#' @param thinning Thinning for internal TMVN Gibbs sampler
#' @param burn_in Burn-in for internal TMVN Gibbs sampler
#' @param jitter Initial diagonal jitter for SPD operations
#' @param max_tries Maximum number of jitter expansions
#' @param rel Logical; if TRUE, scale jitter relative to matrix magnitude
#' @return Initialised latent state matrix
init_Y_state <- function(Y_obs, w, mu, Sigma, cens,
                         ag = "gibbs", thinning = 20, burn_in = 0,
                         jitter = 1e-10, max_tries = 6, rel = TRUE) {
  Y_state <- as.matrix(Y_obs)
  update_Y_cens(
    Y_state = Y_state,
    Y_obs = Y_obs,
    w = w,
    mu = mu,
    Sigma = Sigma,
    cens = cens,
    ag = ag,
    thinning = thinning,
    burn_in = burn_in,
    jitter = jitter,
    max_tries = max_tries,
    rel = rel
  )
}

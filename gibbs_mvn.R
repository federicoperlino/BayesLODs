######### MVN Imputation Model ##############################################
#
# Gibbs sampler for the parametric (multivariate normal) imputation model
# described in Section 3.1 of the paper.
#
# At each outer iteration s the sampler:
#   1. Updates mu     ~ N(...)            [exact draw]
#   2. Updates Sigma  ~ InvWishart(...)   [exact draw]
#   3. Imputes y_i1   ~ TN(mu_{1|2}, Sigma_{1|2}; zeta)
#      using a componentwise Gibbs kernel warm-started from Y[i, , s-1].
#
# The warm-start ensures the inner Gibbs kernel is a valid invariant
# transition for the truncated full conditional, so the overall chain
# targets the correct posterior pi(Y, mu, Sigma | Z).
# See R/utils.R for the imputation kernel implementation.

# -------------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------------

#' Gibbs sampler for the MVN imputation model
#'
#' @param Z       Data matrix (n x q) with NA for censored entries
#' @param mu_0    Prior mean for mu (length q)
#' @param Omega   Prior covariance for mu (q x q)
#' @param v0      Prior degrees of freedom for Sigma
#' @param S_0     Prior scale matrix for Sigma (q x q)
#' @param cens    Censoring limits: either a vector of length q (common LODs)
#'                or a matrix n x q (individual-specific LODs)
#' @param R       Number of Gibbs iterations
#' @param ag      Algorithm for truncated sampling ("gibbs", "rejection", ...)
#' @param burn.in.samples  Burn-in sweeps for the inner TMV sampler (default 0)
#' @param thinning         Thinning for the inner TMV sampler (default 1)
#' @return List with components:
#'   \item{Y}{Array n x q x R of imputed latent data}
#'   \item{mu}{Final mean vector}
#'   \item{Sigma}{Final covariance matrix}
gibbs_mvn <- function(Z, mu_0, Omega, v0, S_0, cens, R,
                      ag = "gibbs",
                      burn.in.samples = 0,
                      thinning = 1) {
  q <- ncol(Z)
  n <- nrow(Z)

  Y <- array(NA_real_, dim = c(n, q, R))

  # --- Initialisation -----------------------------------------------------
  mu    <- rnorm(q)
  Sigma <- diag(runif(1, 0, 2), q)

  # Build initial current state (fill NAs with feasible values)
  Y_curr <- Z
  cens_is_matrix <- is.matrix(cens)

  for (i in 1:n) {
    miss <- is.na(Y_curr[i, ])
    if (any(miss)) {
      cens_i <- if (cens_is_matrix) cens[i, miss] else cens[miss]
      Y_curr[i, miss] <- sapply(cens_i, fix_start_value)
    }
  }

  # First imputation (s = 1)
  Y[, , 1] <- t(sapply(1:n, function(i) {
    tnorm_impute_state(
      obs_row  = Z[i, ],
      curr_row = Y_curr[i, ],
      mu       = mu,
      Sigma    = Sigma,
      cens_row = if (cens_is_matrix) cens[i, ] else cens,
      ag       = ag,
      burn.in.samples = burn.in.samples,
      thinning = thinning
    )
  }))

  # --- Gibbs iterations ---------------------------------------------------
  for (s in 2:R) {
    # 1. Update mu
    ybar <- colMeans(Y[, , s - 1])
    A_n  <- solve(Omega) + n * solve(Sigma)
    b_n  <- solve(Omega) %*% mu_0 + n * solve(Sigma) %*% ybar
    mu   <- as.numeric(rmvnorm(1, solve(A_n) %*% b_n, solve(A_n)))

    # 2. Update Sigma
    Ym   <- sweep(Y[, , s - 1], 2, mu, "-")
    S_mu <- t(Ym) %*% Ym
    S_n  <- solve(S_0 + S_mu)
    Sigma <- solve(rwish(1, v0 + n, S_n)[, , 1])

    # 3. Imputation step: warm-start from Y[, , s-1]
    Y[, , s] <- t(sapply(1:n, function(i) {
      tnorm_impute_state(
        obs_row  = Z[i, ],
        curr_row = Y[i, , s - 1],
        mu       = mu,
        Sigma    = Sigma,
        cens_row = if (cens_is_matrix) cens[i, ] else cens,
        ag       = ag,
        burn.in.samples = burn.in.samples,
        thinning = thinning
      )
    }))

    if (s %% 100 == 0) cat("MVN iteration", s, "/", R, "\n")
  }

  list(Y = Y, mu = mu, Sigma = Sigma)
}

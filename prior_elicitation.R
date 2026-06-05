# =========================================================================
# Prior Elicitation 
# =========================================================================
#
# Standard (weakly informative) prior specification for the MVN, MFM and DPML
# imputation models, as described in Section 4.2 of the paper.
# For the MFM imputation model, the same weakly informative prior specification
# used for the MVN model is adopted.

#' Standard prior elicitation for the LOD imputation models
#'
#' @param model  Character: "MVN", "MFM" or "DPML"
#' @param Z      Data matrix (n x q) with NA for censored entries
#' @param Ek     Expected number of clusters (DPML only, default 1)
#' @param scale  Logical: if TRUE (default), assumes standardised data
#'               (mu_0 = 0, Omega = diag(2.5))
#' @return Named list of prior hyperparameters
standard_prior_elicitation <- function(model, Z, Ek = 1, scale = TRUE) {
  
  # --- DP concentration parameter ----------------------------------------
  c_value <- function(Ek, N) {
    objective <- function(c, Ek, N) {
      abs(sum(c / (1:N - 1 + c)) - Ek)
    }
    optim(par = 1, fn = objective, method = "L-BFGS-B",
          lower = 1e-6, upper = 1e9, N = N, Ek = Ek)$par
  }
  
  q  <- ncol(Z)
  
  # Inverse-Wishart prior for Sigma
  v0 <- q + 4
  S0 <- (v0 - q - 1) * cov(Z, use = "complete.obs")
  
  # Mean prior
  if (scale) {
    m0    <- rep(0, q)
    Omega <- diag(2.5, q)
  } else {
    m0    <- colMeans(Z, na.rm = TRUE)
    Omega <- diag(apply(Z, 2, var, na.rm = TRUE), q)
  }
  
  # Output
  if (model == "MVN" | model == "MFM") {
    list(m0 = m0, Omega = Omega, v0 = v0, S0 = S0)
  } else if (model == "DPML") {
    list(mu_G0 = m0, Sigma_G0 = Omega,
         c = c_value(Ek, nrow(Z)), v0 = v0, S0 = S0)
  } else {
    stop("model must be 'MVN', 'MFM' or 'DPML'")
  }
}

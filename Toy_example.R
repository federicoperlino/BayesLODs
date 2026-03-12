######### Toy Example ######################################################
#
# Minimal working example demonstrating the corrected MVN and DPML
# imputation models with the scatterplot matrix visualisation.

set.seed(1234)

# Source all modules
source("utils.R")
source("gibbs_mvn.R")
source("gibbs_dpml.R")
source("prior_elicitation.R")
source("diagnostics.R")

# -------------------------------------------------------------------------
# 1. Generate data
# -------------------------------------------------------------------------

n <- 200; q <- 8

# AR(1) covariance to induce realistic dependence across variables.
idx <- seq_len(q)
Sigma <- 0.7 ^ abs(outer(idx, idx, "-"))
Y_true <- mvtnorm::rmvnorm(n = n, mean = rep(1, q), sigma = Sigma)

# Standard LODs (common across individuals)
LODs <- apply(Y_true, 2, quantile, probs = 0.25)

# Censor the data
Z_LODs <- sapply(1:q, function(j)
  ifelse(Y_true[, j] <= LODs[j], NA, Y_true[, j]))

# Remove rows that are entirely censored under the common-LOD setting.
keep_std <- !apply(Z_LODs, 1, function(r) all(is.na(r)))
Y_true <- Y_true[keep_std, , drop = FALSE]
Z_LODs <- Z_LODs[keep_std, , drop = FALSE]
n <- nrow(Y_true)

cat("Percentage censored per variable:\n")
print(round(100 * colMeans(is.na(Z_LODs)), 1))

# Individual-specific LODs
IS_LODs <- t(sapply(1:n, function(x)
  apply(Y_true, 2, quantile, probs = 0.25))) + rnorm(n)

Z_IS_LODs_full <- ifelse(Y_true <= IS_LODs, NA, Y_true)

# Remove rows that are entirely NA
keep       <- !apply(Z_IS_LODs_full, 1, function(r) all(is.na(r)))
Z_IS_LODs  <- Z_IS_LODs_full[keep, ]
IS_LODs    <- IS_LODs[keep, ]

# -------------------------------------------------------------------------
# 2. Run MVN model (standard LODs)
# -------------------------------------------------------------------------

R <- 500  # increase for production use

priors_mvn <- standard_prior_elicitation("MVN", Z = Z_LODs)

fit_mvn <- gibbs_mvn(
  Z     = Z_LODs,
  mu_0  = priors_mvn$mu,
  Omega = priors_mvn$Omega,
  v0    = priors_mvn$v0,
  S_0   = priors_mvn$S0,
  cens  = LODs,
  R     = R
)

# -------------------------------------------------------------------------
# 3. Run DPML model (standard LODs)
# -------------------------------------------------------------------------

priors_dpml <- standard_prior_elicitation("DPML", Z = Z_LODs, Ek = 2)

fit_dpml <- gibbs_DPML(
  Z       = Z_LODs,
  mu_0    = priors_dpml$mu_G0,
  Sigma_0 = priors_dpml$Sigma_G0,
  c       = priors_dpml$c,
  v0      = priors_dpml$v0,
  S_0     = priors_dpml$S0,
  cens    = LODs,
  R       = R
)

# -------------------------------------------------------------------------
# 4. Basic checks
# -------------------------------------------------------------------------

cat("\n--- MVN check ---\n")
check_mvn <- check_imputation_LOD(
  Z_orig = Z_LODs,
  Z_imp  = fit_mvn$Y[, , sample((R %/% 2):R, 1)],
  LOD    = LODs
)

cat("\n--- DPML check ---\n")
check_dpml <- check_imputation_LOD(
  Z_orig = Z_LODs,
  Z_imp  = fit_dpml$Y[, , sample((R %/% 2):R, 1)],
  LOD    = LODs
)

# -------------------------------------------------------------------------
# 5. Visualisation: scatterplot matrix with uncertainty encoding
# -------------------------------------------------------------------------

# Pick 4 variables for a readable plot. Change `sel` to use any subset
# and any ordering of the variables.
sel <- 1:4
varnames <- paste0("X", sel)
burn <- R %/% 5
draw_show <- burn + ceiling((R - burn) / 2)

mvn_modes <- c(alpha = "alpha", size = "size",
               hybrid = "both", whisker = "whisker")

for (nm in names(mvn_modes)) {
  pdf(paste0("scatterplot_mvn_", nm, ".pdf"), width = 8, height = 8)
  scatterplot_matrix_ci(
    Y           = fit_mvn$Y[, sel, ],
    Z           = Z_LODs[, sel],
    LOD         = LODs[sel],
    alpha       = 0.05,
    burn        = burn,
    varnames    = varnames,
    position    = "sample",
    draw        = draw_show,
    uncertainty = mvn_modes[[nm]],
    main        = paste("MVN Model -", toupper(nm), "view")
  )
  dev.off()
}

pdf("scatterplot_mvn.pdf", width = 8, height = 8)
scatterplot_matrix_ci(
  Y           = fit_mvn$Y[, sel, ],
  Z           = Z_LODs[, sel],
  LOD         = LODs[sel],
  alpha       = 0.05,
  burn        = burn,
  varnames    = varnames,
  position    = "sample",
  draw        = draw_show,
  uncertainty = "both",
  main        = "MVN Model - scatterplot matrix"
)
dev.off()

pdf("scatterplot_dpml.pdf", width = 8, height = 8)
scatterplot_matrix_ci(
  Y           = fit_dpml$Y[, sel, ],
  Z           = Z_LODs[, sel],
  LOD         = LODs[sel],
  alpha       = 0.05,
  burn        = burn,
  varnames    = varnames,
  position    = "sample",
  draw        = draw_show,
  uncertainty = "both",
  main        = "DPML Model - scatterplot matrix"
)
dev.off()

cat(
  "\nPlots saved: scatterplot_mvn.pdf, scatterplot_dpml.pdf,\n",
  "scatterplot_mvn_alpha.pdf, scatterplot_mvn_size.pdf,\n",
  "scatterplot_mvn_hybrid.pdf, scatterplot_mvn_whisker.pdf\n",
  sep = ""
)

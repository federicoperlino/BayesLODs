# =========================================================================
# Toy example
# =========================================================================
#
# Minimal working example demonstrating the updated MVN, MFM, DPML, and
# fixed-K finite mixture imputation models, including scatterplot matrix
# visualisation.
#
# Section 8 reports imputation error diagnostics, which require access to the
# true latent values and are therefore specific to the toy example setting.

rm(list=ls())

# Source all modules
source("utils.R")
source("gibbs_mvn.R")
source("gibbs_mfm.R")
source("gibbs_dpml.R")
source("prior_elicitation.R")
source("diagnostics.R")
source("dgp.R")

# -------------------------------------------------------------------------
# 0. Global settings
# -------------------------------------------------------------------------

seed_data <- 1     # seed for data generation
seed_mcmc <- 123   # seed for MCMC initialization
R <- 500           # increase for production use

# -------------------------------------------------------------------------
# 1. Generate data 
# -------------------------------------------------------------------------

n <- 300; q <- 7

# AR(1) covariance to induce realistic dependence across variables.
idx <- seq_len(q)
Sigma <- 0.75 ^ abs(outer(idx, idx, "-"))

# Censoring applied to the latent covariates
censoring <- c(0.20, 0.20, 0.20, 0.15, 0.15, 0.10, 0.05) 

# Choose one latent data-generating distribution:
#   1.1 MVN
#   1.2 Mixture
#
# IMPORTANT:
# Run either Section 1.1 or Section 1.2, not both.
# If both are run sequentially, Section 1.2 overwrites the objects
# created in Section 1.1.

# -------------------------------------------------------------------------
# 1.1 Generate data from a MVN distribution
# -------------------------------------------------------------------------

# Mean 
mu_raw <- rep(0, q)

# Simulate data from a MVN model
dat <- simulate_lod_dgp(
  n = n,
  y_distribution = "mvn",
  mu_raw = mu_raw,
  Sigma_Y = Sigma,
  censoring = censoring, 
  seed = seed_data
)

# Standard LODs (common across individuals)
LODs <- dat$LOD_std

# Censor the data
Z_LODs <- dat$Z_out

cat("Percentage censored per variable:\n")
print(round(100 * colMeans(is.na(Z_LODs)), 1))

# -------------------------------------------------------------------------
# 1.2 Generate data from a mixture distribution
# -------------------------------------------------------------------------

# Mixture specification
pi_mix <- c(0.35, 0.65)
mu_components_raw <- rbind(
  comp1 = c(-1.10, -0.95,  0.75, -0.55, -0.35, -0.55, -0.65),
  comp2 = c( 1.00,  0.85, -0.75,  0.65,  0.60,  0.70,  0.80)
)

# Simulate data from a mixture model
dat <- simulate_lod_dgp(
  n = n,
  y_distribution = "mixture",
  pi_mix = pi_mix,
  mu_components_raw = mu_components_raw,
  Sigma_Y = Sigma,
  censoring = censoring, 
  seed = seed_data
)

# Standard LODs (common across individuals)
LODs <- dat$LOD_std

# Censor the data
Z_LODs <- dat$Z_out

cat("Percentage censored per variable:\n")
print(round(100 * colMeans(is.na(Z_LODs)), 1))

# -------------------------------------------------------------------------
# 2. Run MVN model (standard LODs)
# -------------------------------------------------------------------------

set.seed(seed_mcmc)
priors_mvn <- standard_prior_elicitation("MVN", Z = Z_LODs)

start_time <- Sys.time()

fit_mvn <- gibbs_mvn(
  Z     = Z_LODs,
  mu_0  = priors_mvn$m0,
  Omega = priors_mvn$Omega,
  v0    = priors_mvn$v0,
  S_0   = priors_mvn$S0,
  cens  = LODs,
  R     = R
)

runtime_mvn <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

# -------------------------------------------------------------------------
# 3. Run DPML model (standard LODs)
# -------------------------------------------------------------------------

set.seed(seed_mcmc)
priors_dpml <- standard_prior_elicitation("DPML", Z = Z_LODs, Ek = 2)

start_time <- Sys.time()

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

runtime_dpml <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

# -------------------------------------------------------------------------
# 4. Run MFM model (standard LODs)
# -------------------------------------------------------------------------

set.seed(seed_mcmc)
priors_mfm <- standard_prior_elicitation("MFM", Z = Z_LODs)

start_time <- Sys.time()

fit_mfm <- gibbs_MFM(
  Z     = Z_LODs,
  m0    = priors_mfm$m0,
  Omega = priors_mfm$Omega,
  v0    = priors_mfm$v0,
  S0    = priors_mfm$S0,
  cens  = LODs,
  R     = R, 
  c_total = 1  # Smaller values encourage sparse weights. 
               # Larger values favor more balanced weights and 
               # may support more occupied mixture components.
)

runtime_mfm <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

burn <- R %/% 5
post_idx <- (burn + 1):R
table(fit_mfm$K_plus[post_idx]) 
# If K_plus collapses to 1 after burn-in, inspect mixing and consider
# a sensitivity run with a larger c_total, e.g. 2.

# -------------------------------------------------------------------------
# 5. Run Fixed-K Finite Mixture model (standard LODs)
# -------------------------------------------------------------------------

set.seed(seed_mcmc)
priors_mfm <- standard_prior_elicitation("MFM", Z = Z_LODs)

start_time <- Sys.time()

fit_fixed <- gibbs_MFM(
  Z     = Z_LODs,
  m0    = priors_mfm$m0,
  Omega = priors_mfm$Omega,
  v0    = priors_mfm$v0,
  S0    = priors_mfm$S0,
  cens  = LODs,
  R     = R, 
  K_plus_fixed = 2,
  c_total = 1  # Smaller values encourage sparse weights. 
               # Larger values favor more balanced weights.
)

runtime_fixed <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

# -------------------------------------------------------------------------
# 6. Basic checks
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

cat("\n--- MFM check ---\n")
check_mfm  <- check_imputation_LOD(
  Z_orig = Z_LODs,
  Z_imp  = fit_mfm$Y[, , sample((R %/% 2):R, 1)],
  LOD    = LODs
)

cat("\n--- Fixed-K Finite Mixture check ---\n")
check_fixed <- check_imputation_LOD(
  Z_orig = Z_LODs,
  Z_imp  = fit_fixed$Y[, , sample((R %/% 2):R, 1)],
  LOD    = LODs
)

# -------------------------------------------------------------------------
# 7. Visualisation: scatterplot matrix with uncertainty encoding
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


pdf("scatterplot_mfm.pdf", width = 8, height = 8)
scatterplot_matrix_ci(
  Y           = fit_mfm$Y[, sel, ],
  Z           = Z_LODs[, sel],
  LOD         = LODs[sel],
  alpha       = 0.05,
  burn        = burn,
  varnames    = varnames,
  position    = "sample",
  draw        = draw_show,
  uncertainty = "both",
  main        = "MFM Model - scatterplot matrix"
)
dev.off()


pdf("scatterplot_fixed.pdf", width = 8, height = 8)
scatterplot_matrix_ci(
  Y           = fit_fixed$Y[, sel, ],
  Z           = Z_LODs[, sel],
  LOD         = LODs[sel],
  alpha       = 0.05,
  burn        = burn,
  varnames    = varnames,
  position    = "sample",
  draw        = draw_show,
  uncertainty = "both",
  main        = "Fixed-K Finite Mixture Model - scatterplot matrix"
)
dev.off()

cat(
  "\nPlots saved: scatterplot_dpml.pdf, scatterplot_mfm.pdf, scatterplot_fixed.pdf \n",
  "scatterplot_mvn_alpha.pdf, scatterplot_mvn_size.pdf,\n",
  "scatterplot_mvn_hybrid.pdf, scatterplot_mvn_whisker.pdf\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 8. Evaluation: imputation accuracy and posterior uncertainty
# -------------------------------------------------------------------------

# This section is possible only because the toy example stores the true
# uncensored values in dat$Y_true. It is not applicable to real-world data.

Z_obs    <- dat$Z_out
Y_true   <- dat$Y_true
idx_cens <- is.na(Z_obs)

# Posterior summary
post_mvn <- posterior_summary(
  Y     = fit_mvn$Y,
  Z     = Z_obs,
  alpha = 0.05,
  burn  = burn
)

post_dpml <- posterior_summary(
  Y     = fit_dpml$Y,
  Z     = Z_obs,
  alpha = 0.05,
  burn  = burn
)

post_mfm <- posterior_summary(
  Y     = fit_mfm$Y,
  Z     = Z_obs,
  alpha = 0.05,
  burn  = burn
)

post_fixed <- posterior_summary(
  Y     = fit_fixed$Y,
  Z     = Z_obs,
  alpha = 0.05,
  burn  = burn
)

# Helper function for summarising accuracy and posterior uncertainty
# on censored cells only.
summarise_imputation <- function(post, method, Y_true, idx_cens) {
  
  Y_hat <- post$median
  err   <- Y_hat[idx_cens] - Y_true[idx_cens]
  
  data.frame(
    Method = method,
    N_censored = sum(idx_cens),
    MAE = mean(abs(err)),
    RMSE = sqrt(mean(err^2)),
    Coverage = mean(
      Y_true[idx_cens] >= post$lo[idx_cens] &
        Y_true[idx_cens] <= post$hi[idx_cens]
    ),
    Width = mean(post$hi[idx_cens] - post$lo[idx_cens]),
    SD = mean(post$sd[idx_cens])
  )
}

# Summary table
tab <- rbind(
  summarise_imputation(post_mvn,    "MVN",      Y_true, idx_cens),
  summarise_imputation(post_dpml,   "DPML",     Y_true, idx_cens),
  summarise_imputation(post_mfm,    "MFM",      Y_true, idx_cens), 
  summarise_imputation(post_fixed,  "Fixed_K",  Y_true, idx_cens)
)

# Rounding
tab$MAE      <- round(tab$MAE, 4)
tab$RMSE     <- round(tab$RMSE, 4)
tab$Coverage <- round(tab$Coverage, 3)
tab$Width    <- round(tab$Width, 4)
tab$SD       <- round(tab$SD, 4)

# Run time per model in seconds
runtime_sec <- c(
  MVN     = runtime_mvn,
  DPML    = runtime_dpml,
  MFM     = runtime_mfm,
  Fixed_K = runtime_fixed
)

tab$Runtime_sec <- round(runtime_sec[tab$Method], 2)

# Print tab
# Results may vary across simulated datasets and MCMC runs.
# The table is intended as a single-run diagnostic, not as a formal benchmark.
tab


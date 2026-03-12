# Bayesian Imputation for Multivariate Left-Censored Data

R implementation accompanying the paper:

> **A Bayesian Parametric and Nonparametric Approach for the Imputation of Multivariate Left-Censored Data due to Limit of Detection**  
> Federico L. Perlino, Bernardo Nipoti, Paige L. Williams, and Andrea Bellavia.  
> *Statistics in Medicine*, 2025. [DOI: 10.1002/sim.70326](https://doi.org/10.1002/sim.70326)

## Structure

```
utils.R               Shared helpers: samplers, warm-started TMV kernels
gibbs_mvn.R           Parametric (MVN) imputation model
gibbs_dpml.R          Nonparametric (DPML) imputation model
prior_elicitation.R   Standard prior specification (Section 4.2)
diagnostics.R         LOD checks, posterior summaries,
                      scatterplot matrix with uncertainty encoding
Toy_example.R         Minimal working example
```

## What changed

- The codebase has been refactored and organised so that the MVN and
  DPML imputation models are kept in separate scripts, while shared
  sampling helpers live in `utils.R`.
- Diagnostics and visualisation utilities are collected in
  `diagnostics.R`, including the scatterplot matrix with
  pointwise uncertainty encoding.
- The example workflow in `Toy_example.R` has been refreshed and now uses
  a correlated AR(1) covariance structure, so the scatterplots are more
  informative than a purely independent toy dataset.
- Tuning controls for the truncated sampling step are exposed directly in
  the main Gibbs samplers.

## Quick start

```r
source("utils.R")
source("gibbs_mvn.R")
source("gibbs_dpml.R")
source("prior_elicitation.R")
source("diagnostics.R")

priors <- standard_prior_elicitation("MVN", Z = Z_LODs)

fit <- gibbs_mvn(
  Z = Z_LODs, mu_0 = priors$mu, Omega = priors$Omega,
  v0 = priors$v0, S_0 = priors$S0, cens = LODs, R = 1000
)

scatterplot_matrix_ci(
  Y           = fit$Y,
  Z           = Z_LODs,
  LOD         = LODs,
  burn        = 200,
  position    = "sample",
  uncertainty = "both",
  main        = "MVN scatterplot matrix"
)
```

Alternative uncertainty encodings:

- `uncertainty = "alpha"`: transparency.
- `uncertainty = "size"`: point size.
- `uncertainty = "both"`: transparency and point size.
- `uncertainty = "whisker"`: pointwise credible intervals on imputed axes.

Variable selection is controlled by subsetting the columns passed to the
plotting function. For example:

```r
sel <- c(4, 1, 3, 2)

scatterplot_matrix_ci(
  Y        = fit$Y[, sel, ],
  Z        = Z_LODs[, sel],
  LOD      = LODs[sel],
  varnames = paste0("X", sel)
)
```

`Toy_example.R` writes example PDF files to the working directory. These
outputs are ignored by `.gitignore` and are not intended to be tracked.

## Dependencies

- `tmvtnorm` — truncated multivariate normal sampling
- `plyr` — frequency counting in the Pólya urn (DPML model)
- `mvtnorm` — multivariate normal density evaluation (DPML model)

## Citation

```bibtex
@article{perlino2025bayesian,
  title   = {A Bayesian Parametric and Nonparametric Approach for the
             Imputation of Multivariate Left-Censored Data Due to
             Limit of Detection},
  author  = {Perlino, Federico L. and Nipoti, Bernardo and Williams,
             Paige L. and Bellavia, Andrea},
  journal = {Statistics in Medicine},
  volume  = {44},
  pages   = {e70326},
  year    = {2025},
  doi     = {10.1002/sim.70326}
}
```

## Licence

This repository is released under the MIT Licence. See `LICENSE`.

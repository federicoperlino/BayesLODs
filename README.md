# Bayesian Imputation for Multivariate Left-Censored Data Suffering from Limit of Detection

R implementation accompanying the paper:

> **A Bayesian Parametric and Nonparametric Approach for the Imputation of Multivariate Left-Censored Data due to Limit of Detection**  
> Federico L. Perlino, Bernardo Nipoti, Paige L. Williams, and Andrea Bellavia.  
> *Statistics in Medicine*, 2025. [DOI: 10.1002/sim.70326](https://doi.org/10.1002/sim.70326)

## Structure

```
utils.R               Shared helpers: numerical stabilizers, samplers,
                      and warm started TMV kernels
gibbs_mvn.R           Parametric (MVN) imputation model
gibbs_dpml.R          Nonparametric (DPML) imputation model
gibbs_mfm.R           Dynamic mixture of finite mixtures (MFM) and 
                      optional fixed-K finite mixture imputation models
prior_elicitation.R   Standard prior specification (Section 4.2)
diagnostics.R         LOD checks, posterior summaries,
                      and scatterplot matrix with uncertainty encoding
dgp.R                 Data-generating process for the toy example scenarios
Toy_example.R         Minimal working example
```

## What changed

- Numerical stabilizers have been added in `utils.R` to improve robustness. 
  These include matrix symmetrisation, positive-definiteness checks, 
  and small diagonal jitter corrections before matrix inversions and
  Cholesky decompositions.
- The imputation kernels have been updated to handle fully censored rows 
  by sampling from the corresponding truncated multivariate marginal distribution.
- The MVN and DPML Gibbs samplers have been adapted to use the new numerical
  stabilizers and the updated imputation kernels defined in `utils.R`.
- The MFM sampler has been added in `gibbs_mfm.R` as a dynamic mixture of 
  finite mixtures, with posterior updating of the number of mixture components. 
  An optional fixed-K finite mixture mode has also been added. When this option 
  is used, the number of occupied mixture components is specified by the 
  practitioner and the sampler should be interpreted as a finite mixture sampler 
  rather than as a dynamic MFM.
- The prior specification function in `prior_elicitation.R` has been updated to 
  use a common notation across models and to include the prior quantities 
  required by the MFM sampler.
- The data-generating process has been moved to `dgp.R`, which contains the
  functions used to simulate latent data under different scenarios.
- The example workflow in `Toy_example.R` now includes an additional evaluation 
  section comparing the models in terms of imputation accuracy, 
  posterior uncertainty, and runtime. This section relies on the true latent 
  values and is therefore specific to the toy example setting.

## Quick start

```r
source("utils.R")
source("gibbs_mvn.R")
source("gibbs_mfm.R")
source("gibbs_dpml.R")
source("prior_elicitation.R")
source("diagnostics.R")

priors <- standard_prior_elicitation("MVN", Z = Z_LODs)

fit <- gibbs_mvn(
  Z = Z_LODs, mu_0 = priors$m0, Omega = priors$Omega,
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

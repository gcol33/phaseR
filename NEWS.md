# phaseR 1.0.0

First stable release.

## New Features (since 0.5.0)

* `fixef()`: Extract fixed effect estimates
* `residuals()`: Compute response and Pearson residuals
* `logLik()`, `nobs()`: Standard model accessors for AIC/BIC
* `confint()`: Credible intervals for parameters
* `waic()`, `compare_waic()`: Model comparison via WAIC
* `prior()`: Prior specification interface (defaults documented)
* `sim_phaseR_binomial()`: Simulate binomial phase data
* Improved input validation with informative error messages

## Bug Fixes

* Fixed namespace imports for S3 method generics

---

# phaseR 0.5.0

Initial release.

## Features

* Phase-structured models for longitudinal data
* k-phase support (arbitrary number of sequential phases)
* GLM families: Gaussian, Poisson, Binomial
* Random intercepts via `(1 | group)` syntax
* Counterfactual inference: `ate()`, `ite()`, `counterfactual()`
* Native HMC sampler (no external dependencies)
* Comprehensive diagnostics: Rhat, ESS, trace plots

## Core Functions

* `phases()`: Declare ordered phases
* `transition()`: Specify transition models
* `dynamics()`: Specify phase-specific outcome models
* `phase_model()`: Compose complete model
* `fit_phaseR()`: Fit via Hamiltonian Monte Carlo
* `posterior_predict()`: Posterior predictive simulation

## Vignettes

* "Why Conditioning on Realized Change Fails" - demonstrates the estimand problem

# phaseR 2.0.0

Initial CRAN release.

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

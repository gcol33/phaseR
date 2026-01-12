# phaseR 1.0.4

## New Features

* **Laplace approximation backend**: Fast approximate inference for large datasets
  - `backend = "laplace"` option in `fit_phaseR()`
  - Finds posterior mode via L-BFGS-B optimization
  - Standard errors from inverse Hessian
  - Samples drawn from Gaussian approximation
  - Works with all model types: basic, RE, crossed RE, slopes, nested RE

## Internal Changes

* New `R/backend_laplace.R` with full implementation
* `build_neg_log_posterior()` dispatches to correct C++ likelihood
* `sample_mvnorm()` for drawing from multivariate normal
* `MASS_ginv()` for pseudo-inverse without MASS dependency

---

# phaseR 1.0.3

## New Features

* **Nested random effects**: Support for `(1|site/plot)` syntax
  - Automatically expands to `(1|site) + (1|site:plot)`
  - Three-level nesting: `(1|region/site/plot)`
  - Creates interaction variables automatically
  - `validate_nesting()` helper for checking proper nesting structure

## Internal Changes

* New `expand_nested_re()` function for parsing nested syntax
* New `create_nested_interactions()` for creating interaction columns
* Extended regex patterns to handle group names with colons

---

# phaseR 1.0.2

## New Features

* **Random slopes**: Support for `(1 + x | group)` and `(1 + x || group)` syntax
  - Correlated random intercepts and slopes via Cholesky parameterization
  - Uncorrelated random slopes via `||` syntax with diagonal covariance
  - Slope-only RE via `(0 + x | group)` or `(x | group)`
  - Multiple slope variables: `(1 + x + z | group)`
  - LKJ prior on correlation matrix for correlated slopes

* **Enhanced VarCorr()**: Returns full covariance matrix for correlated slopes
  - Standard deviations and correlations from posterior
  - Pretty-print method for variance components

## Internal Changes

* New `parse_re_term()` function for parsing slope terms
* New `build_re_z_matrix()` for constructing RE design matrices
* Added `phase_likelihood_slopes.cpp` and `hmc_sampler_slopes.cpp`
* Extended parameter naming for slope coefficients and Cholesky factors

---

# phaseR 1.0.1

## New Features

* **Multiple crossed random effects**: Support for `(1|site) + (1|year)` syntax
  - Multiple RE terms per model component (transition, dynamics)
  - New C++ sampler for multi-RE models
  - Automatic dispatch between single-RE and multi-RE backends

## Internal Changes

* Extended formula parser for multiple RE terms
* New `build_multi_re()` function for RE structure construction
* Added `phase_likelihood_multi_re.cpp` and `hmc_sampler_multi_re.cpp`

---

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

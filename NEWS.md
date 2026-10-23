# phaseR (development version)

* `prior` now applies to every model family under both backends. The HMC
  backend, and the Laplace backend for GLM, k-phase, crossed random-effect and
  random-slope models, evaluated hard-coded priors and ignored `prior`. Both
  backends now evaluate one log posterior per model (`posterior_target()`),
  whose prior comes from `prior()`. The hard-coded priors differed from the
  documented defaults, so default fits of those models change: residual SDs
  get Half-Normal(0, 1) on the SD instead of Normal(0, 1) on its logarithm,
  random-effect SDs Half-Normal(0, 0.5) instead of Normal(0, 0.5) on the
  logarithm, and correlated random slopes LKJ(2) instead of LKJ(1). Under
  Poisson and binomial families the unused residual SDs now get the sigma
  prior, which makes their posterior proper.

* Each phase of a two-phase model now uses the design matrix of its own
  `dynamics()` formula in the likelihood. The C++ likelihoods computed the
  second phase's mean from the first phase's design matrix, so fits whose two
  `dynamics()` formulas differed estimated the wrong model, and read past the
  end of the matrix when the second formula had more columns. Predictions,
  counterfactuals and WAIC already used the per-phase matrices.

* The HMC backend is the No-U-Turn sampler (Hoffman and Gelman 2014) in
  Stan's multinomial form, replacing fixed-length leapfrog trajectories:
  trajectories double until a U-turn, a divergence or `max_treedepth`, and
  warmup adapts the step size by dual averaging and a diagonal metric from
  the posterior variance in expanding windows. `adapt_delta` and
  `max_treedepth` are passed through `...` of `fit_phaseR()`. One sampler
  now serves every model family. Fits report divergences after warmup,
  transitions at maximum tree depth, E-BFMI and the adapted metric.

* `summary()` computes the effective sample size with Stan's multi-chain
  initial monotone sequence estimator.

* `fit_phaseR()` honours `seed` with the HMC backend: each chain samples from a
  stream seeded by `(seed, chain)` that also draws its initial values, so two
  fits with the same seed return identical draws. With `seed = NULL` the seed is
  drawn from R's random number generator, so `set.seed()` makes fits
  reproducible, and the seed used is stored in the fit as `fit$seed`.

# phaseR 0.8.0

Experimental release with feature parity to ratiod.

## Features

* **Phase-structured models** for longitudinal data
  - k-phase support (arbitrary number of sequential phases)
  - GLM families: Gaussian, Poisson, Binomial
  - Counterfactual inference: `ate()`, `ite()`, `counterfactual()`

* **Random effects**
  - Random intercepts: `(1 | group)`
  - Multiple crossed RE: `(1|site) + (1|year)`
  - Random slopes: `(1 + x | group)` and `(1 + x || group)`
  - Nested RE: `(1|site/plot)`

* **Inference backends**
  - Native HMC/NUTS sampler (no external dependencies)
  - Laplace approximation: `backend = "laplace"`

* **Custom priors**
  - Formula syntax: `prior(beta ~ normal(0, 2.5), sigma ~ half_cauchy(1))`
  - Distributions: `normal()`, `student_t()`, `cauchy()`, `half_normal()`, `half_cauchy()`, `lkj()`

* **Model comparison**
  - WAIC: `waic()`, `compare_waic()`
  - LOO-CV via PSIS: `loo()`, `loo_compare()`

## Core Functions

* `phases()`: Declare ordered phases
* `transition()`: Specify transition models
* `dynamics()`: Specify phase-specific outcome models
* `phase_model()`: Compose complete model
* `fit_phaseR()`: Fit model
* `posterior_predict()`: Posterior predictive simulation
* `fixef()`, `ranef()`, `VarCorr()`: Extract estimates
* `residuals()`, `confint()`: Diagnostics

---

# phaseR 0.5.0

Initial release.

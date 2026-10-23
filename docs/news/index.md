# Changelog

## phaseR (development version)

- `prior` now applies to every model family under both backends. The HMC
  backend, and the Laplace backend for GLM, k-phase, crossed
  random-effect and random-slope models, evaluated hard-coded priors and
  ignored `prior`. Both backends now evaluate one log posterior per
  model
  ([`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)),
  whose prior comes from
  [`prior()`](https://gillescolling.com/phaseR/reference/prior.md). The
  hard-coded priors differed from the documented defaults, so default
  fits of those models change: residual SDs get Half-Normal(0, 1) on the
  SD instead of Normal(0, 1) on its logarithm, random-effect SDs
  Half-Normal(0, 0.5) instead of Normal(0, 0.5) on the logarithm, and
  correlated random slopes LKJ(2) instead of LKJ(1). Under Poisson and
  binomial families the unused residual SDs now get the sigma prior,
  which makes their posterior proper.

- Each phase of a two-phase model now uses the design matrix of its own
  [`dynamics()`](https://gillescolling.com/phaseR/reference/dynamics.md)
  formula in the likelihood. The C++ likelihoods computed the second
  phase’s mean from the first phase’s design matrix, so fits whose two
  [`dynamics()`](https://gillescolling.com/phaseR/reference/dynamics.md)
  formulas differed estimated the wrong model, and read past the end of
  the matrix when the second formula had more columns. Predictions,
  counterfactuals and WAIC already used the per-phase matrices.

- The HMC backend is the No-U-Turn sampler (Hoffman and Gelman 2014) in
  Stan’s multinomial form, replacing fixed-length leapfrog trajectories:
  trajectories double until a U-turn, a divergence or `max_treedepth`,
  and warmup adapts the step size by dual averaging and a diagonal
  metric from the posterior variance in expanding windows. `adapt_delta`
  and `max_treedepth` are passed through `...` of
  [`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md).
  One sampler now serves every model family. Fits report divergences
  after warmup, transitions at maximum tree depth, E-BFMI and the
  adapted metric.

- [`summary()`](https://rdrr.io/r/base/summary.html) computes the
  effective sample size with Stan’s multi-chain initial monotone
  sequence estimator.

- [`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md)
  honours `seed` with the HMC backend: each chain samples from a stream
  seeded by `(seed, chain)` that also draws its initial values, so two
  fits with the same seed return identical draws. With `seed = NULL` the
  seed is drawn from R’s random number generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes fits
  reproducible, and the seed used is stored in the fit as `fit$seed`.

## phaseR 0.8.0

Experimental release with feature parity to ratiod.

### Features

- **Phase-structured models** for longitudinal data
  - k-phase support (arbitrary number of sequential phases)
  - GLM families: Gaussian, Poisson, Binomial
  - Counterfactual inference:
    [`ate()`](https://gillescolling.com/phaseR/reference/ate.md),
    [`ite()`](https://gillescolling.com/phaseR/reference/ite.md),
    [`counterfactual()`](https://gillescolling.com/phaseR/reference/counterfactual.md)
- **Random effects**
  - Random intercepts: `(1 | group)`
  - Multiple crossed RE: `(1|site) + (1|year)`
  - Random slopes: `(1 + x | group)` and `(1 + x || group)`
  - Nested RE: `(1|site/plot)`
- **Inference backends**
  - Native HMC/NUTS sampler (no external dependencies)
  - Laplace approximation: `backend = "laplace"`
- **Custom priors**
  - Formula syntax:
    `prior(beta ~ normal(0, 2.5), sigma ~ half_cauchy(1))`
  - Distributions:
    [`normal()`](https://gillescolling.com/phaseR/reference/normal.md),
    [`student_t()`](https://gillescolling.com/phaseR/reference/student_t.md),
    [`cauchy()`](https://gillescolling.com/phaseR/reference/cauchy.md),
    [`half_normal()`](https://gillescolling.com/phaseR/reference/half_normal.md),
    [`half_cauchy()`](https://gillescolling.com/phaseR/reference/half_cauchy.md),
    [`lkj()`](https://gillescolling.com/phaseR/reference/lkj.md)
- **Model comparison**
  - WAIC:
    [`waic()`](https://gillescolling.com/phaseR/reference/waic.md),
    [`compare_waic()`](https://gillescolling.com/phaseR/reference/compare_waic.md)
  - LOO-CV via PSIS:
    [`loo()`](https://gillescolling.com/phaseR/reference/loo.md),
    [`loo_compare()`](https://gillescolling.com/phaseR/reference/loo_compare.md)

### Core Functions

- [`phases()`](https://gillescolling.com/phaseR/reference/phases.md):
  Declare ordered phases
- [`transition()`](https://gillescolling.com/phaseR/reference/transition.md):
  Specify transition models
- [`dynamics()`](https://gillescolling.com/phaseR/reference/dynamics.md):
  Specify phase-specific outcome models
- [`phase_model()`](https://gillescolling.com/phaseR/reference/phase_model.md):
  Compose complete model
- [`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md):
  Fit model
- [`posterior_predict()`](https://gillescolling.com/phaseR/reference/posterior_predict.md):
  Posterior predictive simulation
- [`fixef()`](https://gillescolling.com/phaseR/reference/fixef.md),
  [`ranef()`](https://gillescolling.com/phaseR/reference/ranef.md),
  [`VarCorr()`](https://gillescolling.com/phaseR/reference/VarCorr.md):
  Extract estimates
- [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html): Diagnostics

------------------------------------------------------------------------

## phaseR 0.5.0

Initial release.

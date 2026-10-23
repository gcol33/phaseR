# phaseR (development version)

* The Laplace backend searches the posterior mode from data-informed
  starting points instead of all-zero parameters, and NUTS chains start
  from them. A phase model's posterior has a mode in which no unit enters a
  phase and that phase's parameters sit at their prior; from zeros, far
  from the data, the optimiser stopped there (60 units x 4 times from
  `sim_phaseR(seed = 8)`: active intercept -0.003 with SD 2.5 under
  Laplace, 8.15 with SD 0.23 under NUTS), and so could a NUTS chain. Each
  starting point fits every phase's regression (Gaussian, Poisson or
  logistic, ridge-stabilised) and every transition's logistic regression to
  a provisional phase assignment: one splits each unit into equal runs, the
  others draw each unit's moves at random from `seed`. Laplace optimises
  from `n_starts` of them (default 5) and keeps the highest log posterior,
  reported in `diagnostics$start_values`; NUTS chain `c` starts at the
  `c`-th. Both are reproducible from `seed`, and R's random number stream
  outside the fit is left untouched.

* Both backends use the exact gradient of the log posterior instead of
  finite differences. The backward recursion over the latent phases gives
  the derivative of each unit's log likelihood with respect to its emission
  and transition terms (the posterior probabilities of each phase and each
  move), and the chain rule carries it through the outcome families, the
  linear predictors, the random-effect parameterisations and the priors.
  Tests compare it with Richardson-extrapolated central differences on
  sixteen models (relative error below 1e-6). On the README example
  (100 units x 5 times, 2 chains x 1000 iterations) NUTS takes 6.5 s
  instead of 67.5 s at the same minimum ESS (721 against 734), and the
  Laplace optimiser gets the gradient through `optim(gr = )`.

* Random-effect terms whose grouping factor varies within a unit, such as
  `(1 | year)` crossed with units, now add the effect of each observation's
  own level. The likelihood used the level of the unit's first observation
  at every time, so a crossed year effect acted as an effect of the first
  year. Random slopes without an intercept, `(0 + x | group)`, now multiply
  their covariate; the likelihood added them as random intercepts.

* All model families evaluate their likelihood through one forward recursion
  over the latent phases (`src/phase_process.h`) and one model description
  (`src/phase_model.h`): K ordered phases, one logistic transition predictor
  per consecutive pair, one dynamics predictor per phase with a Gaussian,
  Poisson or binomial outcome, and random-effect terms on any predictor. The
  five likelihood files that each carried their own copy of the recursion
  are gone; log posteriors of every model kind agree with the previous code
  to 1e-10 at fixed parameters (`tests/testthat/fixtures`).

* Uncorrelated random slopes, `(1 + x || group)`, now give the effects of
  coefficient j the prior Normal(0, sigma_j). The prior read the SD of each
  effect from a position that moved with the effect, past the SD block and
  past the end of the parameter vector.

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

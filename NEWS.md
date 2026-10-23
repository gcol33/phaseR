# phaseR (development version)

* The Laplace backend is faster for models with random effects: on 80
  units with 12 sites it takes 0.37 s instead of 1.61 s for site
  intercepts in both phases and 0.45 s instead of 4.23 s for a correlated
  site slope, with the same mode (to 0.002 posterior SD). The gradient of
  the log marginal now costs one third-derivative pass per colour of the
  effects and one Hessian-vector product, instead of one per colour and
  parameter; the Hessian in the effects is block diagonal over groups of
  effects that share no unit and is factorised block by block, with
  Newton's method run in each block separately; and each starting point is
  first moved to the mode of the fixed effects with the effects at zero,
  so starts that arrive at the same point share one search. The mode of
  the effects is still found from `u = 0` at every parameter value, so the
  approximation remains a function of the parameters.

* A Laplace search that reaches parameters where the mode of the random
  effects is degenerate or not finite stops that start instead of looping
  forever, and the fit reports when this happened from every start or
  when the effects of some groups have several modes at the optimum
  (`diagnostics$multimodal_groups`). `diagnostics$cov` is the covariance of
  the parameters outside the random effects; `diagnostics$se` covers all.

* Correlated random effects, `(1 + x | group)`, are parameterised by the
  log SD of each coefficient and the Cholesky factor of their correlation
  matrix, built from canonical partial correlations (the transform of
  Lewandowski, Kurowicka and Joe 2009 that Stan uses). The coefficients are
  `diag(sigma) L z` with standardised effects `z`; the SDs get the
  random-effect prior (`re ~`) and the correlation matrix exactly the
  `cor ~ lkj(eta)` prior, Jacobians included. The previous parameterisation,
  a free lower-triangular factor of the covariance, combined an LKJ-type
  term on its diagonal with the random-effect prior reflected onto it; that
  was not the LKJ density of any correlation matrix, and the factor's sign
  made its posterior bimodal. Parameters are named
  `log_sigma_<component>_re_<group>_<coef>` and
  `atanh_cpc_<component>_re_<group>_<i>_<j>`.

* `ranef()` and `VarCorr()` read the random effects from the model
  description. `ranef()` returns the coefficients of every level
  (`diag(sigma) L z` for correlated terms, which it returned as the
  standardised `z`) as one data frame with a row per predictor, term,
  level and coefficient, and found the levels of dynamics terms by phase
  name (it looked them up under `phase_0`/`phase_1` and fell back to
  `group_1`, ...). `VarCorr()` returns for every term the SD of each
  coefficient and the posterior mean correlation matrix.

* Models with more than two phases take outcome families and random
  effects like two-phase models. `prepare_stan_data()` built them with
  Gaussian dynamics and the fixed part of each formula, dropping
  `family` and every random-effect term without a message. One data
  preparation now serves every number of phases; it matches each
  `dynamics()` and `transition()` to its phases by name (a two-phase model
  read its dynamics in the order given), uses an intercept-only Gaussian
  model for a phase without `dynamics()`, and rejects a transition that
  skips a phase, which the left-to-right phase process cannot represent.

* Predictions, counterfactuals and treatment effects are computed by the
  model both backends fit (`src/phase_model.h`) instead of R copies of the
  two-phase Gaussian recursion, so they cover every model: k phases,
  Poisson and binomial dynamics and random effects (levels the fit has not
  seen contribute no effect). For two Gaussian phases without random
  effects, `posterior_predict(type = "phase")` and
  `predict(type = "transition")` are unchanged. Other changes:
  - `posterior_predict()` returns rows in the order of `newdata`; it
    returned them sorted by unit and time, so `ite()` matched effects to
    the wrong units when the data were not sorted. `type = "phase"` gains
    `phase`, the phase whose probability is returned (default: the last).
  - `predict(type = "phase")` and `predict(type = "response")` average over
    the posterior draws; they used the transition coefficients at their
    posterior mean, and `"response"` ignored the outcome family (link) and
    the random effects. `predict(type = "transition")` returns one element
    per move for models with more than two phases.
  - `ate()` and `ite()` contrast the expected outcomes under treatment and
    control, integrating over the phases, for the same posterior draws.
    They contrasted two independent outcome simulations, from separately
    sampled draws when `n_sims` was set, which added outcome noise to the
    effect's interval.
  - `fitted()` returns the one-step predictive mean
    E(y_t | y_0, ..., y_(t-1)) of each observation, `residuals()` the
    corresponding residuals (Pearson residuals divide by the predictive
    SD), and `logLik()` the log-likelihood marginalised over the phases at
    the posterior mean. They used a Gaussian density around phase-averaged
    predictions with the average residual SD.

* The Cholesky factor of correlated random effects with three or more
  coefficients is named in the order the likelihood and `VarCorr()` read it
  (by column: `L_..._1_1`, `L_..._2_1`, `L_..._3_1`, `L_..._2_2`, ...). The
  names followed the rows, so `L_..._2_2` labelled the draws of `L_31`.

* `waic()` and `loo()` compute the pointwise log likelihood with the forward
  recursion the backends fit with. The R copy of the recursion they used
  knew only two Gaussian phases without random effects, so for models with
  random effects, Poisson or binomial outcomes or more than two phases it
  evaluated a different model (random effects at zero, every outcome
  Gaussian, two phases).

* The Laplace backend integrates random effects out before it looks for a
  mode. It approximated the joint posterior of all parameters by a Gaussian
  at its mode, but with random effects that mode does not exist: the SDs
  are already parameterised by their logarithm with the Jacobian included,
  yet as a random-effect SD shrinks with its effects at zero, the Normal
  density of the G effects grows like SD^-G while a proper SD prior times
  the Jacobian falls only like SD. The optimiser ran towards zero SD and
  the Hessian there was singular ("Hessian not invertible, using
  pseudo-inverse", a ridge and an eigenvalue shift then made it
  invertible). The marginal posterior of the other parameters, with the
  effects integrated out by the Laplace formula (Tierney and Kadane 1986,
  as in INLA and TMB: the mode of the effects given the parameters, and the
  log determinant of the negative Hessian in the effects there), is bounded
  and has a proper mode. The backend maximises it with its exact gradient
  (third derivatives by forward-mode dual numbers through the exact
  gradient of the log posterior), takes the covariance from the Hessian of
  the log marginal, and draws the effects given the parameters from their
  Gaussian, linearised in the parameters. The ridge, the pseudo-inverse and
  the eigenvalue shift are gone; a Hessian that is not positive definite
  now stops the fit with a message. On a random-intercept model (60 units x
  5 times, 10 sites) Laplace means agree with NUTS within 0.5 posterior SD
  (0.75 for the skewed log SDs of the random effects) and its SDs within
  a factor of 1.25 (1.5). Models without random effects get the same
  Gaussian at the joint mode as before.

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
